"""Assemble manifests, numerical benchmark tables, and a portable ZIP."""
from pathlib import Path
import json,csv,shutil,zipfile,hashlib

ROOT=Path(__file__).resolve().parent;DEST=ROOT/'output/COMSOL_reproduction_20261014';DATA=ROOT/'output/guided_collective'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write(p,value):Path(p).write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')
def csvwrite(path,rows):
    with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

items=[]
for category in ['main','history']:
    rows=read(DEST/f'export_{category}.json')
    assert len(rows)==(84 if category=='main' else 124)
    for row in rows:
        source=Path(row['source']);target=DEST/row['destination']
        assert source.stat().st_mtime_ns==row['source_mtime_ns'],'Original model modified'
        assert target.exists() and row['reopened_solution_empty'] and all(row['cleared_solutions'].values())
        if category=='main':
            record=read(source.parent/'results.json')
            row.update(original_signature=record['signature'],material_status=record['material_status'],finite_voids=record.get('finite_voids'))
            if (source.parent/'geometry/metadata.json').exists():row['pore_metadata']=read(source.parent/'geometry/metadata.json')
        else:
            # COMSOL's meshfilenamehidden is a UI boolean, not an external file.
            row['external_inputs']=[a for a in row.get('external_inputs',[]) if not a['property'].lower().endswith('hidden')]
            for asset in row['external_inputs']:
                assert asset['exists'] and (target.parent/asset['relative_path']).is_file(),asset
        write(target.with_suffix('.parameters.json'),row);items.append(row)
    write(DEST/f'export_{category}.json',rows)
assert len(list(DEST.rglob('*.mph')))==208

manifest=[];mode_rows=[];mi=hi=0
for row in items:
    main=row['category']=='main_models'
    if main:mi+=1;identifier=f'M{mi:03}'
    else:hi+=1;identifier=f'H{hi:03}'
    status='本轮主模型' if main else ('已知孔隙合并旧版 仅追溯' if row['category']=='historical_invalid' else '历史版本 非本轮推荐')
    if not main and not row['cleared_solutions']:status+=' 无求解器序列 可能仅几何草稿'
    manifest.append(dict(编号=identifier,分类=status,模型相对路径=row['destination'],原模型相对路径=row['relative_source'],
                         清除场解后重新打开检查='通过',模型大小MB=round(row['exported_bytes']/1e6,3),
                         频率kHz='/'.join(f'{f/1000:g}' for f in row.get('frequencies_hz',[])),
                         模态选择='0=S0 1=A0' if main else '保留原设置'))
    if not main:continue
    record=read(Path(row['source']).parent/'results.json');name=row['destination']
    L=48 if 'L48' in name else (24 if 'L24' in name else 8)
    if 'checks/' in name:ref_name='checks/'+('reference_fine' if L==8 else f'reference_L{L}_fine')
    else:ref_name='reference' if L==8 else f'reference_L{L}'
    reference=read(DATA/ref_name/'results.json')
    for result in record['results']:
        ref=next(r for r in reference['results'] if r['frequency_hz']==result['frequency_hz'])
        for i,mode in enumerate(['S0','A0']):
            t=complex(result['transmission_real'][i][i],result['transmission_imag'][i][i])
            tr=complex(ref['transmission_real'][i][i],ref['transmission_imag'][i][i]);tn=t/tr
            rmax=max(abs(complex(result['reflection_real'][j][i],result['reflection_imag'][j][i])-complex(ref['reflection_real'][j][i],ref['reflection_imag'][j][i])) for j in range(2))
            mode_rows.append(dict(模型编号=identifier,模型相对路径=name,参考模型=ref_name+'/model.mph',频率Hz=result['frequency_hz'],入射模态=mode,
                                 T同模态实部=t.real,T同模态虚部=t.imag,归一化T实部=tn.real,归一化T虚部=tn.imag,
                                 个体幅值下降百分比=100*(1-abs(tn)),单孔式绝对幅值变化百分比=100*abs(1-abs(tn)),反射变化最大百分比=100*rmax))
csvwrite(DEST/'模型清单.csv',manifest)
assert len(mode_rows)==394
refdir=DEST/'参考数据';refdir.mkdir(exist_ok=True)
csvwrite(refdir/'394个模态工况对照数据.csv',mode_rows)
for name in ['sweep_summary.json','extended_attenuation.json','length_attenuation.json','numerical_checks.json',
             'single_controls_summary.json','threshold_sensitivity.json','study_settings.json','basis_convergence.json']:
    shutil.copy2(DATA/name,refdir/name)
shutil.copy2(ROOT/'config.json',refdir/'示例材料与生成配置.json')
(refdir/'说明.txt').write_text('本目录为复算时比较用的小型数值表，不是COMSOL场解。394行覆盖本轮84个主模型的已计算频率及两种入射。模型清单中的编号与Word附录一致。个体幅值下降可能为负，不予截断。统计集合均值必须对复数T先平均，不能直接平均此表的幅值下降代替。\n',encoding='utf-8')
(DEST/'先读我.txt').write_text('''本包采用 COMSOL 6.2 保存，共208个MPH文件，均已清除求解数据并重新打开核对。

建议先阅读《薄板孔隙群数值实验流程与复核手册.docx》第5节。
优先使用 main_models 中的84个本轮主模型。
mode_id：0=S0，1=A0。
freq_id：按每个模型参数说明选择有效编号；f0自动更新，不要覆盖其公式。
三频模型为1=50、2=100、3=200 kHz；双频为1=100、2=200 kHz；单频为1=200 kHz。
内嵌模态剖面仅对应本轮示例材料和厚度。改变材料或厚度后应重算模态剖面。

模型保留几何、网格、物理场、材料、边界条件、研究与求解器设置。
场解、数值求值缓存和结果表已清空；原始含场解文件没有被覆盖。
参考数据目录是供核对用的数值表，不包含场解向量。

historical_models 中122个文件保留早期设置，部分属于几何草稿或尚未充分验证的版本。
historical_invalid 中2个文件有已知孔隙合并问题，仅作历史追溯，不用于本轮统计。
历史模型不自动获得频率和模态选择功能。详见模型清单.csv。
历史三维模型的assets文件夹应随模型一起复制；重新导入时如提示路径问题，请在导入节点选择该文件夹中的对应文件。

所有文件为另外保存的副本。复算建议先另存工作副本，再点击计算。
两种模式的抽样复算与原接收系数相对差小于1e-10；该检查只确认打包前后的一致性，不表示科学验证的全部条件均已通过。
''',encoding='utf-8')
summary=dict(models=208,main_models=84,historical_models=122,known_invalid_historical_models=2,
             all_reopened_without_solution=True,original_files_unchanged_by_mtime=True,benchmark_rows=394,
             model_bytes=sum((DEST/r['destination']).stat().st_size for r in items),
             smoke_check=read(DEST/'recompute_smoke_check.json'))
write(DEST/'交付检查.json',summary)
archive=ROOT/'output/COMSOL无结果模型与实验手册_20261014.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1,allowZip64=True) as z:
    for file in sorted(DEST.rglob('*')):
        if file.is_file() and not file.name.endswith(('.lock','.mphbin','.recovery')):
            z.write(file,Path(DEST.name)/file.relative_to(DEST))
with zipfile.ZipFile(archive) as z:
    assert sum(n.endswith('.mph') for n in z.namelist())==208
    assert z.testzip() is None
print(json.dumps(dict(summary,zip_path=str(archive),zip_bytes=archive.stat().st_size),ensure_ascii=False,indent=2),flush=True)
