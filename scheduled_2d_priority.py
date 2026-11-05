"""One-time 05:20 packaging of existing 2D primary models, with no new solves."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import argparse,csv,json,shutil,zipfile,traceback

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'output/COMSOL_reproduction_20261014'
DEST=ROOT/'output/COMSOL_2D_priority_20261014_0520'
ARCHIVE=ROOT/'output/二维优先复核资料_20261014_0520.zip'
STATE=ROOT/'output/二维定时整理执行状态.json'
BJT=timezone(timedelta(hours=8))
SCHEDULED=datetime(2026,10,14,5,20,tzinfo=BJT)


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def state(**values):
    STATE.write_text(json.dumps(dict(time_beijing=datetime.now(BJT).isoformat(),**values),ensure_ascii=False,indent=2),encoding='utf-8')


def validate():
    rows=read(SOURCE/'export_main.json')
    assert len(rows)==84
    for row in rows:
        path=SOURCE/row['destination']
        assert path.is_file() and path.stat().st_size==row['exported_bytes']
        assert row['category']=='main_models' and row['reopened_solution_empty']
        dimensions=[d for c in row.get('structure',{}).values() for d in c.get('geometry_dimensions',{}).values()]
        # The first exported smoke case predates the optional structure record;
        # its source signature explicitly identifies the same 2D plate layout.
        assert dimensions==[2] or (not dimensions and len(row['original_signature']['layout'])==3)
    assert len(read(SOURCE/'export_history.json'))==124
    assert (SOURCE/'薄板孔隙群数值实验流程与复核手册.docx').is_file()
    return rows


def run():
    if datetime.now(BJT)<SCHEDULED:raise RuntimeError('Not yet 05:20 Beijing time; no packaging performed.')
    rows=validate();state(status='running',scheduled_beijing=SCHEDULED.isoformat())
    DEST.mkdir(parents=True,exist_ok=True)
    shutil.copytree(SOURCE/'main_models',DEST/'main_models',dirs_exist_ok=True)
    shutil.copytree(SOURCE/'参考数据',DEST/'参考数据',dirs_exist_ok=True)
    for name in ['薄板孔隙群数值实验流程与复核手册.docx','主模型参数索引.csv','export_main.json','recompute_smoke_check.json']:
        shutil.copy2(SOURCE/name,DEST/name)
    with (SOURCE/'模型清单.csv').open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f);fields=reader.fieldnames;records=[r for r in reader if r['编号'].startswith('M')]
    with (DEST/'二维模型清单.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(records)
    (DEST/'先读我_二维优先.md').write_text('''# 二维优先复核资料

本包只包含已经检查过的84个二维主模型，不重新整理历史版本，不启动新仿真。
历史122个文件及2个已知问题旧版本已经在原完整包中，不需要重复整理；原文件未覆盖。
三维模型本次不复制、不求解，留作后续优化和探索。

1. 先阅读Word手册第5节，从main_models/reference/model.mph开始。
2. mode_id为0时S0入射，为1时A0入射；freq_id按各文件的说明选择，f0由公式更新。
3. 先复核无孔与单孔，再复核孔群、长度对照及统计衰减。
4. 参数、求解器与原网格已保留，MPH不含求解场；参考数据为核对用的小型数值表。

Word沿用已逐页检查的完整手册。其中第10节的historical_models和historical_invalid
指原完整包中的历史目录，并非本二维包缺失文件。请以本包的二维模型清单.csv为准。

后续二维验证重点：随机样本收敛、100 kHz A0网格检查、吸收区长度与轮廓分辨率检查。
当前参考衰减系数尚未全部通过稳定性判据，不能作为已收敛材料常数。
''',encoding='utf-8')
    (DEST/'三维后续探索计划.md').write_text('''# 三维后续探索计划

当前状态为待开展，不表示已完成新的三维验证。

- 首先以二维结果完成网格、随机样本和边界误差排查，再选少量有代表性的参数进入三维。
- 使用三维随机场生成独立三维孔隙，保留三维尺寸、形貌与分布统计，不能将二维孔形直接拉伸替代。
- 补充实测三维材料或铺层参数，再检查导波模态、方向性及三维散射与二维趋势的差异。
- 继续区分透射、反射与模态转换的检测判据，并给出数值误差和随机样本不确定性。

现有三维历史文件仅保留在原完整包中，本任务不重复整理或运行。
''',encoding='utf-8')
    assert len(list(DEST.rglob('*.mph')))==84
    temporary=ARCHIVE.with_suffix('.zip.partial')
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1,allowZip64=True) as z:
        for path in sorted(DEST.rglob('*')):
            if path.is_file():z.write(path,Path(DEST.name)/path.relative_to(DEST))
    with zipfile.ZipFile(temporary) as z:
        assert sum(n.endswith('.mph') for n in z.namelist())==84
        assert z.testzip() is None
    temporary.replace(ARCHIVE)
    state(status='complete',models=84,archive=str(ARCHIVE),archive_bytes=ARCHIVE.stat().st_size,
          original_package_unchanged=True,history_reprocessed=False,new_comsol_solves=0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--validate-only',action='store_true');a=p.parse_args()
    if a.validate_only:
        print(json.dumps(dict(validated_models=len(validate()),scheduled_beijing=SCHEDULED.isoformat(),packaging_started=False),ensure_ascii=False))
    else:
        try:run()
        except Exception as e:
            state(status='failed',error=str(e),traceback=traceback.format_exc());raise
