"""Audit actual pore-free references and display the normalized baseline explicitly."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from guided_statistics import normalized_transmission,coherent_interval

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'output/guided_collective'
OUT=ROOT/'output/baseline_comparison';OUT.mkdir(exist_ok=True)
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
settings=read(DATA/'study_settings.json')
items=read(ROOT/'output/COMSOL_reproduction_20261014/export_main.json')
audit=[]
for item in items:
    name=item['destination'].removeprefix('main_models/').removesuffix('/model.mph')
    length=48 if 'L48' in name else (24 if 'L24' in name else 8)
    ref=('checks/'+('reference_fine' if length==8 else f'reference_L{length}_fine')) if name.startswith('checks/') else ('reference' if length==8 else f'reference_L{length}')
    current=read(DATA/name/'results.json');base=read(DATA/ref/'results.json')
    a=current['signature'];b=base['signature']
    assert b['geometry_hash']==hashlib.sha256(b'').hexdigest()
    assert a['material']==b['material'] and a['layout']==b['layout'] and a['mesh']==b['mesh']
    assert a.get('port_positions_m',settings['pre_positions']+settings['post_positions'])==b.get('port_positions_m',settings['pre_positions']+settings['post_positions'])
    assert set(a['frequencies']).issubset(b['frequencies'])
    for f in a['frequencies']:
        assert np.max(abs(normalized_transmission(base,base,f)-1))<1e-12
    audit.append(dict(case=name,reference=ref,baseline_has_no_pores=True,
                      matched_material_layout_receivers_background_mesh=True,frequencies_hz=a['frequencies']))

plot=[]
for length in [8,24,48]:
    ref=DATA/('reference' if length==8 else f'reference_L{length}')
    values=[]
    for seed in range(101,113):
        source=DATA if seed<=106 else DATA/'extended_samples'
        values.append(normalized_transmission(read(source/f'd80/L{length}/seed_{seed}/results.json'),read(ref/'results.json'),200000.))
    values=np.array(values)
    for index,mode in enumerate(['S0','A0']):
        result=coherent_interval(values[:,index])
        ci=result['ci95_amplitude_reduction']
        plot.append(dict(mode=mode,length_mm=length,reference=ref.relative_to(DATA).as_posix(),
                         baseline_normalized_percent=100.,porous_coherent_transmission_percent=100*(1-result['amplitude_reduction']),
                         ci95_percent=[100*(1-ci[1]),100*(1-ci[0])]))

payload=dict(normalization='t_m = T_m(porous) / T_m(pore-free)',threshold_fraction=.001,
             threshold_role='Decision threshold after normalization; not a reference model and not porosity.',
             checked_cases=len(audit),case_reference_mapping=audit,
             fit_metadata_note='reference_subtracted in the generic fit records only in-function division; these guided fits receive transmissions already normalized upstream.',
             demonstration=dict(frequency_hz=200000,diameter_group_um=80,porosity=.015,samples=12,rows=plot))
(OUT/'基板比较核对.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':11})
fig,axes=plt.subplots(1,2,figsize=(12,5.3),sharey=True)
for ax,mode,color in zip(axes,['S0','A0'],['#2563a5','#16806a']):
    rows=[r for r in plot if r['mode']==mode]
    x=np.array([r['length_mm'] for r in rows]);y=np.array([r['porous_coherent_transmission_percent'] for r in rows]);ci=np.array([r['ci95_percent'] for r in rows])
    ax.plot(x,[100]*3,'s-',color='#343e48',lw=2,label='无孔基板：归一化基准 100%')
    ax.axhline(99.9,color='#bf6324',ls='--',lw=1.4,label='下降 0.1% 的示例判据线')
    ax.errorbar(x,y,yerr=[y-ci[:,0],ci[:,1]-y],fmt='o-',color=color,lw=2,capsize=4,label='含孔板：相干集合幅值与 95% 区间')
    ax.set(title=f'{mode} 模态',xlabel='孔隙段长度（毫米）',xticks=x,ylim=(99.54,100.055),xlim=(4,52))
    ax.grid(alpha=.2);ax.legend(loc='lower left',fontsize=9)
axes[0].set_ylabel('相对于无孔基板的透射幅值（%）')
fig.suptitle('实际比较基准是无孔基板，0.1% 是判断变化的阈值',fontsize=16,y=.99)
fig.text(.5,.915,'已有数据重绘：200 千赫 · 80 微米孔径组 · 孔隙段面积率 1.5% · 12 个随机样本',ha='center',fontsize=11)
fig.text(.07,.025,'无孔基板的 100% 来自同频同模态的参考归一化，不表示参考计算没有数值误差。',fontsize=10,color='#54616d')
fig.subplots_adjust(top=.81,bottom=.16,left=.08,right=.98,wspace=.12)
fig.savefig(OUT/'无孔基板与孔群透射对比.png',dpi=180,bbox_inches='tight')
plt.close(fig)
(OUT/'基板比较说明.md').write_text('''# 当前比较基准核对

本轮84个主模型的匹配关系已核对。比较基准为相同材料、外形、频率、模态、接收位置及背景网格级别的无孔基板。
归一化采用 t = T含孔 / T无孔；T为模态分离并按实际入射归一化后的同模态透射系数。
无孔基板自身归一化为1，即100%；幅值下降图中的基板因此对应0%。
0.1%仅是归一化后的示例检测阈值，和基板定义、孔隙率不是同一个量。

8毫米孔隙段用reference；24毫米用reference_L24；48毫米用reference_L48。
这些无孔基板的总尺寸均约338毫米×1毫米，名称中的L24和L48用于表示匹配的接收截面位置，不表示板长。
加密网格模型对应checks中的reference_fine、reference_L24_fine、reference_L48_fine。

原始拟合JSON中的reference_subtracted字段仅表示通用拟合函数内部是否再次执行参考除法。
本轮导波结果在调用拟合函数之前已经完成无孔参考归一化，不能用该字段为false推断没有使用基板。
数值结果和模型均未因本次说明而重新求解或更改。
''',encoding='utf-8')
print(json.dumps(dict(audited_cases=len(audit),demonstration=plot),ensure_ascii=False,indent=2))
