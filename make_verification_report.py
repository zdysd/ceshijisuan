"""Summarize fresh COMSOL verification; preserve unsuccessful checks explicitly."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
plt.rcParams['axes.unicode_minus']=False

ROOT=Path(__file__).resolve().parent


def load(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    out=ROOT/'output/revalidation'
    audit=load(out/'verification.json')
    transmission=load(out/'incident_normalized_transmission.json')['results']
    fine=load(out/'transmission_mesh_check.json')
    corrected=ROOT/'output/revalidation_corrected_domain'
    gates=load(corrected/'gates.json')
    with (corrected/'validation.csv').open(encoding='utf-8-sig') as stream:
        domain=list(csv.DictReader(stream))
    cfg=load(ROOT/'config.json')
    corrected_forward=dict(gates,frequency_hz=5e6,
                           alpha_range_np_m=[min(float(r['alpha_np_m']) for r in domain),max(float(r['alpha_np_m']) for r in domain)],
                           warning='Absolute cumulative tolerance does not establish a stable material attenuation coefficient.')
    (corrected/'forward_gates.json').write_text(json.dumps(corrected_forward,indent=2),encoding='utf-8')
    legacy_geometry=[]
    for row in audit['coefficient_convergence']:
        record=load(ROOT/'output/final_validation'/row['case']/'pores/measurement.json')
        legacy_geometry.append(dict(case=row['case'],requested_porosity=record['porosity'],
                                    actual_porosity=record['actual_vector_porosity'],
                                    passed=abs(record['actual_vector_porosity']-record['porosity'])<=cfg['geometry_porosity_absolute_tolerance']))
    # Retire the incomplete legacy gate without changing its simulation data.
    old_gate_path=ROOT/'output/final_validation/forward_gates.json'
    old=load(old_gate_path)
    old.setdefault('numerical_pass_before_geometry_audit',old['numerical_pass'])
    old['gates']['geometry_porosity']=dict(passed=all(x['passed'] for x in legacy_geometry))
    old['numerical_pass']=all(x['passed'] for x in old['gates'].values())
    old['relative_domain_coefficient_passed']=False
    old['superseded_by']='output/revalidation_corrected_domain/gates.json'
    old_gate_path.write_text(json.dumps(old,indent=2),encoding='utf-8')
    aggregate=dict(reproduction=audit,incident_normalized_transmission=transmission,
                   transmission_mesh_check=fine,legacy_geometry_check=legacy_geometry,
                   corrected_domain_gates=gates,
                   collective_response_supported_in_local_periodic_example=all(r['group_change_exceeds_largest_single'] for r in transmission),
                   relative_transmission_mesh_checks=[dict(frequency_hz=r['frequency_hz'],case=r['case'],passed=r['relative_reduction_change']<=.05) for r in fine],
                   experimental_undetectability_verified=False,centimetre_extrapolation_validated=False)
    (out/'validated_summary.json').write_text(json.dumps(aggregate,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    x=np.arange(len(transmission)); width=.24
    for offset,key,label in [(-1,'median_single_amplitude_reduction','中位尺寸单孔'),
                             (0,'largest_single_amplitude_reduction','最大尺寸单孔'),
                             (1,'coherent_group_amplitude_reduction','孔群相干平均响应')]:
        values=np.array([r[key] for r in transmission])*100
        axes[0].bar(x+offset*width,values,width,label=label)
        for xi,yi in zip(x+offset*width,values):axes[0].annotate(f'{yi:.3g}',(xi,yi),xytext=(0,4),textcoords='offset points',ha='center',fontsize=8)
    for i,r in enumerate(transmission):
        mean=r['coherent_group_amplitude_reduction']*100
        lo,hi=np.array(r['amplitude_reduction_ci95'])*100
        axes[0].errorbar(i+width,mean,yerr=[[mean-lo],[hi-mean]],color='black',capsize=3)
    axes[0].set(xticks=x,xticklabels=['5 MHz','10 MHz'],ylabel='入射波归一化后的幅值下降（%）',title='粗网格对照：孔群95%自助重采样置信区间')
    axes[0].legend(fontsize=8);axes[0].grid(axis='y',alpha=.2)
    values=np.array([float(r['alpha_np_m']) for r in domain])*20/np.log(10)/1000
    axes[1].bar(np.arange(len(domain)),values,color=['#2674a6','#2674a6','#2674a6','#b95034','#b95034'])
    axes[1].set(xticks=np.arange(len(domain)),xticklabels=['基准模型','加密网格','完整轮廓','长度加倍','宽度加倍'],ylabel='局部拟合衰减系数（dB/mm）',title='几何修正后：仍需验证域尺寸收敛')
    axes[1].tick_params(axis='x',labelsize=8);axes[1].grid(axis='y',alpha=.2)
    fig.savefig(out/'verification_comparison.png',dpi=180);plt.close(fig)
    comparison_rows=[]
    for r in transmission:
        ci=np.array(r['amplitude_reduction_ci95'])*100
        comparison_rows.append(f"| {r['frequency_hz']/1e6:g} MHz | {100*r['median_single_amplitude_reduction']:.4f}% | {100*r['largest_single_amplitude_reduction']:.4f}% | {100*r['group_seed101_amplitude_reduction']:.3f}% | {100*r['coherent_group_amplitude_reduction']:.3f}% [{ci[0]:.3f}, {ci[1]:.3f}] |")
    mesh_rows=[]
    for r in fine:
        mesh_rows.append(f"| {r['frequency_hz']/1e6:g} MHz / {r['case']} | {100*r['coarse_amplitude_reduction']:.5f}% | {100*r['fine_amplitude_reduction']:.5f}% | {100*r['absolute_reduction_change']:.5f}个百分点 | {100*r['relative_reduction_change']:.2f}% |")
    domain_rows=[]
    base=float(domain[0]['alpha_np_m'])
    for r in domain:
        alpha=float(r['alpha_np_m']);change=abs(alpha-base)/max(abs(base),abs(alpha))
        domain_rows.append(f"| {r['case']} | {100*float(r['actual_porosity']):.6f}% | {alpha:.6f} | {100*change:.2f}% |")
    text=f'''# 数值复核与单孔／孔群对照

**局部二维示例支持通过整体透射响应表征孔隙群；衰减系数的域尺寸收敛尚不能据此宣告完成。** 用户尚无实验噪声／分辨率判据，故本次不判定单孔在实验中不可检测。材料卡含示例补全参数，周期横向边界也是模型假设。

## 本次实际验证

- 重新打开20个随机孔隙体波模型及两个无孔参考，从位移与法向应力重新积分提取前行波；5、10 MHz的相干衰减仍为0.078716、0.512426 dB/mm，逐样本系数、集合系数及500次Bootstrap区间与原记录的差值均为0。
- 重新打开三维模型，确认1426440个自由度、299751个四面体，重新提取的探针相对差值为0；一维等效模型的保存场差值也为0。此项验证保存解的完整性，不替代三维网格／随机集合收敛。
- 科学计算14项检查通过，增加了“近零累计误差通过但系数相对精度未通过”的回归断言。
- 新求解中位尺寸和最大尺寸单孔，沿用seed 101孔群中的原轮廓和位置；新求解加密网格以及修正孔隙率容差后的域尺寸对照。

## 对照量与结果

在孔隙段前x=0.15 mm处提取入射前行波，在孔隙段后x=1.45、1.55、1.65 mm处接收；这些截面均无孔。定义

`T = (A_forward,out / A_forward,in)_pores / (A_forward,out / A_forward,in)_reference`

表中为`100*(1-|T|)`，三个出口截面取平均；孔群集合先平均复数T。它是相干位移幅值变化，不是全部能量损失。若只用固定载荷下的出口幅值比较，激励端反射会改变实际入射波，甚至造成相反趋势；本表使用入射波归一化后的结果。

| 频率 | 中位尺寸单孔 | 最大单孔 | 同一132孔样本 | 10样本相干孔群均值及95%区间 |
|---|---:|---:|---:|---:|
{chr(10).join(comparison_rows)}

seed 101的中位、最大等效孔径分别为9.688、33.809 μm；示例介质体波波长为480.54 μm（5 MHz）和240.27 μm（10 MHz）。相关长度20、10 μm并不是孔的半轴。单孔与孔群的总孔隙率不同，因此这个对照证明整体响应的累积变化，不能单独证明聚集效应、非线性相互作用或多重散射占主导。

![复核结果]({(out/'verification_comparison.png').as_posix()})

## 透射量的网格检查

| 频率／案例 | 原网格幅值下降 | 加密网格幅值下降 | 绝对变化 | 相对变化 |
|---|---:|---:|---:|---:|
{chr(10).join(mesh_rows)}

很小的单孔响应必须连同网格敏感性报告；不能把一个很小的正数直接当成已经辨识的可检测信号。此处没有用绝对0.01 Np门限替单孔响应保证5%相对精度。

**中位尺寸单孔在两个频率下的幅值变化均随加密翻号，因此粗网格约0.01%的变化不能作为可靠的物理衰减值。** 微小负值保留用于数值误差诊断，不截断成零。5 MHz最大单孔的相对变化约12.7%，也尚未达到5%目标；10 MHz最大单孔约0.32%。5 MHz同一132孔样本的整体幅值下降由0.854%变为0.843%，相对变化约1.28%，通过此项检查。10 MHz孔群尚无本次独立加密对照，不能把5 MHz网格结论自动用于10 MHz。

## 修正后的孔隙率与域尺寸复算

旧0.2 μm轮廓简化的base、mesh_fine、length_long、width_wide模型实际孔隙率偏离名义1.5%超过设定的0.01个百分点容差。原验收遗漏了这项，现已加入检查，旧`forward_gates.json`已标记不通过并指向新复算。原数值解仍保留用于追溯。

| 案例 | 实际孔隙率 | 局部α（Np/m） | 相对新base的系数变化 |
|---|---:|---:|---:|
{chr(10).join(domain_rows)}

新几何孔隙率容差通过：**{gates['gates']['geometry_porosity']['passed']}**；近零累计误差等数值门限通过：**{gates['numerical_pass']}**；长度／宽度的5%相对系数条件通过：**{gates['relative_domain_coefficient_passed']}**。改变域尺寸也改变有限随机样本，后续应对每种域尺寸使用匹配孔隙统计的独立集合，不能把单次变化完全归因于边界误差。

两个频率的Bootstrap区间相对半宽约24%、27%；前5个与10个样本的点估计接近不等于95%区间已缩窄到5%。目前的等效模型只用于原1 mm局部示例，厘米级传播与A0/S0衰减外推均未获验证。

## 复现

依次使用ASCII路径COMSOL运行环境执行`verify_outputs.py`、`verify_transmission.py`、`run_study.py validate --output output/revalidation_corrected_domain`、`validate_transmission_mesh.py`；用项目`.venv`运行`make_verification_report.py`生成本报告。首次核验会重新求解单孔模型；完整域复算需要额外时间。

机器可读证据：`validated_summary.json`；原始复信号：`reextracted_signals.npz`、`incident_normalized_transmission.npz`；新求解的单孔`.mph`位于本目录，域复算模型位于`../revalidation_corrected_domain/`。
'''
    (out/'REPORT.md').write_text(text,encoding='utf-8')


if __name__=='__main__':main()
