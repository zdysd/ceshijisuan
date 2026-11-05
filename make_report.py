"""Render a report from measured outputs; never fill missing runs with estimates."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
plt.rcParams['axes.unicode_minus']=False


def load(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    root=Path('output');fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    summaries=[]
    for folder,label,color in [('final_ensemble_5MHz','5 MHz','#2674a6'),('final_ensemble_10MHz','10 MHz','#b95034')]:
        summary=load(root/folder/'forward_summary.json');summaries.append(summary)
        with (root/folder/'forward_receivers.csv').open(encoding='utf-8-sig') as stream:rows=list(csv.DictReader(stream))
        x=np.array([float(r['x_m']) for r in rows]);a=np.array([complex(float(r['real']),float(r['imag'])) for r in rows])
        ref=np.array([complex(float(r['reference_real']),float(r['reference_imag'])) for r in rows]);ratio=a/ref
        logs=np.log(abs(ratio));x=(x-x[0])*1000
        axes[0,0].plot(x,logs-logs[0],'o',color=color,label=label)
        fitted=np.polyfit(x,logs,1);axes[0,0].plot(x,np.polyval(fitted,x)-np.polyval(fitted,x[0]),color=color)
    axes[0,0].set(xlabel='传播距离（mm）',ylabel='前行波幅值比的自然对数变化',title='10个随机样本的相干前行波');axes[0,0].legend();axes[0,0].grid(alpha=.25)
    for i,s in enumerate(summaries):
        value=s['fit']['alpha_db_mm'];lo,hi=np.array(s['bootstrap_ci_95_np_m'])*8.685889638/1000
        axes[0,1].errorbar([s['frequency_hz']/1e6],[value],yerr=[[value-lo],[hi-value]],fmt='o',capsize=5)
    axes[0,1].set(xlabel='频率（MHz）',ylabel='幅值衰减系数（dB/mm）',title='示例材料：自助重采样95%置信区间');axes[0,1].grid(alpha=.25)
    hold=load(root/'held_out/verification.json');d=np.load(root/'held_out/comparison.npz')
    axes[1,0].plot((d['x_m']-d['x_m'][0])*1000,abs(d['observed']),'o',label='未参与标定的随机样本')
    axes[1,0].plot((d['x_m']-d['x_m'][0])*1000,abs(d['prediction']),label='标定后的等效传播模型')
    axes[1,0].set(xlabel='传播距离（mm）',ylabel='归一化前行波幅值',title=f'未参与标定的局部模型：幅值误差{100*hold["ensemble_max_magnitude_error"]:.2f}%');axes[1,0].legend();axes[1,0].grid(alpha=.25)
    storage=load(root/'comsol_transient_initial/storage_compaction.json')
    random=next(r for r in storage if r['model'].endswith('random_pores.mph'))
    axes[1,1].bar(['保留全部有限元场','末帧有限元场＋全部接收历史'],[random['bytes_before']/1e6,random['bytes_after']/1e6],color=['#888888','#2674a6'])
    axes[1,1].set(ylabel='模型文件大小（MB）',title='完整保留601帧接收历史')
    fig.savefig(root/'simulation_summary.png',dpi=180);plt.close(fig)
    gate_path=root/'revalidation_corrected_domain/forward_gates.json'
    if not gate_path.exists():gate_path=root/'final_validation/forward_gates.json'
    meshgates=load(gate_path);scalar=load(root/'equivalent_bulk_5MHz/verification.json')
    macro=load(root/'equivalent_bulk_5MHz/parameters.json');three=load(root/'comsol_3d_mesh_based/results.json')
    time_c=np.load(root/'time_coarse/no_pore_transient.npz');time_f=np.load(root/'time_fine/no_pore_transient.npz')
    interpolated=np.array([[np.interp(time_c['times_s'],time_f['times_s'],s) for s in component.T] for component in time_f['displacement']]).transpose(0,2,1)
    time_error=np.linalg.norm(interpolated-time_c['displacement'])/np.linalg.norm(interpolated)
    guide=load(root/'comsol_guided_corrected/forward_results.json');air=load(root/'air_comparison_corrected/comparison.json')
    rows=[]
    for s in summaries:
        fit=s['fit'];ci=np.array(s['bootstrap_ci_95_np_m'])*8.685889638/1000
        rows.append(f'| {s["frequency_hz"]/1e6:g} MHz | {fit["alpha_db_mm"]:.5f} | [{ci[0]:.5f}, {ci[1]:.5f}] | {fit["r_squared"]:.4f} |')
    text=f'''# 本次实现与示例仿真核验

已实现论文的高斯随机场、排序二值化孔隙模型，并在 COMSOL 6.2 中完成局部体波频域／瞬态、A0/S0、空腔与空气对照、小规模三维体网格模型，以及单方向标量传播替代模型。14项科学计算检查通过。

后续单孔／孔群透射对照、保存解复提取以及孔隙率验收修正见[数值复核报告](revalidation/REPORT.md)。存在该报告时，以其中的修正验收为准。

**所有结果使用明确标注的示例补全材料卡；不是用户实际 CFRP 的标定结果。** 未获得原失败 .mph 文件，因此本次验证的是新建模型，不能断言已定位原模型不收敛的原因。

![本次数值结果](simulation_summary.png)

## 可复现输出

- `sample/`：二维随机场、像素／矢量孔隙、实际形貌、相关性、分散／聚集对照。
- `sample_3d/`、`comsol_3d_refined/random_geometry/`：三维随机场、封闭孔隙曲面、体积校准和体网格。
- `calibration_demo/`：同时匹配水平／垂直尺寸；演示目标40、20 μm，所得平均40.281、19.375 μm。
- `final_ensemble_5MHz/`、`final_ensemble_10MHz/`：每组10个独立孔隙模型、无孔参考、前行波提取和95%置信区间。
- `comsol_guided_corrected/`：A0/S0参考与随机孔隙模型、前后方向分离与非模态残差。
- `final_validation/`：轮廓、网格、局部长度／宽度检查。
- `held_out/`：未参与标定的211、212号随机种子。
- `equivalent_bulk_5MHz/equivalent_P_5MHz.mph`：71自由度、单方向／单频率的相干P波标量模型。
- `comsol_3d_mesh_based/random_pores.mph`：不经过CAD拟合的三维局部固体力学模型。

## 关键优化及证据

1. 扰动标准差不再做冗余迭代。固定随机场的正比例缩放不改变排序孔形；孔隙率使用比例输入。
2. 实际矢量面积／三维曲面体积单独校准，不把像素／体素计数直接当作COMSOL孔隙率。
3. 各向异性方向阻抗明确使用全局坐标系。默认边界坐标系会把阻抗分量放错方向；无孔体波基准幅值变异系数由约23%降低至0.00087%。
4. **前行波与反射波必须分开。** 对体波联合位移、法向应力，使用参考介质阻抗投影；用COMSOL切线单元积分，避免等间距采样跨越微孔时的积分误差。
5. 导波联合两个位移和两个牵引应力分量，分别投影S0/A0的前行／后行分支，同时保存未解释的场残差。
6. 三维STL转CAD会产生细小面、薄缝和退化网格。本次改用TetGen构造保留孔面和孔腔的体网格，再作为COMSOL网格几何导入，绕过CAD曲面拟合。保留13个封闭孔腔，实际孔隙率{100*three['actual_porosity']:.6f}%，299751个四面体，{three['dofs']}个自由度，求解约{three['solve_seconds']:.1f}秒。它是小规模求解验证，不是三维总体衰减已收敛的证明。
7. 全部接收历史先输出NPZ，再仅保留最后一个有限元场及导数。随机孔瞬态模型从{random['bytes_before']/1e9:.3f} GB降至{random['bytes_after']/1e6:.1f} MB；601帧接收历史和最后一帧位移均重新打开核对。

## 体波示例数据

孔隙率名义1.5%；相关长度20、10 μm；单向示例材料；两个频率复用同一组随机孔隙统计。指标为相干前行波的幅值衰减，不能等同于全部能量吸收。

| 频率 | 衰减 dB/mm | Bootstrap 95% CI | 指数拟合 R² |
|---|---:|---:|---:|
{chr(10).join(rows)}

该结果支持本示例中提高频率后整体响应增强；两个频率点不足以建立普适幂次律，也不能把二维结果直接当成真实三维CFRP结果。

## 验收与实用边界

- 无孔瞬态内部步长减半，接收波形相对L2差异{100*time_error:.3f}%。
- 5 MHz前行波的网格、轮廓变化分别约{100*meshgates['gates']['mesh_fine']['relative_change']:.2f}%、{100*meshgates['gates']['contour_full']['relative_change']:.2f}%，均小于5%。
- **长度／宽度变化下，小非零衰减系数的相对变化仍很大。** 当前仅满足0.01 Np的近零累计误差门限，不能宣称衰减系数本身已达到5%统计代表性。
- 两个未见随机样本的平均幅值预测最大误差{100*hold['ensemble_max_magnitude_error']:.2f}%，复幅值最大误差{100*hold['ensemble_max_complex_amplitude_error']:.2f}%；验证跨度仅{hold['tested_distance_span_m']*1000:.2f} mm。
- 标量模型相对于自身解析传播解，复幅值最大误差{100*scalar['max_complex_amplitude_error']:.4f}%，衰减相对误差{100*scalar['relative_alpha_error']:.4f}%。71自由度对应约28万自由度的显式局部案例，但**只预测标定方向的平均波，不能复现微孔散射场、尾波或所有弹性模态**。
- A0/S0随机孔模型的非模态残差约{100*next(r for r in guide if r['mode']=='A0' and r['case']=='random_pores')['nonmodal_residual_max']:.1f}%、{100*next(r for r in guide if r['mode']=='S0' and r['case']=='random_pores')['nonmodal_residual_max']:.1f}%；当前不足以批准导波宏观衰减外推。需增加随机集合，并考虑在均匀材料截面接收、改变孔隙段长度提取传输衰减。
- 半径20 μm单孔、5 MHz时，空腔与空气模型的最大幅值差异{100*air['max_relative_amplitude_difference']:.4f}%；不覆盖任意尺寸的孔腔共振。
- 体波与导波演示采用了不同孔径统计，不能直接比较其数值作为同一孔隙群的检测性能。

## 后续正式研究

代码已提供0%、0.1%、0.5%、1%、2%、3%、4%、5%的扫描、种子复用、实际形貌统计、形貌参数入口、同孔形分布对照和失败验收。尚未运行所有材料／形貌／分布组合，也未生成未经验证的大试件定量预测。正式研究需替换完整材料卡、明确铺层和孔径分布；先建立统计代表性，再开展全量扫描。

模型改用不同孔隙率、形貌、分布或频率后需要重新辨识传播参数。当前替代模型的有效范围是5 MHz、P波、x方向及本示例孔隙统计；禁止自动映射到A0/S0或全张量本构。

建模依据：用户指定的张翔（2011）论文；COMSOL 6.2本机API文档及[瞬态波指南](https://www.comsol.com/support/knowledgebase/1118)。三维体网格使用[MeshPy/TetGen官方接口](https://documen.tician.de/meshpy/tri-tet.html)。本报告由脚本根据实际输出自动整理。
'''
    (root/'REPORT.md').write_text(text,encoding='utf-8')
    validation=dict(unit_checks_passed=14,time_step_relative_l2=float(time_error),scalar_model=scalar,
                    heldout=hold,three_dimensional=three,mesh_and_contour_gates=meshgates)
    (root/'verification_summary.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')


if __name__=='__main__':main()
