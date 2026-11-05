"""Generate Chinese figures and an auditable report from completed native solves."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from guided_collective import OUT,read
from guided_statistics import summarize,normalized_transmission,coherent_interval
from guided_extend_statistics import summarize as summarize_extended,family_path
from guided_length_study import summarize as summarize_lengths

plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'figure.dpi':150})
COLORS=['#2563a5','#c76121','#37866a']


def save(fig,name):
    fig.savefig(OUT/name,bbox_inches='tight');plt.close(fig)


def main():
    sweep=summarize(OUT);length=summarize_lengths(6);extended=summarize_extended(12)
    assert all(r['samples']==12 for r in extended),'Twelve seed families are required.'
    checks=read(OUT/'numerical_checks.json');reference=read(OUT/'reference/results.json')
    single=[]
    for d in [20,40,80]:
        paths=list((OUT/f'd{d}').glob('single_*/results.json'))+list((OUT/f'd{d}').glob('global_maximum_*/results.json'))
        assert len(list((OUT/f'd{d}').glob('global_maximum_*/results.json')))==3
        for f in [50000.,100000.,200000.]:
            for j,mode in enumerate(['S0','A0']):
                transmissions=[];reflections=[]
                for p in paths:
                    record=read(p);transmissions.append((abs(1-abs(normalized_transmission(record,reference,f)[j])),str(p.parent.relative_to(OUT))))
                    a=next(r for r in record['results'] if r['frequency_hz']==f)
                    b=next(r for r in reference['results'] if r['frequency_hz']==f)
                    R=np.array(a['reflection_real'])+1j*np.array(a['reflection_imag'])-np.array(b['reflection_real'])-1j*np.array(b['reflection_imag'])
                    reflections.append(float(np.max(abs(R[:,j]))))
                worst,case=max(transmissions)
                error=max(r['mesh_magnitude_error'] for r in checks['transmission_checks'] if r['case']==f'd{d}_maximum' and r['frequency_hz']==f and r['mode']==mode)
                single.append(dict(diameter_group_um=d,frequency_hz=f,mode=mode,max_transmission_change=float(worst),worst_case=case,
                    max_reflected_modal_change=max(reflections),representative_mesh_difference=error,
                    transmission_screen_with_representative_mesh_margin=bool(worst+error<.001),
                    note='Mesh margin comes from representative maximum, not a certified error bound for every translated shape.'))
    (OUT/'single_controls_summary.json').write_text(json.dumps(single,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),sharey=True)
    for j,mode in enumerate(['S0','A0']):
        ax=axes[j]
        for d,c in zip([20,40,80],COLORS):
            rows=[r for r in sweep if r['mode']==mode and r['diameter_group_um']==d]
            x=np.array([r['frequency_hz']/1000 for r in rows]);y=np.array([r['amplitude_reduction']*100 for r in rows])
            ci=np.array([r['ci95_amplitude_reduction'] for r in rows])*100
            ax.plot(x,y,'o-',color=c,label=f'孔径组 {d} 微米');ax.fill_between(x,ci[:,0],ci[:,1],color=c,alpha=.13)
        ax.axhline(.1,color='#a22',ls='--',label='示例阈值 0.1%');ax.axhline(0,color='.6',lw=.7)
        ax.set(title=f'{mode} 模态：8 毫米孔隙段',xlabel='激励频率（千赫）',xticks=[50,100,200]);ax.grid(alpha=.2)
    axes[0].set_ylabel('相干集合透射幅值下降（%）');axes[1].legend(fontsize=8)
    fig.suptitle('孔径与频率扫描：每组 6 个随机样本，阴影为 95% 自助法区间');fig.tight_layout();save(fig,'孔径频率扫描.png')

    fig,axes=plt.subplots(1,2,figsize=(11,4.3),sharey=True)
    contrasts=[]
    for j,mode in enumerate(['S0','A0']):
        values=[]
        for L in [8,24,48]:
            ref=read(OUT/('reference' if L==8 else f'reference_L{L}')/'results.json')
            signals=np.array([normalized_transmission(read(family_path(L,seed)/'results.json'),ref,200000.)[j] for seed in range(101,113)])
            summary=coherent_interval(signals);contrasts.append(dict(mode=mode,length_mm=L,**summary));values.append(signals)
        values=np.array(values).T;x=np.array([8,24,48]);ax=axes[j]
        for row in values:ax.plot(x,100*(1-abs(row)),color='.7',alpha=.45,lw=.8)
        rows=[r for r in contrasts if r['mode']==mode];means=np.array([r['amplitude_reduction'] for r in rows])*100
        ci=np.array([r['ci95_amplitude_reduction'] for r in rows])*100
        ax.errorbar(x,means,yerr=[means-ci[:,0],ci[:,1]-means],fmt='o-',color=COLORS[j],capsize=4,label='集合均值与 95% 区间')
        ax.axhline(.1,color='#a22',ls='--',label='示例阈值 0.1%');ax.set(title=f'{mode} 模态',xlabel='孔隙段长度（毫米）',xticks=x);ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes[0].set_ylabel('归一化透射幅值下降（%）');fig.suptitle('80 微米孔径组、200 千赫：12 组样本的长度验证\n灰线为单个样本，集合均值保留复数相位');fig.tight_layout();save(fig,'孔群长度与统计区间.png')

    fig,axes=plt.subplots(1,2,figsize=(11,4.3),sharey=True)
    for j,mode in enumerate(['S0','A0']):
        ax=axes[j]
        for d,c in zip([20,40,80],COLORS):
            rows=[r for r in single if r['mode']==mode and r['diameter_group_um']==d]
            ax.plot([r['frequency_hz']/1000 for r in rows],[r['max_transmission_change']*100 for r in rows],'o-',color=c,label=f'{d} 微米组：透射')
            ax.plot([r['frequency_hz']/1000 for r in rows],[r['max_reflected_modal_change']*100 for r in rows],'x--',color=c,label=f'{d} 微米组：反射')
        for threshold in [.01,.1,1]:ax.axhline(threshold,color='.5',ls=':',lw=.8)
        ax.set(yscale='log',title=f'{mode} 入射：单孔对照中的最大变化',xlabel='激励频率（千赫）',xticks=[50,100,200]);ax.grid(alpha=.2)
    axes[0].set_ylabel('归一化模态幅值变化（%）');axes[1].legend(fontsize=7,ncol=2)
    fig.suptitle('透射低于阈值，不能推出反射也低于阈值\n点线分别为 0.01%、0.1%、1% 示例阈值');fig.tight_layout();save(fig,'单孔透射与反射判据.png')

    allrecords=[]
    for p in OUT.rglob('results.json'):
        if 'invalid_cad_before_repair' not in p.parts:allrecords.append((p,read(p)))
    rawrows=[r for _,record in allrecords for r in record['results']]
    residual=max(r['uniform_port_residual_max'] for r in rawrows)
    energy=max(abs(v-1) for r in rawrows for v in r['energy_balance_by_incident_mode'])
    baseline=read(OUT/'baseline_mode_validation.json');fe=read(OUT/'fe_port_check.json');basis=read(OUT/'basis_convergence.json')
    thresholds=[]
    for r in extended:
        mode=r['mode'];control=next(c for c in single if c['diameter_group_um']==80 and c['frequency_hz']==200000 and c['mode']==mode)
        numerical=max(c['mesh_magnitude_error'] for c in checks['transmission_checks'] if c['case']=='d80_L48' and c['frequency_hz']==200000 and c['mode']==mode)
        for t in [.0001,.001,.01]:
            thresholds.append(dict(mode=mode,threshold=t,single_below=control['max_transmission_change']<t,
                group_ci_lower_above=r['longest_segment_contrast']['ci95_amplitude_reduction'][0]>t,
                group_ci_lower_minus_representative_mesh_above=r['longest_segment_contrast']['ci95_amplitude_reduction'][0]-numerical>t,
                reflected_single_below=control['max_reflected_modal_change']<t))
    (OUT/'threshold_sensitivity.json').write_text(json.dumps(thresholds,indent=2),encoding='utf-8')
    def yes(v):return '通过' if v else '未通过'
    lines=['# 薄板孔群效应：二维 COMSOL 验证结果',
      '', '本轮预定的频率/孔径扫描、单孔对照、模态分离和统计复核已完成。结果支持在指定透射幅值观测量下观察到累积孔群响应；不支持“单孔在所有信号通道均不可检测”，也尚未达到衰减系数随样本数变化小于 5% 的稳定性目标。',
      '', '## 模型和实际计算范围',
      '', '- 二维平面应变板，厚 1 mm，总长约 338 mm；孔隙起点约 92 mm，吸收区从约 154 mm 开始。上下表面及孔壁自由，孔腔由几何扣除形成。吸收区较长是为控制低频反射。',
      '- 50、100、200 kHz 分别进行 S0、A0 频域谐波激励。本轮不是时域脉冲仿真，因此没有时步收敛结论。',
      '- 孔径组为等效面积直径约 20、40、80 μm，组内有分布，不是三个固定直径。8 mm 孔隙段分别含 377、94、24 个孔；各孔径组均有 6 个随机样本，孔隙段面积孔隙率 1.5%。',
      '- 80 μm 组补充 24、48 mm 孔隙段，分别约 71、141 孔；100 kHz 使用 6 组样本，200 kHz 扩展到 12 组，三种长度各 12 个模型。',
      '- 检查代表性最小、中位、最大单孔，并将最大孔移至上下表面附近。另从本轮全部已生成的组样本中选取实际最大孔，在板中部和上下表面附近复核。三位置不是对所有形貌、取向与位置的穷尽证明。',
      '- 孔形来源是高斯随机场生成的不规则轮廓库；按直径缩放后重新无重叠放置。保留了轮廓，不保留原随机场的位置自相关。',
      '- 材料是明确标注的示例补全 CFRP 参数，零材料本底阻尼、理想连续固体；不能作为实测铺层或三维孔隙的定量预测。二维孔腔对应面外延伸缺陷。',
      f'- 共保留 {len(allrecords)} 个有效算例记录，{len(rawrows)*2} 次已记录模态求解。每个 model.mph 只保存最后一次频率的 A0 解；所有已计算频率及两种入射的端口数据在 signals.npz 和 results.json 中。',
      '', '## 模态分离的验证',
      '', '在孔隙段前后各取三个完整无孔截面，联合位移和牵引应力分离 S0/A0 的前行、后行波。各截面先补偿传播相位，再平均；两次独立模态激励组成入射/出射矩阵，消除实际入射中混合模态及源端反馈的影响。',
      '', f'- 无孔参考相速度与半解析有限元模态的最大相对差：{max(r["relative_phase_velocity_difference"] for r in baseline)*100:.5f}%。',
      f'- 厚度模态基函数从 120 单元加密至 240 单元，最大波长变化为 {max(r["relative_wavelength_difference"] for r in basis["rows"])*100:.5f}%；第一高阶模态截止频率约 {basis["first_higher_cutoff_hz"]/1e6:.3f} MHz，高于本次 0.05～0.2 MHz 频段。',
      f'- 全部有效记录的最大无孔端口拟合残差：{residual*100:.3f}%；残差是场拟合量，不等同于接收幅值误差。',
      f'- 无孔参考最大远端反向入射幅值比：{max(r["right_incoming_ratio_max"] for r in reference["results"])*100:.4f}%。',
      f'- 模态功率平衡最大偏差：{energy*100:.4f}%。',
      f'- 代表性 A0/200 kHz 场用 COMSOL 截面积分的互易投影独立复核；与同一直接透射观测量的最小二乘投影最大差：{max(r["method_magnitude_difference"] for r in fe)*100:.5f} 个百分点。该复核不把单次激励比值与去嵌入的双入射矩阵混为同一观测量。',
      '', '因此，原先在含孔截面取样导致的大残差问题已修正；本轮模态分离有无孔频散、前后行波、互易投影及功率平衡的交叉检查。验证范围是当前均匀示例板与频段。',
      '', '## 孔群响应与统计衰减',
      '', '主观测量为同一模态的归一化复透射系数 T（除以匹配无孔参考）。幅值下降为 1−|平均T|；同时记录每个样本的 1−|T|，区分集合相位相消和单个试件的幅值损失。负下降值原样保留。',
      '', '对 8、24、48 mm 孔隙段拟合 ln|平均T(L)| = b−αL，保留截距 b。不同长度按同一随机种子重新生成，灰线不是同一试件逐渐延长。α 是本测试区间内的表观相干衰减；三个接收截面用来稳定模态提取，不能用孔外无孔段的幅值斜率代替孔隙衰减。95% 区间按完整随机种子族自助重采样 1000 次，仅反映样本统计不确定性。',
      '', '| 80 μm 组，200 kHz，12 样本 | S0 | A0 |','|---|---:|---:|']
    a,b=extended
    metrics=[('48 mm 集合幅值下降',lambda r:f'{r["longest_segment_contrast"]["amplitude_reduction"]*100:.3f}%'),
             ('上述下降的 95% 区间',lambda r:' ～ '.join(f'{v*100:.3f}%' for v in r['longest_segment_contrast']['ci95_amplitude_reduction'])),
             ('α（Np/m）',lambda r:f'{r["fit"]["alpha_np_m"]:.5f}'),
             ('α 的 95% 区间（Np/m）',lambda r:' ～ '.join(f'{v:.5f}' for v in r['ci95_alpha_np_m'])),
             ('逐样本幅值定义的 α（Np/m）',lambda r:f'{r["typical_fit"]["alpha_np_m"]:.5f}'),
             ('三长度指数拟合 R²',lambda r:f'{r["fit"]["r_squared"]:.4f}'),
             ('样本数 6→12 的 α 变化',lambda r:f'{r["sample_count_relative_change"]*100:.1f}%'),
             ('样本数变化 ≤5%',lambda r:yes(r['sample_count_5percent_passed']))]
    lines += [f'| {label} | {fn(a)} | {fn(b)} |' for label,fn in metrics]
    lines += ['', 'α 的正区间说明本轮存在随孔隙段长度累积的统计幅值损失；样本数稳定性未通过意味着这些 α 仍是有区间的局部估计，不能作为已收敛的材料常数。该结果也没有证明非线性协同作用或独有的多重散射增益。',
      '', '![孔群长度与统计区间](孔群长度与统计区间.png)', '', '![孔径频率扫描](孔径频率扫描.png)',
      '', '## 单孔条件与阈值敏感性',
      '', '0.1% 是用户选定的示例阈值，并非实验噪声。透射只比较幅值，不能据此判断相位、反射、转换模态或成像可检测性。',
      '', '| 80 μm 组、200 kHz，单孔对照的最不利观测 | S0 入射 | A0 入射 |','|---|---:|---:|']
    controls=[next(c for c in single if c['diameter_group_um']==80 and c['frequency_hz']==200000 and c['mode']==m) for m in ['S0','A0']]
    for label,key in [('透射幅值变化','max_transmission_change'),('反射模态幅值变化','max_reflected_modal_change')]:
        lines.append(f'| {label} | {controls[0][key]*100:.4f}% | {controls[1][key]*100:.4f}% |')
    lines += ['', '| 模态 | 示例阈值 | 全部已测单孔透射低于阈值 | 孔群 95% 下界扣除代表性网格差后高于阈值 | 单孔反射也低于阈值 |','|---|---:|---|---|---|']
    lines += [f'| {r["mode"]} | {r["threshold"]*100:g}% | {yes(r["single_below"])} | {yes(r["group_ci_lower_minus_representative_mesh_above"])} | {yes(r["reflected_single_below"])} |' for r in thresholds]
    lines += ['', '0.1% 下得到的是“单孔透射幅值变化较小、长孔群透射损失较大”的窗口。由于最大单孔的反射超过阈值，严格的“单孔完全不可检测而孔群可检测”条件未满足。也不能用孔径小于波长替代实际可检测性判断。',
      '', '![单孔透射与反射判据](单孔透射与反射判据.png)',
      '', '## 数值误差与未通过项',
      '', '背景最大网格从 150 μm 缩至 75 μm，孔边最小尺寸从代表孔径的 1/30 缩至 1/60，并对细网格设置 1 nm 绝对几何修复容差。下表是同一 seed101 的三长度拟合检查，不是整个随机集合的网格误差上界。',
      '', '| 频率 | 模态 | 代表性 α 的网格变化 | 5% 检查 |','|---|---|---:|---|']
    lines += [f'| {r["frequency_hz"]/1000:g} kHz | {r["mode"]} | {r["relative_alpha_change"]*100:.2f}% | {yes(r["relative_alpha_mesh_passed"])} |' for r in checks['alpha_checks']]
    lines += ['', '- 100 kHz 的 A0 网格检查未通过；200 kHz 的两模态在代表样本上通过。网格检查不能代替随机样本收敛。',
      '- 全局最大孔及其移动位置没有逐例做成对加密；敏感性表中的误差余量来自代表性加密对照，属于实测数值差估计，不是严格误差界。0.01% 级别的判断尤其不能视为充分验证。',
      '- 发现默认相对几何容差合并邻近孔：20 μm/seed106、40 μm/seed105 两例已保留旧输出并以 1 nm 绝对容差重新计算；扩展序列遇到同类问题时自动重建。当前统计不包含旧失效输出。',
      '- 尚未完成轮廓逐级加密、所有随机样本的网格收敛、独立吸收层长度变化、三维薄板与真实铺层验证。因此没有把当前结果标记为全项验证通过。',
      '', '## 在 COMSOL 中复核',
      '', '可打开 `d80/L48/seed_101/model.mph` 查看长孔群的 A0/200 kHz 解；打开 `reference_L48/model.mph` 查看相同截面配置的无孔参考。原始孔轮廓及参数在各算例的 `geometry/`。模型中的源函数按模态设置；简单改激励频率并不会自动刷新正确的模态形状，扫频请运行配套脚本或同步更新源函数。',
      '', '复现入口：`guided_collective.py`（基础扫描）、`guided_length_study.py`（长度）、`guided_extend_statistics.py`（12 样本）、`guided_global_maximum.py`（实际最大孔）、`guided_numerical_checks.py`（数值检查）、`make_guided_collective_report.py`（本报告）。已完成的匹配算例可复用缓存。',
      '', '所有数值均来自本地 COMSOL 6.2 实际求解。完整机器可读结果：`sweep_summary.json`、`extended_attenuation.json`、`numerical_checks.json`、`single_controls_summary.json`、`threshold_sensitivity.json`、`fe_port_check.json`、`geometry_audit.json`。',
      '', '下一步若要把 α 用作定量反演参数，应先替换完整实测材料卡，并围绕当前候选频段增加独立样本、完成集合网格及吸收层验证。当前结果足以支持开展透射孔群实验，但不足以承诺真实设备中单孔不可见或低频必然能检出任意低孔隙率。']
    dimensions=['## 实际孔径与波长','','| 名义孔径组 | 全部孔群的实际等效直径范围 |','|---|---:|']
    for d in [20,40,80]:
        paths=list(OUT.glob(f'd{d}/L*/seed_*/geometry/metadata.json'))+list(OUT.glob(f'extended_samples/d{d}/L*/seed_*/geometry/metadata.json'))
        metadata=[read(p) for p in paths]
        dimensions.append(f'| {d} μm | {min(r["diameter_min_um"] for r in metadata):.2f}～{max(r["diameter_max_um"] for r in metadata):.2f} μm |')
    dimensions+=['','| 频率 | S0 波长 | A0 波长 |','|---|---:|---:|']
    for f in [50000.,100000.,200000.]:
        waves=[next(r for r in basis['rows'] if r['mode']==m and r['frequency_hz']==f)['wavelength_m']*1000 for m in ['S0','A0']]
        dimensions.append(f'| {f/1000:g} kHz | {waves[0]:.3f} mm | {waves[1]:.3f} mm |')
    dimensions+=['','孔径远小于波长只表明本例处于亚波长尺度；是否可检测仍按各通道幅值与示例阈值判断。','']
    position=lines.index('## 模态分离的验证');lines[position:position]=dimensions
    (OUT/'验证结果_中文.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('Written',OUT/'验证结果_中文.md')


if __name__=='__main__':main()
