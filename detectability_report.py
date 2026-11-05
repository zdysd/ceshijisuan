"""Scientific CSV/JSON/Markdown and figures, explicitly separated by evidence stage."""
import os
import numpy as np
from detectability_200kHz_L48 import OUT,ROOT,MODES,TAU,read,write,write_csv,stamp,cp_lower

def render():
    p0=read(OUT/'p0_audit.json');protocol=read(OUT/'protocol.json');models=read(OUT/'p1_models.json')
    state=read(OUT/'worker_state.json') if (OUT/'worker_state.json').exists() else dict(stage='imported_snapshot_not_running')
    resource_policy=read(OUT/'resource_policy.json') if (OUT/'resource_policy.json').exists() else {}
    stages={}
    for stage in ['p1','p2','p3']:
        path=OUT/f'{stage}_validation.json';stages[stage]=read(path) if path.exists() else dict(status='pending',passed=None)
    if state['stage'] in ['imported_snapshot_not_running','paused_for_shutdown']:
        for gate in stages.values():
            if gate.get('status') in ['running','running_waiting_for_resources']:gate['status']='paused_snapshot'
    status=dict(updated_at=stamp(),P0=dict(status='complete',passed=True),P1=stages['p1'],P2=stages['p2'],P3=stages['p3'],
        P4=dict(status='requires_experimental_data',passed=None,experimental_detection_claim_allowed=False),worker=state,
        interpretation='Fixed illustrative threshold in a specified numerical geometry distribution; never experimental POD')
    write(OUT/'validation.json',status)
    rows=read(OUT/'per_case_metrics.json')
    additions=[]
    for path in sorted((OUT/'metrics').glob('*.json')) if (OUT/'metrics').exists() else []:additions.extend(read(path))
    if additions:
        rows=[r for r in rows if r['category'] in ['development_group','legacy_single']]+additions
        write(OUT/'per_case_metrics.json',rows);write_csv(rows)
    lines=['# 200 kHz、S0/A0、48 mm 验证报告','',f'更新时间：{status["updated_at"]}','',
        f'当前运行状态：{state["stage"]}。P0 已完成；后续阶段只按实际验收记录标记，不把排队或预备脚本视为已经求解。','',
        '主阈值固定为无量纲 0.001，即相对参考幅值的 0.1%。loss=1−|T| 保留符号；score=|1−|T||。所有 T 均由双入射矩阵去混合后，除以同端口无孔参考。','',
        '|模态|开发集逐样本原始越阈数|开发数据描述性单侧95%下限|相干集合下降|单样本下降均值|旧单孔最大 score|',
        '|---|---:|---:|---:|---:|---:|']
    for s in p0['development_statistics']:
        lines.append(f'|{s["mode"]}|{s["raw_score_successes"]}/12|{100*s["one_sided95_lower_descriptive"]:.3f}%|{100*s["coherent_amplitude_loss"]:.6f}%|{100*s["mean_individual_loss"]:.6f}%|{100*s["strongest_legacy_single_score"]:.6f}%|')
    lines+=['','P1 completed entries indicate numerical acceptance, not necessarily threshold exceedance. Read per_case_metrics for each mode and specimen.','', '11/12 的描述性下限约66.1%，不能宣称已达到90%的可靠越阈目标。旧单孔使用其原匹配参考；其孔纵向位置和端口不是新48 mm控制配置，所以只作为开发背景。低响应 seed102 完整保留。','',
        '## 数值证据与缺项','',
        '|代表性 seed101 的既有48 mm网格对照|幅值差（百分点）|非零损失相对变化|相位差（°）|',
        '|---|---:|---:|---:|']
    for r in p0['mesh_checks']:lines.append(f'|{r["mode"]}|{100*r["magnitude_delta"]:.6g}|{100*r["relative_loss_change"]:.4g}%|{r["phase_error_deg"]:.6g}|')
    lines+=['','上表只覆盖 seed101 的该网格对照，不能充当所有模型的误差上界。完整 u_num 尚需逐模型网格、吸收区和提取敏感性；缺少项时不填零，不给保守检出结论。',
        '原有限元模型通常只存最后频率的 A0 场；S0 先复用原始端口数组，不能用 A0 场替代。8 mm 已通过瞬态结果仅复用实现与设置，不作为48 mm验收证据。',
        '形貌上游仅可追溯到固定折线库；本轮限定为固定折线几何条件下的验证，不把加密插点称为独立轮廓收敛。','',
        '## 固定独立验证规则','',
        f'新样本已冻结为 {protocol["independent_seeds"][0]}–{protocol["independent_seeds"][-1]}，共30个。P1/P2通过之前仅生成几何，不计算或查看这些样本的 T。',
        f'精确单侧95%二项下限：30/30 为 {100*cp_lower(30,30):.3f}%，29/30 为 {100*cp_lower(29,30):.3f}%。预设目标为下限≥90%；不追加种子、换样本或挑选模态。',
        'S0/A0 分开验收，数值误差使分类待定时保守计为未成功。旧12个开发样本不进入固定30个独立验证分母。','',
        '## 阶段状态','']
    for name in ['P0','P1','P2','P3','P4']:lines.append(f'- {name}：{status[name]["status"]}；通过标记 {status[name].get("passed")}。')
    for name,gate in stages.items():
        if gate.get('passed') is not False:continue
        lines+=['',f'## {name.upper()} 未通过项与诊断','',
            f'已记录原因：{gate.get("reason","预设验收或独立验证目标未达到")}。',
            f'[完整阶段证据]({(OUT/(name+"_validation.json")).as_posix()})']
        failed=gate.get('failed')
        if failed:
            lines+=['',f'失败模型：{failed.get("case_id",failed.get("case"))}；保留到细化级 {failed.get("level","见证据")}。',
                '```json',__import__('json').dumps(failed,ensure_ascii=False,indent=2),'```']
        lines+=['',
            '证据只能确认哪些指标未满足预设要求，不能单凭失败标记唯一归因。若网格差占主导，核对孔边/厚度离散和CAD孔数、面积；若吸收区差占主导，在孔区与端口不动的条件下检查反向入射和端部反射；若提取差占主导，检查厚度采样与同一双入射观测量的互易积分。',
            '若280/560 μs窗口或半时间步失败，核对实际到达顺序、完整尾波及频谱相位；达到记录上限不再延长。若数值检查通过而透射裕量仍不足，应如实判为低于阈值或待定，而非修改阈值。',
            '依赖阶段停止，全部版本和失败样本保留。任何下一轮方案需重新预先固定，不用追加种子把本轮结果凑成通过。']
    if stages['p3'].get('statistics'):
        lines+=['','## 固定30样本的独立数值验证','',
            '|模态|保守成功数/30|精确单侧95%下限|达到90%目标|','|---|---:|---:|---|']
        for v in stages['p3']['statistics']:lines.append(f'|{v["mode"]}|{v["successes"]}/30|{100*v["one_sided95_lower"]:.4f}%|'+('通过' if v['passed'] else '未通过')+'|')
        for i,mode in enumerate(MODES):
            below=[r['seed'] for r in stages['p3']['completed'] if r['classifications'][i]=='below']
            pending=[r['seed'] for r in stages['p3']['completed'] if r['classifications'][i]=='indeterminate']
            lines.append(f'{mode}：低于阈值样本 {below}；待定且按未成功计的样本 {pending}。')
        lines+=['','数值敏感性中吸收区项来自代表孔群的已测包络，逐新样本实测网格及采样项；该覆盖限制必须随结论保留，不是严格数学误差界。独立数值验证不等于真实实验POD。',
            '![独立样本与数值敏感性](independent_thresholds.png)','',
            '![独立数值越阈统计](independent_probability.png)']
    if state.get('error'):lines+=['',f'停止／等待原因：{state["error"]}']
    lines+=['','## 实验适用性（P4）','',
        '缺少实测噪声、重复拆装/耦合数据、传感器电压及孔径/频响模型、温度与增益变化、实测材料卡/铺层/阻尼、三维孔形及孔隙真值。四场变量理想投影不能直接视为单个传感器电压。','',
        '实验应先在独立无孔校准数据上锁定阈值，目标误报率≤5%；再用未参与校准的数据估计误报率及其单侧95%上限。孔群实验POD目标为单侧95%下限≥90%；单孔低检出能力需另设例如单侧95%上限≤10%的目标，两者使用相同设备、平均次数和检测规则。',
        '实测阈值高于0.1%时必须按真实阈值重判；指标或设备改变后重新校准并独立验证。没有这些数据时，不添加合成白噪声冒充实验验证。','',
        '本轮最终结论最多涉及固定二维示例材料、孔形生成分布和示例阈值下的数值越阈能力。反射、转换和相位另列，不能将单孔透射低于阈值写成所有通道不可检测。','',
        '[逐模型数值](per_case_metrics.csv) · [完整机器可读指标](per_case_metrics.json) · [复用与来源哈希](reuse_manifest.json) · [冻结协议](protocol.json)','',
        '![开发集逐样本与阈值](development_thresholds.png)','',
        '![单孔背景与开发集统计区间](detection_evidence.png)']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    plan=['# 执行范围与实际状态','',f'更新时间：{stamp()}','',f'状态：{state["stage"]}','',
        ('旧P波任务的进程、窗口上限、断点和依赖关系保持原样。48 mm导波单会话串行、最多2核；仅在连续无用户输入达到'
         +str(resource_policy.get('require_user_idle_minutes',0))+'分钟且物理内存满足门限时启动。与旧P波并行仅限'
         +', '.join(resource_policy.get('overlap_stages',['各阶段']))+'；其余阶段等待旧队列完成。冻结的科学判据不变。'
         if resource_policy.get('allow_concurrent_with_protected_p') else
         '旧P波任务的进程、窗口上限、断点和依赖关系保持原样。新COMSOL作业等待旧重计算释放资源后串行执行，最多2核；轻量审计和几何生成不启动JVM。'),'',
        '|阶段|范围|实际状态|', '|---|---|---|',
        '|P0|原始双入射矩阵重建、兼容性、12样本表、误差证据|已完成|',
        f'|P1|{len(models)}个几何（含参考）、初始base/fine/absorber变体，401/801/1601端口提取|{stages["p1"]["status"]}|',
        f'|P2|至多4几何；17频点、网格/吸收区、双模态瞬态及半步；280→最多560 μs|{stages["p2"]["status"]}|',
        f'|P3|冻结30种子，48 mm、200 kHz，至少成对网格，固定分母30|{stages["p3"]["status"]}|',
        '|P4|真实噪声、传感器、材料及实验清单|已列清单，实验待数据|','',
        f'P1初始几何变体上限 {len(models)*3}，双模态求解上限 {len(models)*6}（未扣除可复用结果）。P3至少60个孔群网格变体、120个中心频率模态求解，匹配参考优先复用。开始各批次时另记录实测首例耗时，不将此前P波耗时当作导波FD耗时。','',
        '## P1模型与来源','', '|模型|类型|几何来源|', '|---|---|---|']
    for m in models:plan.append(f'|{m["case_id"]}|{m["kind"]}|{m.get("geometry_path") or "匹配无孔矩形"}|')
    plan+=['','代表孔群选择：'+', '.join(map(str,p0['representative_seeds']))+'。中间值按排序后的较低中位秩取值，规则已冻结。',
        '每个失败设置最多追加两级针对性细化；失败后写 needs_adjustment、保留全部版本并阻止依赖阶段。u_num为有限控制差异的保守估计，不是严格数学界。',
        '运行入口：项目 .venv 执行 detectability_200kHz_L48.py audit；ASCII COMSOL Python 执行 detectability_200kHz_L48.py run。原脚本不被修改。']
    (OUT/'PLAN_AND_STATUS.md').write_text('\n'.join(plan)+'\n',encoding='utf-8')
    plots(p0,rows)
    if stages['p3'].get('statistics'):independent_plots(stages['p3'])

def independent_plots(gate):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4),sharey=True)
    for i,(ax,mode) in enumerate(zip(axes,MODES)):
        rows=gate['completed'];x=[v['seed'] for v in rows];y=[100*v['loss'][i] for v in rows];u=[100*v['u_num'][i] for v in rows]
        ax.errorbar(x,y,yerr=u,fmt='o',capsize=2,label='单样本损失及有限数值敏感性估计')
        ax.axhline(.1,color='red',ls='--',label='冻结阈值0.1%');ax.set(title=mode,xlabel='独立 seed');ax.grid(alpha=.25);ax.legend(fontsize=8)
    axes[0].set_ylabel('透射幅值下降（%）');fig.tight_layout();fig.savefig(OUT/'independent_thresholds.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4))
    for i,v in enumerate(gate['statistics']):
        point=v['successes']/30;lower=v['one_sided95_lower'];ax.errorbar(i,point*100,yerr=[[100*(point-lower)],[0]],fmt='o',capsize=7)
        ax.text(i+.08,point*100,f'{v["successes"]}/30; L={100*lower:.2f}%')
    ax.axhline(90,color='red',ls='--');ax.set(xticks=[0,1],xticklabels=MODES,xlim=(-.3,1.65),ylim=(0,105),ylabel='保守成功率及单侧95%下限（%）',title='固定30样本：数值越阈概率，非实验POD')
    fig.tight_layout();fig.savefig(OUT/'independent_probability.png',dpi=160);plt.close(fig)

def plots(p0,rows):
    os.environ['MPLCONFIGDIR']=str(OUT/'_mplconfig')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),sharey=True)
    for ax,mode in zip(axes,MODES):
        selected=[r for r in rows if r['category']=='development_group' and r['mode']==mode]
        ax.plot([r['seed'] for r in selected],[100*r['loss'] for r in selected],'o-',label='逐样本有符号下降')
        measured=[r for r in selected if 'measured_mesh_difference' in r]
        for r in measured:ax.errorbar(r['seed'],100*r['loss'],yerr=100*r['measured_mesh_difference'],fmt='none',capsize=5,color='black',label='仅该例已测网格差')
        ax.axhline(100*TAU,color='red',ls='--',label='固定主阈值 0.1%')
        ax.set(title=mode,xlabel='开发 seed（非独立验证集）');ax.grid(alpha=.25);ax.legend(fontsize=8)
    axes[0].set_ylabel('透射幅值下降（%）')
    fig.suptitle('48 mm / 200 kHz：总 u_num 尚未齐备，不能给出保守检出结论')
    fig.tight_layout();fig.savefig(OUT/'development_thresholds.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.3))
    scores=[100*s['strongest_legacy_single_score'] for s in p0['development_statistics']]
    axes[0].bar(MODES,scores,color=['#2876b5','#d58622']);axes[0].axhline(.1,color='red',ls='--')
    axes[0].set(title='旧单孔最大 score（原端口，仅作背景）',ylabel='幅值变化（%）')
    for i,s in enumerate(p0['development_statistics']):
        point=s['raw_score_successes']/12;lower=s['one_sided95_lower_descriptive']
        axes[1].errorbar(i,100*point,yerr=[[100*(point-lower)],[0]],fmt='o',capsize=6)
        axes[1].text(i+.08,100*point,f'{s["raw_score_successes"]}/12')
    axes[1].axhline(90,color='red',ls='--',label='预设独立验证下限目标90%')
    axes[1].set(xticks=[0,1],xticklabels=MODES,xlim=(-.3,1.4),ylim=(0,103),ylabel='原始越阈比例及单侧95%下限（%）',title='开发数据描述性区间（不是独立验证结论）')
    axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(axis='y',alpha=.25)
    fig.tight_layout();fig.savefig(OUT/'detection_evidence.png',dpi=160);plt.close(fig)
