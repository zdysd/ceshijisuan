"""Resume bounded P window/time refinements without altering accepted Lamb data."""
import hashlib,json,shutil
from copy import deepcopy
from pathlib import Path
import run_transient as rt
import run_refined_pipeline as pipeline
from cfrp_pores.transient import bulk_adjustment,next_bulk_window,bulk_window_limit


PAIRS=[(v,c) for v in ['base','time'] for c in ['reference','group']]


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def checked(path):
    path=Path(path).resolve()
    if not path.is_relative_to(rt.ROOT.resolve()):raise RuntimeError('Bulk trial path is outside the workspace.')
    return path


def write_limit_diagnosis(root,gate,rows,capped):
    """Report measured failures; physical causes remain hypotheses to verify."""
    selected=[r for r in rows if r['kind']=='bulk' and r['variant'] in ['base','time']]
    checks=[c for c in gate['checks'] if c['kind']=='bulk']
    rt.write_json(root/'bulk_limit_diagnosis.json',dict(capped_cases=capped,
        controls=checks,windows=selected,thresholds=gate['thresholds'],
        physical_root_cause_confirmed=False,additional_window_solves_allowed=False))
    limits=gate['thresholds']
    lines=['# P 波窗口上限后的失败诊断','',
        '参考上限 50 μs，孔群上限 100 μs。下列指标来自已保存的完整本轮对照；停止增加观察时间，不放宽验收门限。','',
        '|算例|步长对照|实际记录 μs|窗口复透射变化|尾部位移平方积分占比|',
        '|---|---|---:|---:|---:|']
    for r in selected:
        lines.append(f'|{r["case"]}|{r["variant"]}|{r["recorded_end_time_s"]*1e6:.6g}|{r["window_delta_max"]:.6g}|{r["tail_energy_fraction"]:.6g}|')
    lines+=['',f'正式窗口门限 {limits["absolute_amplitude"]:g}；尾部门限 {limits["tail_energy"]:g}。尾部指标不是材料能量吸收率。','']
    for c in checks:
        lines.append(f'## {c["case"]} 的完整对照证据')
        lines+=['','```json',json.dumps(dict(passed=c['passed'],frequency_band_consistency=c.get('frequency_band_consistency'),
            half_step=c.get('differences',{}).get('time'),arrival_checks=c.get('arrival_checks')),ensure_ascii=False,indent=2),'```','']
    lines+=['## 原因判断及调整方案','',
        '已证实的失败项目以表格和完整对照为准。仅凭窗口不收敛，不能确认是边界回波、孔群散射长尾或信号提取误差。以下是需要用现有记录区分的原因，尚不是已确证的物理结论。','',
        '1. 先用已保存的截面位移和牵引应力检查各时窗的入射／反射／透射分量与尾波到达顺序，结合传播距离与群速度判断是否存在端部回波污染；若证据支持，调整吸收区或端部位置，并做单因素反射对照。',
        '2. 对相同完整记录核对频带、输出采样及窗口端点，分别检查原始复透射和无孔归一化结果；若问题来自提取，修正提取流程后重新分析已有记录，避免重复有限元求解。',
        '3. 若窗口稳定而半时间步或频域一致性失败，优先检查时间离散的幅值／相位误差；根据实际失败频点提出定向时间步调整，仍保持各自窗口上限。',
        '4. 若前述检查排除了回波和提取问题，且孔群尾波仍显著，则本轮只能报告该窗口内的响应，不能宣布已经获得收敛的透射损失或长度衰减常数。','',
        'A0／S0 已通过结果保留。P 未通过时不展开长度／随机样本统计。可信范围限于示例材料、二维孔腔和当前周期单元。',
        '',f'[全部失败指标与阈值]({(root/"bulk_limit_diagnosis.json").resolve().as_posix()})']
    (root/'BULK_LIMIT_DIAGNOSIS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def prepare_case(source,destination,action):
    """Copy immutable completed data; only extend its window or relabel its role."""
    source=checked(source);destination=checked(destination)
    record=deepcopy(read(source/'run.json'))
    if record['status']!='complete':raise RuntimeError('Bulk refinement source is incomplete.')
    s=record['signature']['layout']
    if s['kind']!='bulk' or s.get('profile')!='refined_v6':raise RuntimeError('Unsupported bulk refinement source.')
    original=deepcopy(record['signature']);fingerprint=record['fingerprint']
    if action=='promote_time':
        if s['variant']!='time':raise RuntimeError('Only an independently solved half-step case can be promoted.')
        s.update(variant='base',time_step_refinement=2*s.get('time_step_refinement',1.))
    elif action=='extend_window':
        target=next_bulk_window(s['end_time'],bulk_window_limit(record['signature']['case']))
        if target is None:raise RuntimeError('P observation window has reached its case-specific limit.')
        s.update(end_time=target,window_extension_factor=2*s.get('window_extension_factor',1.))
        record.update(status='checkpoint_fork',runs=[])
        for key in ['active_mode','active_recovery_directory','error']:record.pop(key,None)
    elif action!='copy':raise ValueError(action)
    record['fingerprint']=hashlib.sha256(json.dumps(record['signature'],sort_keys=True).encode()).hexdigest()
    provenance=dict(action=action,source=str(source),source_fingerprint=fingerprint,
                    target_initial_fingerprint=record['fingerprint'],original_signature=original)
    if (destination/'run.json').exists():
        previous=read(destination/'bulk_refinement_provenance.json')
        if previous!=provenance:raise RuntimeError('Existing bulk trial was created from different settings.')
        return
    ignore=(shutil.ignore_patterns('run.json','P_history.npz','P_solve.json','P_solve_progress.json',
                                  'P_fields.npz','P_block_validation.json','spectra.npz')
            if action=='extend_window' else shutil.ignore_patterns('run.json'))
    shutil.copytree(source,destination,dirs_exist_ok=True,ignore=ignore)
    if action=='extend_window':
        progress=destination/'blocks/P/progress.json';state=read(progress)
        state.update(window_converged=False,stable_windows=0)
        rt.write_json(progress,state)
    rt.write_json(destination/'bulk_refinement_provenance.json',provenance)
    rt.write_json(destination/'run.json',record)


def commit_trial(root,plan):
    """Recoverable directory swap: the old complete case tree stays archived."""
    active=checked(root/'bulk');trial=checked(Path(plan['trial'])/'bulk')
    archive=checked(root/'bulk_history'/f'attempt_{plan["attempt"]:03d}')
    if trial.exists():
        if active.exists():
            if archive.exists():raise RuntimeError('Both active and archived bulk trees exist before trial installation.')
            archive.parent.mkdir(parents=True,exist_ok=True);active.rename(archive)
        trial.rename(active)
    for v,c in PAIRS:
        case=active/v/f'L8/seed_101/{c}';actual=read(case/'run.json')
        if actual['status']!='complete' or actual['fingerprint']!=plan['fingerprints'][f'{v}/{c}']:
            raise RuntimeError('Installed bulk trial does not match its complete-case manifest.')
        provenance=case/'bulk_refinement_provenance.json'
        if provenance.exists():
            details=read(provenance);source=Path(details['source'])
            if source.is_relative_to(active):
                archived=archive/source.relative_to(active)
                if read(archived/'run.json')['fingerprint']!=details['source_fingerprint']:
                    raise RuntimeError('Archived refinement source does not match its provenance.')
                details['archived_source']=str(archived);rt.write_json(provenance,details)


def refine(client,root,status,max_adjustments=3,cores=2):
    """At most three measured adjustments; unexplained failures remain explicit."""
    root=Path(root);state_path=root/'bulk_refinement_state.json'
    state=read(state_path) if state_path.exists() else dict(attempt=0,stage='evaluate')
    while True:
        if state['stage']=='install':
            commit_trial(root,state)
            state.update(stage='evaluate');rt.write_json(state_path,state)
        if state['stage']=='evaluate':
            rows,results=rt.report(root);gate=rt.validation(root,rows,results)
            rt.statistics(root,rows);rt.append_validation_report(root,gate)
            from refresh_transient_validation import main as refresh
            refresh(root)
            decision=bulk_adjustment(gate,rows)
            if decision['action']=='accepted':return gate
            if decision['action']=='extend_window':
                capped=[pair for pair in decision['cases'] if next_bulk_window(
                    read(root/'bulk'/pair[0]/f'L8/seed_101/{pair[1]}/run.json')['signature']['layout']['end_time'],bulk_window_limit(pair[1])) is None]
                if capped:
                    write_limit_diagnosis(root,gate,rows,capped)
                    rt.write_json(root/'bulk_window_stop.json',dict(
                        reason='Case-specific hard limit reached; diagnose the failed complete controls before further changes.',
                        failed_cases=decision['cases'],capped_cases=capped,window_limits_s={c:bulk_window_limit(c) for c in ['reference','group']},
                        automatic_window_extension_allowed=True,extension_stopped_at_limit=True))
                    raise RuntimeError('P window checks failed at the reference 50 us / group 100 us limits; no longer window will be solved. Diagnose the retained complete controls.')
            if state['attempt']>=max_adjustments:
                raise RuntimeError('P controls still fail after three measured window/time adjustments; inspect all retained trials before further changes.')
            number=state['attempt']+1;trial=checked(root/'bulk_trials'/f'attempt_{number:03d}')
            state=dict(attempt=number,stage='prepare',trial=str(trial),decision=decision)
            rt.write_json(state_path,state)
        trial=Path(state['trial']);action=state['decision']['action']
        status('bulk_refinement',attempt=state['attempt'],action=action,trial=str(trial))
        if state['stage']=='prepare':
            shutil.copytree(root/'spatial_checks',trial/'spatial_checks',dirs_exist_ok=True)
            for variant,case in PAIRS:
                destination=trial/'bulk'/variant/f'L8/seed_101/{case}'
                if action=='refine_time':
                    if variant=='base':prepare_case(root/'bulk/time'/f'L8/seed_101/{case}',destination,'promote_time')
                else:
                    operation='extend_window' if [variant,case] in state['decision']['cases'] else 'copy'
                    prepare_case(root/'bulk'/variant/f'L8/seed_101/{case}',destination,operation)
            state.update(stage='solve');rt.write_json(state_path,state)
        for variant,case in PAIRS:
            folder=trial/'bulk'/variant/f'L8/seed_101/{case}'
            record=read(folder/'run.json') if (folder/'run.json').exists() else None
            if record and record['status']=='complete':continue
            config=pipeline.case_args(cores,'refined_v6')
            if record:
                s=record['signature']['layout'];config.frames=record['signature']['frames']
            else:
                s=read(trial/'bulk/base'/f'L8/seed_101/{case}/run.json')['signature']['layout']
                config.frames=case=='group'
            config.bulk_step_factor=s.get('time_step_refinement',1.)
            config.bulk_window_factor=s.get('window_extension_factor',1.)
            config.bulk_window_end_time=s['end_time']
            status('bulk_refinement',attempt=state['attempt'],action=action,variant=variant,case=case,
                   step_factor=config.bulk_step_factor,window_factor=config.bulk_window_factor)
            rt.run_case(client,trial,'bulk',case,8.,101,variant,config)
        state.update(stage='install',fingerprints={f'{v}/{c}':read(trial/'bulk'/v/f'L8/seed_101/{c}/run.json')['fingerprint'] for v,c in PAIRS})
        rt.write_json(state_path,state)
