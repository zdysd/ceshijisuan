"""Durable continuation: native receiver sampling, longer window, then gated statistics.

The existing jobs are left running. Their completed, immutable outputs are reused.
No threshold is relaxed and no statistical sweep starts on a failed numerical gate.
"""
import hashlib,json,msvcrt,os,shutil,time
from datetime import datetime,timezone
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import run_transient as rt
import run_refined_pipeline as pipeline
from watch_refined_results import alive

OLD=rt.ROOT/'output/transient_converged'
PROBE=rt.ROOT/'output/transient_refined/guided_sampling'
OUT=rt.ROOT/'output/transient_verified'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def status(stage,**details):
    value=dict(updated_at=datetime.now(timezone.utc).isoformat(),stage=stage,**details)
    rt.write_json(OUT/'continuation_state.json',value)
    lines=['# 瞬态验证续算状态','',f'更新时间：{value["updated_at"]}',f'当前阶段：{stage}','',
        '本目录使用每个内部时间步输出的 Lamb 波接收历史，观察窗延长到 280 μs，并独立计算时间步减半对照。',
        '体波数据只在原队列完成后导入，仍须通过全部原有门限。任何未通过的检查都会阻止长度和随机样本扩展。',
        '实际位移动画必须来自求解器时间帧；衰减统计不能冒充已完成或已收敛。','',
        '临时和恢复目录：E:/comsol_runtime/transient_matched/。',
        '续算入口：工作区 resume_transient.ps1；交接说明：CONTINUE_TRANSIENT.md。']
    if details:lines+=['','```json',json.dumps(details,ensure_ascii=False,indent=2),'```']
    for name,label in [('REPORT.md','当前波形与数值验收'),('LATEST_VALIDATION.md','最新验证快照'),
                       ('BULK_LIMIT_DIAGNOSIS.md','P 波达到上限后的失败指标与调整方案'),
                       ('ATTENUATION_REPORT.md','衰减与样本区间'),('validation.json','完整验收数据')]:
        path=OUT/name
        if path.exists():lines+=['',f'[{label}]({path.as_posix()})']
    (OUT/'CONTINUATION_STATUS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    result_lines=list(lines)
    if stage=='complete':
        result_lines += ['', '数值验收和长度／随机样本计算已完成。样本稳定性与指数拟合仍以衰减报告中的通过标记和置信区间为准。',
            '可信范围限定为示例材料、二维孔腔及当前周期单元；同模态透射损失包含散射和模态转换，不等于真实 CFRP 的材料吸收。']
    else:result_lines += ['', '全部数值验收及统计扩展尚未完成，未完成结果不解释为已收敛的材料衰减常数。']
    for kind,mode in [('guided','S0'),('guided','A0'),('bulk','P')]:
        folder=OUT/kind/'base/L8/seed_101/group';movie=folder/f'{mode}_displacement.mp4'
        validation=folder/f'{mode}_video_validation.json'
        if movie.exists() and validation.exists() and read(validation).get('full_decode_verified'):
            result_lines += ['',f'![{mode} 实际位移动画]({movie.as_posix()})']
    for name,label in [('latest_transient_comparison.png','导波瞬态与频域对比'),('comparison_summary.png','体波与导波对比'),
                       ('attenuation_length_fits.png','衰减随孔隙段长度的拟合')]:
        path=OUT/name
        if path.exists():result_lines += ['',f'![{label}]({path.as_posix()})']
    temporary=OUT/f'RESULTS.{os.getpid()}.tmp';temporary.write_text('\n'.join(result_lines)+'\n',encoding='utf-8');temporary.replace(OUT/'RESULTS.md')
    print(json.dumps(value,ensure_ascii=False),flush=True)


def wait_for_cases(root,folders,process_file,stage):
    previous=None
    while True:
        pending=[str(p) for p in folders if not (root/p/'run.json').exists() or read(root/p/'run.json').get('status')!='complete']
        if not pending:return
        if pending!=previous:status(stage,pending=pending);previous=pending
        process=read(process_file) if process_file.exists() else {}
        if not alive(process.get('launcher_pid',0),process.get('started')):
            raise RuntimeError(f'Source calculation stopped with incomplete cases. Resume it first: {process_file}')
        time.sleep(30)


def fork_native_prefix(source,destination):
    """Change only the bookkeeping variant and final time; reuse u AND velocity."""
    if (destination/'run.json').exists():return
    record=read(source/'run.json')
    if record['status']!='complete':raise RuntimeError('Cannot fork an unfinished source case.')
    s=record['signature']['layout']
    if s['kind']!='guided' or s['variant']!='sampling' or s['output_dt']!=s['dt']:
        raise RuntimeError('The source must be the completed native-output guided case.')
    for path in [source,destination]:
        if not path.resolve().is_relative_to(rt.ROOT.resolve()):raise RuntimeError('Copy is outside the workspace.')
    destination.mkdir(parents=True,exist_ok=True)
    shutil.copytree(source/'blocks',destination/'blocks',dirs_exist_ok=True)
    for pattern in ['geometry.json','*_harmonic*.npz','*_boundary_selections.json']:
        for path in source.glob(pattern):shutil.copy2(path,destination/path.name)
    provenance=dict(source_case=str(source.resolve()),source_fingerprint=record['fingerprint'],
                    imported_blocks={m:read(source/'blocks'/m/'progress.json')['completed'] for m in rt.MODES['guided']},
                    original_end_time_s=s['end_time'],new_end_time_s=2*s['end_time'],
                    method='Preserve completed endpoints, displacement and velocity; continue after the original final short block.')
    signature=record['signature'];signature['layout']=dict(s,variant='base',end_time=2*s['end_time'])
    record.update(signature=signature,status='checkpoint_fork',runs=[],
                  fingerprint=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest())
    for key in ['active_mode','active_recovery_directory','error']:record.pop(key,None)
    rt.write_json(destination/'continuation_provenance.json',provenance)
    rt.write_json(destination/'run.json',record)


def reuse_verified_movie_fields():
    """A change of receiver output cadence leaves the manually stepped FE trajectory unchanged."""
    src=OLD/'guided/base/L8/seed_101/group';dst=OUT/'guided/base/L8/seed_101/group'
    old=read(src/'run.json')['signature'];new=read(dst/'run.json')['signature']
    ignored={'variant','output_dt','end_time'}
    if old['material']!=new['material'] or old['geometry_hash']!=new['geometry_hash'] or \
       {k:v for k,v in old['layout'].items() if k not in ignored}!={k:v for k,v in new['layout'].items() if k not in ignored}:
        raise RuntimeError('Movie source uses a different physical or internal-step model.')
    C,rho=rt.validate_material(new['material']);stress=2*np.pi*new['layout']['f0']*np.sqrt(rho*C[0,0])
    scale=np.array([1,1,stress,stress])[None,None,:,None];checks=[]
    for mode in rt.MODES['guided']:
        a=dict(np.load(src/f'{mode}_history.npz'));b=dict(np.load(dst/f'{mode}_history.npz'))
        indices=np.rint(a['times']/new['layout']['output_dt']).astype(int)
        np.testing.assert_allclose(a['times'],b['times'][indices],rtol=1e-12,atol=1e-17)
        error=float(np.linalg.norm((a['values']-b['values'][indices])/scale)/np.linalg.norm(a['values']/scale))
        if error>1e-8:raise RuntimeError('Shared histories differ; original movie fields cannot be reused.')
        shutil.copy2(src/f'{mode}_fields.npz',dst/f'{mode}_fields.npz')
        checks.append(dict(mode=mode,shared_history_relative_error=error,source_fields=str(src/f'{mode}_fields.npz'),
                           movie_end_time_s=float(np.load(src/f'{mode}_fields.npz')['times'][-1]),
                           receiver_validation_end_time_s=float(b['times'][-1])))
    rt.write_json(dst/'movie_field_provenance.json',dict(checks=checks,
        note='Actual original FE frames of the identical manually stepped model, verified at every shared receiver sample. Movie shows the first 140 us; receiver validation extends to 280 us.'))


def main():
    OUT.mkdir(parents=True,exist_ok=True);rt.configure_scratch(OUT)
    lock=(OUT/'.continuation.lock').open('a+b');lock.seek(0)
    try:msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
    except OSError:raise RuntimeError('Another continuation worker already owns this output directory.')
    try:
        runtime=OUT/'_runtime';runtime.mkdir(exist_ok=True)
        sources=['continue_transient_validation.py','refine_bulk_transient.py','run_transient.py','run_refined_pipeline.py','cfrp_pores/transient.py','cfrp_pores/comsol.py']
        rt.write_json(runtime/f'continuation_code_{os.getpid()}.json',dict(pid=os.getpid(),
            started=datetime.now(timezone.utc).isoformat(),
            sha256={p:hashlib.sha256((rt.ROOT/p).read_bytes()).hexdigest() for p in sources}))
        if (OUT/'continuation_state.json').exists() and read(OUT/'continuation_state.json')['stage']=='complete':
            print('Already complete; no recomputation.',flush=True);return
        wait_for_cases(PROBE,[Path(f'guided/sampling/L8/seed_101/{c}') for c in ['reference','group']],
                       OLD/'sampling_process.json','waiting_for_sampling_controls')
        check=rt.ROOT/'output/transient_refined/window_extension_validation_bulk_p3.json'
        if not check.exists() or not read(check)['passed']:raise RuntimeError('Independent extended-window restart check has not passed.')
        shutil.copytree(OLD/'spatial_checks',OUT/'spatial_checks',dirs_exist_ok=True)
        config=pipeline.case_args(2,'refined_v6');config.native_guided_output=True;config.duration_factor=2.;config.frames=False
        for case in ['reference','group']:
            fork_native_prefix(PROBE/f'guided/sampling/L8/seed_101/{case}',OUT/f'guided/base/L8/seed_101/{case}')
        guide_cases=[OUT/f'guided/{v}/L8/seed_101/{c}' for v in ['base','time'] for c in ['reference','group']]
        cached_guide=all((p/'run.json').exists() and read(p/'run.json').get('status')=='complete' for p in guide_cases)
        existing_gate=read(OUT/'validation.json') if (OUT/'validation.json').exists() else {}
        accepted_guide=[c for c in existing_gate.get('checks',[]) if c['kind']=='guided' and c['case'] in ['reference','group']]
        guide_already_valid=cached_guide and len(accepted_guide)==2 and all(c['passed'] for c in accepted_guide)
        # Complete cached cases still undergo their signature checks below.
        # Keep the waiting worker free of an idle JVM during the long P solves.
        client=None if cached_guide else rt.start_session(OUT,2)
        for variant in ['base','time']:
            for case in ['reference','group']:
                status('checking_cached_guided_results' if cached_guide else 'guided_native_output_and_long_window',
                       variant=variant,case=case,end_time_us=280)
                rt.run_case(client,OUT,'guided',case,8.,101,variant,config)
        from refresh_transient_validation import main as refresh
        gate=existing_gate if guide_already_valid else pipeline.summarize(OUT)
        if not guide_already_valid:refresh(OUT)
        controls=[c for c in gate['checks'] if c['kind']=='guided' and c['case'] in ['reference','group']]
        if len(controls)!=2 or not all(c['passed'] for c in controls):
            raise RuntimeError('Native-output, 280 us guided controls did not all pass; inspect validation.json. Statistical expansion remains blocked.')
        if not guide_already_valid or not (OUT/'guided/base/L8/seed_101/group/movie_field_provenance.json').exists():
            reuse_verified_movie_fields();pipeline.summarize(OUT,True)
        bulk=[Path(f'bulk/{v}/L8/seed_101/{c}') for v in ['base','time'] for c in ['reference','group']]
        imported=OUT/'bulk_import.json'
        if not imported.exists():
            wait_for_cases(OLD,bulk,OLD/'pipeline_process.json','waiting_for_bulk_transient_controls')
            for folder in bulk:shutil.copytree(OLD/folder,OUT/folder,dirs_exist_ok=True)
            rt.write_json(imported,dict(source=str(OLD),cases=[str(p) for p in bulk],
                source_fingerprints={str(p):read(OLD/p/'run.json')['fingerprint'] for p in bulk}))
        if client is None:
            process=read(OLD/'pipeline_process.json') if (OLD/'pipeline_process.json').exists() else {}
            if alive(process.get('launcher_pid',0),process.get('started')):
                status('waiting_for_source_release',reason='All P cases are saved; wait for the source JVM to finish reporting and release memory.')
                while alive(process.get('launcher_pid',0),process.get('started')):time.sleep(10)
            client=rt.start_session(OUT,4)
        from refine_bulk_transient import refine
        worker_cores=4 if cached_guide else 2
        gate=refine(client,OUT,status,cores=worker_cores)
        gate=pipeline.summarize(OUT,True);refresh(OUT)
        if not gate['ready_for_ensemble']:raise RuntimeError('Combined guided/P numerical acceptance failed; no length/random expansion was started.')
        status('numerical_acceptance_passed_starting_length_and_seed_controls')
        args=SimpleNamespace(profile='refined_v6',cores=worker_cores,samples=12,lengths_mm=[8.,16.,24.])
        rows=pipeline.ensemble(client,OUT,args)
        pipeline.summarize(OUT,True);rt.audit_saved_models(client,OUT);refresh(OUT)
        status('complete',sample_stability_all_passed=all(r['sample_stability_passed'] for r in rows),
               attenuation_fits_all_accepted=all(r['exponential_fit_passed'] and r['positive_attenuation_identified'] for r in rows))
    except Exception as error:
        status('needs_adjustment',error=str(error));raise
    finally:
        lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()


if __name__=='__main__':
    import sys
    sys.stdout.reconfigure(encoding='utf-8');sys.stderr.reconfigure(encoding='utf-8');main()
