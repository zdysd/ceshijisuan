"""Resume numerical validation before any length/seed attenuation sweep.

Representative waveforms/movies use COMSOL's transient solver. Only after
transient/frequency consistency and spatial/time controls pass, center-frequency
statistics reuse the validated linear frequency-domain operator.
"""
import argparse,hashlib,json,time,sys
from pathlib import Path
from types import SimpleNamespace
import run_transient as rt
import numpy as np
from guided_statistics import attenuation_interval,coherent_interval
from cfrp_pores.attenuation import fit_attenuation


def arguments():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['convergence','ensemble','report','all'])
    parser.add_argument('--output',type=Path,default=rt.ROOT/'output/transient_converged')
    parser.add_argument('--cores',type=int,default=4)
    parser.add_argument('--samples',type=int,default=12)
    parser.add_argument('--profile',choices=['refined_v4','refined_v5','refined_v6'],default='refined_v4')
    parser.add_argument('--lengths-mm',type=float,nargs='+',default=[8.,16.,24.])
    parser.add_argument('--wait-for-jobs',type=Path)
    args=parser.parse_args()
    if args.samples<6 or len(set(args.lengths_mm))<3:parser.error('Use >=6 seeds and >=3 distinct lengths.')
    return args


def wait_for_existing_jobs(root,path):
    """Wait on actual Windows process handles, avoiding duplicate solves and PID reuse."""
    if path is None:return
    import ctypes
    from ctypes import wintypes
    from datetime import datetime
    kernel=ctypes.windll.kernel32;kernel.OpenProcess.restype=ctypes.c_void_p
    kernel.WaitForSingleObject.argtypes=[ctypes.c_void_p,wintypes.DWORD]
    kernel.CloseHandle.argtypes=[ctypes.c_void_p]
    kernel.GetProcessTimes.argtypes=[ctypes.c_void_p]+[ctypes.POINTER(wintypes.FILETIME)]*4
    jobs=json.loads(path.read_text(encoding='utf-8-sig'));pending={}
    try:
        for job in jobs:
            handle=kernel.OpenProcess(0x100000|0x1000,0,int(job['pid']))
            if not handle:continue
            created,exited,cpu,user=[wintypes.FILETIME() for _ in range(4)]
            if not kernel.GetProcessTimes(handle,*[ctypes.byref(v) for v in [created,exited,cpu,user]]):
                kernel.CloseHandle(handle);raise RuntimeError('Cannot verify an existing solver process.')
            actual=((created.dwHighDateTime<<32)|created.dwLowDateTime)/1e7-11644473600
            expected=datetime.fromisoformat(job['created_utc'].replace('Z','+00:00')).timestamp()
            if abs(actual-expected)>1:
                kernel.CloseHandle(handle);continue
            pending[int(job['pid'])]=handle
        while pending:
            for pid,handle in list(pending.items()):
                if kernel.WaitForSingleObject(handle,0)==0:
                    kernel.CloseHandle(handle);del pending[pid]
            if pending:
                rt.write_json(root/'pipeline_progress.json',dict(stage='waiting_for_existing_solves',pending_pids=list(pending)))
                time.sleep(5)
    finally:
        for handle in pending.values():kernel.CloseHandle(handle)


def case_args(cores,profile='refined_v4'):
    return SimpleNamespace(profile=profile,diameter=80.,porosity=.015,seed=101,length_mm=8.,
        bulk_frequency=5e6,guided_frequency=2e5,duration_factor=1.,frames=False,cores=cores)


def summarize(root,render=False):
    rows,results=rt.report(root,render);gate=rt.validation(root,rows,results)
    rt.statistics(root,rows);rt.append_validation_report(root,gate)
    return gate


def converge(client,root,args):
    config=case_args(args.cores,args.profile)
    for kind in ['guided','bulk']:
        rt.write_json(root/'pipeline_progress.json',dict(stage='spatial_and_spectral_checks',kind=kind,profile=args.profile))
        rt.spatial_checks(client,root,kind,config)
        rt.spatial_checks(client,root,kind,config,band=True)
    controls=[]
    for kind in ['guided','bulk']:
        base=root/'spatial_checks'/kind
        def values(variant,case):return np.load(base/variant/case/'spectra_band.npz')['T']
        for variant in (['mesh','absorber'] if kind=='guided' else ['mesh','absorber','width']):
            for case in ['reference','group']:
                a=np.diagonal(values('base',case),axis1=1,axis2=2);b=np.diagonal(values(variant,case),axis1=1,axis2=2)
                if case=='group':
                    a=a/np.diagonal(values('base','reference'),axis1=1,axis2=2)
                    b=b/np.diagonal(values(variant,'reference'),axis1=1,axis2=2)
                magnitude=float(np.max(abs(abs(a)-abs(b))));phase=float(np.max(abs(np.angle(a/b)))*180/np.pi)
                relative=float(np.max(abs(abs(a)-abs(b))/np.maximum(np.maximum(abs(1-abs(a)),abs(1-abs(b))),1e-12)))
                controls.append(dict(kind=kind,variant=variant,case=case,magnitude_delta=magnitude,phase_error_deg=phase,
                    relative_loss_change=relative,passed=magnitude<=2e-4 and phase<=1 and (case=='reference' or relative<=.05)))
    rt.write_json(root/'spatial_band_convergence.json',controls)
    from refresh_transient_validation import main as refresh_snapshot
    refresh_snapshot(root)
    if not all(c['passed'] for c in controls):raise RuntimeError('Full-band spatial controls failed; transient/statistical expansion is deferred.')
    for kind in ['guided','bulk']:
        for variant in ['base','time']:
            for case in ['reference','group']:
                rt.write_json(root/'pipeline_progress.json',dict(stage='transient_validation',kind=kind,variant=variant,case=case,profile=args.profile))
                # Keep the finer bulk trajectory too, so a later time-step
                # refinement never requires replaying a solved case for movies.
                config.frames=case=='group' and (variant=='base' or kind=='bulk')
                rt.run_case(client,root,kind,case,8.,101,variant,config)
    gate=summarize(root)
    refresh_snapshot(root)
    rt.write_json(root/'pipeline_progress.json',dict(stage='validation_complete',passed=gate['ready_for_ensemble']))
    if not gate['ready_for_ensemble']:
        raise RuntimeError('Numerical gates failed. No length/random-seed sweep has been started; inspect validation.json.')


def center_case(client,root,kind,length_mm,seed,case,variant='base',profile='refined_v4'):
    cfg=json.loads((rt.ROOT/'config.json').read_text());C,rho=rt.validate_material(cfg)
    f0=5e6 if kind=='bulk' else 2e5
    s=rt.layout(kind,length_mm,f0,C,rho,variant,profile)
    plate=rt.Plate(s['width'],C[np.ix_([0,1,5],[0,1,5])],rho,s['safe_elements'])
    pores,metadata=rt.grouped_pores(80.,length_mm,seed,.015)
    if case=='reference':pores=[]
    folder=root/'ensemble'/kind/variant/f'L{length_mm:g}'/f'seed_{seed}'/case
    folder.mkdir(parents=True,exist_ok=True)
    signature=dict(material=cfg,layout={k:v for k,v in s.items() if k not in ['dt','output_dt','end_time','block_duration','adaptive_window','minimum_window']},
        seed=seed,case=case,geometry_hash=hashlib.sha256(b''.join(p.wkb for p in pores)).hexdigest())
    settings=folder/'settings.json'
    if settings.exists() and json.loads(settings.read_text())!=signature:raise RuntimeError('Ensemble cache settings differ.')
    rt.write_json(settings,signature)
    if (folder/'spectra.npz').exists():return dict(np.load(folder/'spectra.npz'))
    model=None
    try:
        for mode in rt.MODES[kind]:
            if (folder/f'{mode}_harmonic.npz').exists():continue
            if model is None:
                minimum=min(s['mesh']/30,min((np.sqrt(4*p.area/np.pi)*1e-6/(60 if variant=='mesh' else 30) for p in pores),default=s['mesh']/30))
                model=rt.create_2d(client,'ensemble_'+mode,cfg,pores,size_m=(s['length'],s['width']),pore_offset_m=(s['start'],0),
                    mesh_m=s['mesh'],mesh_min_m=minimum,guided=kind=='guided',lateral_periodic=kind=='bulk',
                    sponge_start_m=s['sponge'],repair_tolerance_m=1e-9,**rt.builder_options(s))
                if int(model.java.component('comp1').geom('geom1').getNFiniteVoids())!=len(pores):raise RuntimeError('Ensemble CAD lost cavities.')
                measure=model.java.component('comp1').geom('geom1').measureFinal()
                measure.selection().geom('geom1',2);measure.selection().all()
                actual=s['length']*s['width']-float(measure.getArea());target=sum(p.area for p in pores)*1e-12
                if abs(actual-target)/(s['roi_length']*s['width'])>cfg['geometry_porosity_absolute_tolerance']:
                    raise RuntimeError('Ensemble CAD porosity differs from the common pore definition.')
                rt.write_json(folder/'geometry_check.json',dict(actual_porosity=actual/(s['roi_length']*s['width']),target_porosity=target/(s['roi_length']*s['width'])))
                rt.write_json(folder/'boundary_selections.json',model._boundary_selection_checks)
            else:
                j=model.java;j.component('comp1').physics('solid').feature().remove('source')
                for name in ['ur','ui','vr','vi']:
                    if name in list(j.component('comp1').func().tags()):j.component('comp1').func().remove(name)
                for tag in list(j.study().tags()):j.study().remove(tag)
                for tag in list(j.sol().tags()):j.sol().remove(tag)
                for tag in list(j.result().dataset().tags()):j.result().dataset().remove(tag)
            print('CENTER CASE',kind,variant,length_mm,seed,case,mode,flush=True)
            rt.source_and_study(model,s,plate,mode);rt.frequency_checkpoint(model,s,folder,mode)
    finally:
        if model is not None:client.remove(model)
    histories=[dict(np.load(folder/f'{m}_harmonic.npz')) for m in rt.MODES[kind]]
    bank=rt.modal_bank(plate,[f0],histories[0]['y']) if kind=='guided' else None
    result=rt.analyze_ports(histories,[f0],kind,(C,rho),bank,spectra_override=[h['values'][None] for h in histories])
    np.savez_compressed(folder/'spectra.npz',**result)
    rt.write_json(folder/'metadata.json',dict(pore_metadata=metadata,method='Center-frequency COMSOL Solid Mechanics; transient operator validated separately.'))
    return result


def normalized(group,reference):
    return np.diag(group['T'][0])/np.diag(reference['T'][0])


def ensemble(client,root,args):
    gate=json.loads((root/'validation.json').read_text())
    if not gate['ready_for_ensemble']:raise RuntimeError('Both wave families must pass numerical validation before ensemble work.')
    for kind in ['guided','bulk']:
        record=json.loads((root/kind/'base/L8/seed_101/group/run.json').read_text())
        expected=args.profile if kind=='bulk' else 'refined'
        if record['status']!='complete' or record['signature']['layout'].get('profile')!=expected:
            raise RuntimeError('The requested numerical profile has not passed transient validation.')
    lengths=sorted(set(args.lengths_mm));length_checks=[];alpha_checks=[]
    for kind in ['guided','bulk']:
        calibration={'base':[],'mesh':[]}
        for length in lengths:
            rt.write_json(root/'pipeline_progress.json',dict(stage='length_numerical_validation',kind=kind,length_mm=length))
            values={}
            for variant in ['base','mesh']:
                ref=center_case(client,root,kind,length,101,'reference',variant,profile=args.profile)
                group=center_case(client,root,kind,length,101,'group',variant,profile=args.profile)
                values[variant]=normalized(group,ref)
                calibration[variant].append(values[variant])
            a,b=values['base'],values['mesh']
            magnitude=abs(abs(a)-abs(b));phase=abs(np.angle(a/b))*180/np.pi
            relative=magnitude/np.maximum(np.maximum(abs(1-abs(a)),abs(1-abs(b))),1e-12)
            check=dict(kind=kind,length_mm=length,magnitude_delta=magnitude.tolist(),phase_error_deg=phase.tolist(),relative_loss_change=relative.tolist(),
                passed=bool(max(magnitude)<=2e-4 and max(phase)<=1 and max(relative)<=.05))
            length_checks.append(check);rt.write_json(root/'length_numerical_checks.json',length_checks)
            if not check['passed']:raise RuntimeError('A longer specimen failed its mesh check; random expansion has not started.')
        for index,mode in enumerate(rt.MODES[kind]):
            fits={variant:fit_attenuation(np.array(lengths)*1e-3,np.array(signals)[:,index]) for variant,signals in calibration.items()}
            a,b=fits['base']['alpha_np_m'],fits['mesh']['alpha_np_m']
            relative=abs(a-b)/max(abs(a),abs(b),1e-12)
            check=dict(kind=kind,mode=mode,representative_seed=101,base_fit=fits['base'],mesh_fit=fits['mesh'],relative_alpha_change=relative,passed=relative<=.05)
            alpha_checks.append(check);rt.write_json(root/'length_alpha_convergence.json',alpha_checks)
            if not check['passed']:raise RuntimeError('The representative length-fit attenuation failed mesh convergence; random expansion has not started.')
    for kind in ['guided','bulk']:
        for length in lengths:
            reference=center_case(client,root,kind,length,101,'reference',profile=args.profile)
            for seed in range(101,101+args.samples):
                rt.write_json(root/'pipeline_progress.json',dict(stage='random_seed_expansion',kind=kind,length_mm=length,seed=seed,samples=args.samples))
                center_case(client,root,kind,length,seed,'group',profile=args.profile)
    return ensemble_report(root,args)


def ensemble_report(root,args):
    rows=[];lengths=sorted(set(args.lengths_mm));seeds=list(range(101,101+args.samples))
    for kind in ['guided','bulk']:
        data=[]
        for seed in seeds:
            sample=[]
            for length in lengths:
                folder=root/'ensemble'/kind/'base'/f'L{length:g}'
                group=dict(np.load(folder/f'seed_{seed}/group/spectra.npz'))
                reference=dict(np.load(folder/'seed_101/reference/spectra.npz'))
                sample.append(normalized(group,reference))
            data.append(sample)
        data=np.array(data);np.savez_compressed(root/f'{kind}_ensemble_signals.npz',signals=data,seeds=seeds,lengths_m=np.array(lengths)*1e-3)
        for index,mode in enumerate(rt.MODES[kind]):
            signals=data[:,:,index];item=attenuation_interval(np.array(lengths)*1e-3,signals)
            first=root/'ensemble'/kind/'base'/f'L{lengths[0]:g}'/'seed_101/group/spectra.npz'
            with np.load(first) as spectral:wavelength=float(spectral['wavelengths'][0,index])
            half=attenuation_interval(np.array(lengths)*1e-3,signals[:len(seeds)//2])
            alpha=item['fit']['alpha_np_m'];prior=half['fit']['alpha_np_m']
            change=abs(alpha-prior)/max(abs(alpha),abs(prior),1e-30)
            item.update(kind=kind,mode=mode,seeds=seeds,method='Center-frequency ensemble; verified against transient response before expansion.',
                frequency_hz=5e6 if kind=='bulk' else 2e5,material_constant_validated=False,
                host_wavelength_m=wavelength,nominal_pore_diameter_m=80e-6,nominal_diameter_to_host_wavelength=80e-6/wavelength,
                sample_count_relative_change=change,sample_stability_passed=change<=.05,
                alpha_db_m=alpha*8.685889638,ci95_alpha_db_m=(np.array(item['ci95_alpha_np_m'])*8.685889638).tolist(),
                by_length=[dict(length_m=length*.001,**coherent_interval(signals[:,i])) for i,length in enumerate(lengths)])
            rows.append(item)
    rt.write_json(root/'attenuation_center_frequency.json',rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,(kind,mode) in zip(axes,[('bulk','P'),('guided','S0'),('guided','A0')]):
        row=next(v for v in rows if v['kind']==kind and v['mode']==mode)
        data=np.load(root/f'{kind}_ensemble_signals.npz');index=rt.MODES[kind].index(mode)
        values=abs(data['signals'][:,:,index].mean(axis=0));x=data['lengths_m']
        ax.plot(x*1000,np.log(values),'o',label='Coherent mean')
        ax.plot(x*1000,row['fit']['log_amplitude_intercept']-row['fit']['alpha_np_m']*x,label='Fit with intercept')
        ax.set(title=f'{mode}: alpha={row["fit"]["alpha_np_m"]:.4g} Np/m',xlabel='Porous length (mm)',ylabel='ln |mean T|')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle(f'{args.samples} matched seed families; illustrative material; conditional on 2D geometry')
    fig.tight_layout();fig.savefig(root/'attenuation_length_fits.png',dpi=160);plt.close(fig)
    lines=['# 更新后的长度与随机样本衰减结果','',
        '先通过瞬态时间步、观察窗、频谱一致性和空间离散验证，再以同一线性固体力学频域算子完成中心频率的长度／随机样本扫描。动画来自独立瞬态求解，频域结果未冒充瞬态时间帧。',
        '本结果仅针对示例材料和二维孔腔。体波为 1 mm 横向周期单元；重复同一单元的宽度检查属于数值验证，不证明无限随机体的代表体积已收敛。吸收层不属于拟合孔隙段。','',
        '|模态|频率 MHz|无孔基体波长 mm|名义孔径／波长|', '|---|---:|---:|---:|']
    for row in rows:
        lines.append(f'|{row["mode"]}|{row["frequency_hz"]/1e6:g}|{row["host_wavelength_m"]*1000:.6g}|{row["nominal_diameter_to_host_wavelength"]:.6g}|')
    lines+=['','孔径为名义 80 μm 组，孔隙率 1.5%，孔段长度为 '+', '.join(f'{v:g}' for v in lengths)+' mm。频率及孔径／波长比不同，不能仅按幅值下降直接排序体波与导波的普遍灵敏度。','',
        '|模态|频率 MHz|样本数|相干 alpha Np/m|95% 区间 Np/m|典型单样本 alpha Np/m|样本加倍变化|指数拟合 R²|',
        '|---|---:|---:|---:|---|---:|---:|---:|']
    for row in rows:
        ci=row['ci95_alpha_np_m'];lines.append(f'|{row["mode"]}|{row["frequency_hz"]/1e6:g}|{row["samples"]}|{row["fit"]["alpha_np_m"]:.6g}|[{ci[0]:.6g}, {ci[1]:.6g}]|{row["typical_fit"]["alpha_np_m"]:.6g}|{100*row["sample_count_relative_change"]:.3g}%|{row["fit"]["r_squared"]:.5g}|')
    lines+=['','采用 ln|T|=b-alpha L，保留截距；Bootstrap 重采样单位为跨所有长度的完整种子族。区间反映当前有限样本的不确定性，不涵盖材料参数未知与二维建模误差。',
        '随机扩展前，还检查了每个长度的代表样本，以及三长度拟合斜率的网格收敛；记录见 length_numerical_checks.json 和 length_alpha_convergence.json。',
        '样本加倍变化超过 5% 时，仅报告有限样本估计；置信区间跨零时不宣称识别出正衰减；指数拟合未通过时不把斜率解释为统一的长度无关衰减常数。',
        '相干衰减与典型单样本幅值衰减分列，均包含散射／反射／转换效应，不等于材料吸收。','',
        f'![长度拟合]({(root/"attenuation_length_fits.png").resolve().as_posix()})']
    (root/'ATTENUATION_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return rows


def main():
    args=arguments();root=args.output.resolve();root.mkdir(parents=True,exist_ok=True);rt.configure_scratch(root)
    client=None
    try:
        wait_for_existing_jobs(root,args.wait_for_jobs)
        if args.stage!='report':
            rt.write_json(root/'pipeline_progress.json',dict(stage='starting',requested_stage=args.stage,profile=args.profile))
            client=rt.start_session(root,args.cores)
        if args.stage in ['all','convergence']:converge(client,root,args)
        if args.stage in ['all','ensemble']:ensemble(client,root,args)
        if args.stage in ['all','report']:
            summarize(root,True)
            if (root/'attenuation_center_frequency.json').exists():ensemble_report(root,args)
            if client is not None:rt.audit_saved_models(client,root)
            from refresh_transient_validation import main as refresh_snapshot
            refresh_snapshot(root)
        statistics_path=root/'attenuation_center_frequency.json'
        statistics=json.loads(statistics_path.read_text()) if statistics_path.exists() else None
        rt.write_json(root/'pipeline_progress.json',dict(stage='complete',requested_stage=args.stage,
            statistical_expansion_completed=statistics is not None,
            sample_stability_all_passed=all(row['sample_stability_passed'] for row in statistics) if statistics else None,
            attenuation_fits_all_accepted=all(row['exponential_fit_passed'] and row['positive_attenuation_identified'] for row in statistics) if statistics else None))
    except Exception as error:
        rt.write_json(root/'pipeline_progress.json',dict(stage='needs_adjustment',error=str(error),requested_stage=args.stage));raise


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');sys.stderr.reconfigure(encoding='utf-8');main()
