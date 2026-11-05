"""Gated 48 mm bandwidth/transient checks and the fixed independent cohort."""
import json,shutil,time
from pathlib import Path
import numpy as np
from scipy.signal import hilbert
import detectability_200kHz_L48 as d
import detectability_worker as w
import run_transient as rt
from cfrp_pores.portability import runtime_root,atomic_savez,validate_npz

def transient_case(client,model,level,half,duration):
    s0=w.settings();protocol=d.read(d.OUT/'protocol.json');freq=np.array(protocol['transient']['frequencies_hz'])
    if duration>560e-6+1e-15:raise RuntimeError('Guided record limit is 560 us')
    dt=protocol['transient']['initial_dt_s']/(2 if half else 1)
    s=dict(kind='guided',profile='detectability_L48_v1',f0=d.F0,fmin=float(freq[0]),fmax=float(freq[-1]),
        roi_length=.048,start=s0['start'],width=.001,sponge=s0['sponge'],length=s0['length'],absorb_length=s0['absorb_length'],
        positions=s0['pre_positions']+s0['post_positions'],mesh=.00015/2**level,dt=dt,output_dt=dt,
        end_time=duration,block_duration=30e-6,port_samples=401,safe_elements=120,damping_factor=3.)
    slowest=min(s0['plate'].at_frequency(f,m)['group_velocity_m_s'] for f in freq for m in d.MODES)
    first_arrival=max(s['positions'])/slowest+5/d.F0
    if duration<first_arrival:raise RuntimeError('Initial record does not contain the slowest target packet')
    folder=d.OUT/'transient'/model['case_id']/f'level_{level}'/('half' if half else 'base')/f'window_{round(duration*1e6)}us'
    folder.mkdir(parents=True,exist_ok=True)
    sig=dict(layout=s,material=s0['config'],geometry_hash=model['geometry_hash'],source='Five-cycle Hann with full quadrature; actual two-launch decomposition',
             source_code_sha256=d.digest(d.ROOT/'run_transient.py'))
    d.freeze(folder/'signature.json',sig)
    if (folder/'complete.json').exists():return analyze_transient(folder)
    if duration>280e-6:
        earlier=folder.parent/'window_280us';provenance=folder/'prefix_provenance.json'
        if not provenance.exists():
            if not (earlier/'complete.json').exists():raise RuntimeError('Only a complete 280 us prefix can be extended')
            old=d.read(earlier/'signature.json');expected=dict(sig,layout=dict(s,end_time=280e-6))
            if old!=expected:raise RuntimeError('Extension would change an already solved physical/numerical setting')
            shutil.copytree(earlier/'blocks',folder/'blocks',dirs_exist_ok=True)
            d.write(provenance,dict(source=str(earlier),signature_sha256=d.digest(earlier/'signature.json'),method='Preserve complete u and UDot; append absolute-time blocks'))
    for mode in d.MODES:
        history_path=folder/f'{mode}_history.npz'
        if history_path.exists():
            validate_npz(history_path,required=('times','values','y','positions'))
            continue
        w.wait_resources('P2');w.state('solving',stage_name='P2_transient',case=model['case_id'],mode=mode,half_step=half,duration_us=duration*1e6)
        model_native=None
        try:
            progress=folder/'blocks'/mode/'progress.json'
            if progress.exists():
                p=d.read(progress);checkpoint=progress.parent/p.get('checkpoint_file','restart.mph')
                if p['completed'] and not checkpoint.exists():raise RuntimeError('Missing committed native checkpoint')
                model_native=client.load(str(checkpoint.resolve()))
                model_native.java.param().set('tend',f'{duration}[s]')
            else:
                model_native,geometry=w.build(client,s0,model,level)
                d.write(folder/'geometry_check.json',geometry)
                rt.source_and_study(model_native,s,s0['plate'],mode)
            actual=str(client.java.getRecoveryDir(model_native.java.tag()))
            if not Path(actual).resolve().is_relative_to(runtime_root()):raise RuntimeError('Recovery escaped configured runtime directory')
            d.write(folder/'progress.json',dict(mode=mode,status='running',recovery_directory=actual))
            started=time.perf_counter();history=rt.solve_in_blocks(model_native,s,folder,mode,False)
            atomic_savez(history_path,**history)
            d.write(folder/f'{mode}_timing.json',dict(seconds=time.perf_counter()-started,actual_end_time_s=float(history['times'][-1])))
        finally:
            if model_native is not None:client.remove(model_native)
    d.write(folder/'complete.json',dict(signature_sha256=d.signature(sig),first_arrival_bound_s=first_arrival,record_complete=True,
        history_sha256={f'{m}_history.npz':d.digest(folder/f'{m}_history.npz') for m in d.MODES}))
    return analyze_transient(folder)

def analyze_transient(folder):
    complete=d.read(folder/'complete.json')
    if complete['signature_sha256']!=d.signature(d.read(folder/'signature.json')):raise RuntimeError('Transient signature mismatch')
    for name,expected in complete.get('history_sha256',{}).items():
        if d.digest(folder/name)!=expected:raise RuntimeError('Transient history changed')
    for mode in d.MODES:validate_npz(folder/f'{mode}_history.npz',required=('times','values','y','positions'))
    s=d.read(folder/'signature.json')['layout'];settings=w.base_settings();histories=[dict(np.load(folder/f'{m}_history.npz')) for m in d.MODES]
    f=np.linspace(s['fmin'],s['fmax'],17);bank=rt.modal_bank(settings['plate'],f,histories[0]['y'])
    result=rt.analyze_ports(histories,f,'guided',(settings['C'],settings['rho']),bank)
    trim=rt.analyze_ports(histories,f,'guided',(settings['C'],settings['rho']),bank,end_fraction=.9)
    tail=max(rt.tail_fraction(h) for h in histories);window=float(abs(result['T']-trim['T']).max())
    arrivals=[]
    for h,mode in zip(histories,d.MODES):
        component=0 if mode=='S0' else 1;envelope=abs(hilbert(h['values'][:,:3,component,-1],axis=0));peaks=h['times'][np.argmax(envelope,axis=0)]
        measured=(h['positions'][2]-h['positions'][0])/(peaks[2]-peaks[0]) if peaks[2]>peaks[0] else None
        expected=settings['plate'].at_frequency(d.F0,mode)['group_velocity_m_s']
        arrivals.append(dict(mode=mode,measured=measured,expected=expected,relative_error=None if measured is None else abs(measured/expected-1)))
    quality=dict(tail=tail,window_delta=window,residual=float(np.max(result['residual'])),right=float(np.max(abs(result['right']))),
        power_balance_error=float(np.max(abs(result['power_balance']-1))),arrival=arrivals,
        actual_end_time_s=float(histories[0]['times'][-1]),window_passed=tail<=.0001 and window<=.0002)
    np.savez_compressed(folder/'spectra.npz',**result);d.write(folder/'quality.json',quality)
    return dict(folder=folder,result=result,quality=quality)

def td_normalized(a,ref):
    return np.diagonal(a['result']['T'],axis1=1,axis2=2)/np.diagonal(ref['result']['T'],axis1=1,axis2=2)

def run_p2(client,records):
    all_models={v['case_id']:v for v in d.read(d.OUT/'p1_models.json')}
    p0=d.read(d.OUT/'p0_audit.json');weak=f'group_seed_{p0["representative_order"][0]}'
    groups=[v for v in records if v['kind']=='group' and v['case_id']!=weak]
    worst=max(groups,key=lambda v:max(v['u_num']))['case_id']
    single=max((v for v in records if v['kind']=='single'),key=lambda v:max(np.array(v['score'])+v['u_num']))['case_id']
    models=[all_models[k] for k in ['reference',weak,worst,single]]
    level=max(v['level'] for v in records if v['case_id'] in [m['case_id'] for m in models])
    selection=dict(models=models,level=level,weakest_rule='Minimum min(signed S0 loss, signed A0 loss) in development',
                   different_adverse_rule='Largest measured P1 u_num among other representative groups',strongest_single_rule='Maximum score+u_num among P1 controls')
    d.freeze(d.OUT/'p2_models.json',selection)
    freq=d.read(d.OUT/'protocol.json')['transient']['frequencies_hz'];reference=models[0]
    d.write(d.OUT/'p2_workload.json',dict(geometries=4,frequencies=17,baseline_band_modal_solves=136,
        additional_mesh_absorber_modal_solves=136,transient_modal_solves=16,initial_window_us=280,maximum_window_us=560,
        note='Bounds before exact center-frequency/cache reuse; one entire-window extension is allowed.'))
    base={m['case_id']:w.fd_case(client,m,level,frequencies=freq,stage='P2_band') for m in models}
    checks=[]
    for model in [reference,all_models[worst]]:
        fine=w.fd_case(client,model,level+1,frequencies=freq,stage='P2_band_mesh')
        fref=w.fd_case(client,reference,level+1,frequencies=freq,stage='P2_band_mesh')
        absorb=w.fd_case(client,model,level,True,frequencies=freq,stage='P2_band_absorber')
        aref=w.fd_case(client,reference,level,True,frequencies=freq,stage='P2_band_absorber')
        a=w.normalized(base[model['case_id']],base['reference'])
        if model['kind']=='reference':
            a=np.diagonal(base['reference'][401]['T'],axis1=1,axis2=2)
            b=np.diagonal(fine[401]['T'],axis1=1,axis2=2);c=np.diagonal(absorb[401]['T'],axis1=1,axis2=2)
        else:b=w.normalized(fine,fref);c=w.normalized(absorb,aref)
        item=dict(case=model['case_id'],mesh=w.compare(a,b,model['kind']=='group'),absorber=w.compare(a,c,model['kind']=='group'))
        item['passed']=item['mesh']['passed'] and item['absorber']['passed'] and all(w.quality(x,model['kind']=='reference') for x in [base[model['case_id']],fine,absorb])
        checks.append(item);d.write(d.OUT/'p2_band_checks.json',checks)
        if not item['passed']:
            d.write(d.OUT/'p2_validation.json',dict(status='needs_adjustment',passed=False,band_checks=checks,reason='Full-band spatial control failed; transient/independent stages stopped.'))
            raise RuntimeError('P2 frequency-band spatial check failed')
    result=None
    for duration in [280e-6,560e-6]:
        td={(m['case_id'],half):transient_case(client,m,level,half,duration) for m in models for half in [False,True]}
        if any(not x['quality']['window_passed'] for x in td.values()):
            d.write(d.OUT/'p2_window_checks.json',dict(duration_s=duration,cases={f'{k[0]}_{k[1]}':v['quality'] for k,v in td.items()}))
            if duration<560e-6:continue
            d.write(d.OUT/'p2_validation.json',dict(status='needs_adjustment',passed=False,reason='Record window reached 560 us; no further extension.',windows={f'{k[0]}_{k[1]}':v['quality'] for k,v in td.items()}))
            raise RuntimeError('P2 reached its 560 us observation cap')
        result=[]
        for m in models:
            a=td_normalized(td[(m['case_id'],False)],td[('reference',False)])
            b=td_normalized(td[(m['case_id'],True)],td[('reference',True)])
            f=w.normalized(base[m['case_id']],base['reference'])
            if m['kind']=='reference':
                a=np.diagonal(td[('reference',False)]['result']['T'],axis1=1,axis2=2)
                b=np.diagonal(td[('reference',True)]['result']['T'],axis1=1,axis2=2)
                f=np.diagonal(base['reference'][401]['T'],axis1=1,axis2=2)
            checks_t=w.compare(a,b,m['kind']=='group');checks_f=w.compare(a,f,m['kind']=='group')
            q=[td[(m['case_id'],half)]['quality'] for half in [False,True]]
            channel_quality=all(x['residual']<=.01 and x['power_balance_error']<=.002 and (m['kind']!='reference' or x['right']<=.001) and
                all(y['relative_error'] is not None and y['relative_error']<=.05 for y in x['arrival']) for x in q)
            # The center endpoint classification must survive time/FD controls.
            classes=[(abs(1-abs(x[8]))>d.TAU).tolist() if m['kind']=='single' else (1-abs(x[8])>d.TAU).tolist() for x in [a,b,f]]
            consistent=classes[0]==classes[1]==classes[2]
            result.append(dict(case=m['case_id'],time_step=checks_t,frequency=checks_f,quality=q,center_raw_classification_consistent=consistent,
                passed=bool(checks_t['passed'] and checks_f['passed'] and channel_quality and consistent)))
        break
    passed=all(x['passed'] for x in result)
    d.write(d.OUT/'p2_validation.json',dict(status='complete' if passed else 'needs_adjustment',passed=passed,band_checks=checks,transient_checks=result,
        coverage='Selected representative geometries only, not all 30 independent specimens'))
    if not passed:raise RuntimeError('P2 transient or center-classification check failed')

def run_p3(client,records):
    protocol=d.read(d.OUT/'protocol.json');models=d.read(d.OUT/'frozen_geometry_manifest.json')
    reference=d.read(d.OUT/'p1_models.json')[0];initial=max(r['level'] for r in records)
    absorber_envelope=np.max([r['controls']['absorber']['magnitude_delta'] for r in records if r['kind']=='group'],axis=0)
    outputs=[]
    d.write(d.OUT/'p3_workload.json',dict(samples=30,initial_group_mesh_variants=60,initial_modal_solves=120,maximum_extra_refinement_levels=2))
    for entry in models:
        m=dict(case_id=f'independent_seed_{entry["seed"]}',kind='group',seed=entry['seed'],geometry_path=entry['path'],geometry_hash=entry['geometry_hash'])
        final=None
        for level in range(initial,initial+3):
            a=w.fd_case(client,m,level-1,stage='P3');b=w.fd_case(client,m,level,stage='P3')
            r0=w.fd_case(client,reference,level-1,stage='P3');r1=w.fd_case(client,reference,level,stage='P3')
            coarse=w.normalized(a,r0);fine=w.normalized(b,r1);dense=w.normalized(b,r1,801)
            mesh=w.compare(coarse,fine,True);ports=w.compare(fine,dense,True)
            u=np.array(mesh['magnitude_delta'])+ports['magnitude_delta']+absorber_envelope
            loss=1-abs(fine[0]);classes=[d.classify(float(v),float(e)) for v,e in zip(loss,u)]
            numeric=bool(mesh['passed'] and ports['passed'] and max(u)<=.0002 and w.quality(a) and w.quality(b) and w.quality(r0,True) and w.quality(r1,True))
            final=dict(seed=entry['seed'],case_id=m['case_id'],level=level,loss=loss.tolist(),u_num=u.tolist(),classifications=classes,
                numerical_passed=numeric,mesh=mesh,ports=ports,absorber_coverage='Envelope of measured representative P1 groups; not a certified bound for every specimen')
            d.write(d.OUT/'p3_checks'/m['case_id']/f'level_{level}.json',final)
            if numeric and 'indeterminate' not in classes:break
        if not final['numerical_passed']:
            d.write(d.OUT/'p3_validation.json',dict(status='needs_adjustment',passed=False,completed=outputs,failed=final,
                reason='Numerical acceptance failed after bounded refinement; fixed cohort remains incomplete.'))
            raise RuntimeError('P3 specimen numerical validation failed')
        outputs.append(final)
        d.write(d.OUT/'metrics'/f'p3_{m["case_id"]}.json',w.model_rows(m,b,r1,u,'Own mesh+sampling; representative absorber envelope (conditional scope)','independent_group'))
        d.write(d.OUT/'p3_validation.json',dict(status='running',passed=None,completed=outputs,fixed_denominator=30))
    stats=[]
    for i,mode in enumerate(d.MODES):
        k=sum(r['classifications'][i]=='above' for r in outputs);lower=d.cp_lower(k,30)
        stats.append(dict(mode=mode,successes=k,total=30,one_sided95_lower=lower,passed=lower>=.9,
            interpretation='Numerical threshold exceedance probability conditional on frozen geometry distribution and stated finite error estimates, not experimental POD'))
    passed=all(x['passed'] for x in stats)
    d.write(d.OUT/'p3_validation.json',dict(status='complete' if passed else 'complete_target_failed',passed=passed,completed=outputs,statistics=stats))

def run(client,records):
    run_p2(client,records)
    run_p3(client,records)
