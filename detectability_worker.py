"""Serial native frequency solves and P1 gates; protects the existing P queue."""
import os,time,json,hashlib,msvcrt,shutil,traceback,ctypes
from functools import lru_cache
import numpy as np
import detectability_200kHz_L48 as d
import guided_collective as gc
import run_transient as rt
from watch_refined_results import alive
from cfrp_pores.dispersion import Plate,validate_material
from cfrp_pores.pores import load_polygons
from cfrp_pores.portability import project_path,atomic_savez,validate_npz

def state(stage,**details):
    d.write(d.OUT/'worker_state.json',dict(updated_at=d.stamp(),stage=stage,**details))
    print(json.dumps(dict(stage=stage,**details),ensure_ascii=True),flush=True)

def protected_jobs():
    records=[]
    protected=d.ROOT/'output/transient_verified'
    stop=protected/'shutdown_stop_status.json'
    if (protected/'PAUSED_FOR_SHUTDOWN.json').exists() or (stop.exists() and d.read(stop).get('stage') in ['waiting_for_checkpoint','stopping_after_committed_checkpoint','ready_for_shutdown']):
        return [dict(pid=None,metadata=str(stop),reason='User shutdown/pause remains active; no automatic new solve')]
    for path in [d.ROOT/'output/transient_converged/pipeline_process.json',d.ROOT/'output/transient_verified/continuation_process.json']:
        if path.exists():
            v=d.read(path)
            if alive(v.get('launcher_pid',0),v.get('started')):records.append(dict(pid=v['launcher_pid'],created=v['started'],metadata=str(path)))
    if not records:
        restart=d.ROOT/'output/transient_converged/checkpoint_restart.json'
        if restart.exists() and d.read(restart).get('stage') in ['waiting_for_committed_block','restarting_after_commit']:
            return [dict(pid=None,metadata=str(restart),reason='Protected checkpoint restart in progress')]
        continuation=protected/'continuation_state.json'
        if continuation.exists() and d.read(continuation).get('stage') not in ['complete','needs_adjustment']:
            return [dict(pid=None,metadata=str(continuation),reason='Protected source has not reached a terminal state; inspect before releasing resources')]
    return records

def free_physical_gib():
    class MemoryStatus(ctypes.Structure):
        _fields_=[('dwLength',ctypes.c_ulong),('dwMemoryLoad',ctypes.c_ulong),
            ('ullTotalPhys',ctypes.c_ulonglong),('ullAvailPhys',ctypes.c_ulonglong),
            ('ullTotalPageFile',ctypes.c_ulonglong),('ullAvailPageFile',ctypes.c_ulonglong),
            ('ullTotalVirtual',ctypes.c_ulonglong),('ullAvailVirtual',ctypes.c_ulonglong),
            ('ullAvailExtendedVirtual',ctypes.c_ulonglong)]
    status=MemoryStatus();status.dwLength=ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        raise OSError('Cannot measure available physical memory')
    return status.ullAvailPhys/1024**3

def user_idle_minutes():
    class LastInput(ctypes.Structure):
        _fields_=[('cbSize',ctypes.c_uint32),('dwTime',ctypes.c_uint32)]
    value=LastInput();value.cbSize=ctypes.sizeof(value)
    if not ctypes.windll.user32.GetLastInputInfo(ctypes.byref(value)):
        raise OSError('Cannot measure time since last user input')
    return ((ctypes.windll.kernel32.GetTickCount()-value.dwTime)&0xffffffff)/60000

def apply_process_priority():
    policy_path=d.OUT/'resource_policy.json'
    policy=d.read(policy_path) if policy_path.exists() else {}
    if policy.get('below_normal_process_priority'):
        kernel=ctypes.WinDLL('kernel32',use_last_error=True)
        kernel.GetCurrentProcess.restype=ctypes.c_void_p
        kernel.SetPriorityClass.argtypes=[ctypes.c_void_p,ctypes.c_uint32]
        kernel.SetPriorityClass.restype=ctypes.c_int
        if not kernel.SetPriorityClass(kernel.GetCurrentProcess(),0x4000):
            raise OSError(ctypes.get_last_error(),'Cannot set the new guided worker below normal priority')

def resource_decision(jobs,initial=False,stage='P1'):
    policy_path=d.OUT/'resource_policy.json'
    policy=d.read(policy_path) if policy_path.exists() else {}
    if any(job.get('pid') is None for job in jobs):
        return False,'Protected P-wave pause or restart is pending',None
    if jobs and not policy.get('allow_concurrent_with_protected_p',False):
        return False,'Protected P-wave process still runs',None
    overlap=policy.get('overlap_stages')
    if jobs and overlap is not None and not any(stage.startswith(name) for name in overlap):
        return False,f'{stage} waits until the protected P-wave queue finishes',None
    minimum_idle=policy.get('require_user_idle_minutes')
    if minimum_idle is not None and user_idle_minutes()<float(minimum_idle):
        return False,f'Computer must be idle for {minimum_idle} minutes',None
    threshold_key='minimum_free_gib_before_session' if initial else 'minimum_free_gib_before_solve'
    threshold=policy.get('minimum_free_gib_before_solve_with_protected_p',policy.get(threshold_key)) if jobs and not initial else policy.get(threshold_key)
    if threshold is None:return True,None,None
    available=free_physical_gib()
    if available<float(threshold):return False,f'Available RAM below {threshold} GiB safeguard',available
    return True,None,available

def wait_resources(stage='P1',initial=False):
    last=None
    while True:
        jobs=protected_jobs()
        allowed,reason,available=resource_decision(jobs,initial,stage)
        if allowed:return
        key=json.dumps(dict(jobs=jobs,reason=reason),sort_keys=True)
        if key!=last:
            state('waiting_resources',pending_stage=stage,protected_jobs=jobs,available_ram_gib=available,error=reason)
            from detectability_report import render
            render();last=key
        time.sleep(30)

def verify_inputs():
    p=d.read(d.OUT/'protocol.json')
    for name,key in [('single_controls.json','single_controls_sha256'),('frozen_geometry_manifest.json','geometry_manifest_sha256'),('independent_seeds.json','seed_manifest_sha256')]:
        if d.digest(d.OUT/name)!=p[key]:raise RuntimeError(f'Frozen input changed: {name}')
    if d.digest(d.ROOT/'output/sample/geometry.json')!=p['source_geometry_sha256']:raise RuntimeError('Upstream contour library changed')
    for g in d.read(d.OUT/'frozen_geometry_manifest.json'):
        if d.digest(g['path'])!=g['file_sha256']:raise RuntimeError('Independent frozen geometry changed')
    manifest=d.OUT/'execution_manifest.json'
    if manifest.exists():
        for name,value in d.read(manifest)['engine_sha256'].items():
            if d.digest(d.ROOT/name)!=value:raise RuntimeError(f'Execution code changed after queue registration: {name}')

@lru_cache(maxsize=1)
def base_settings():
    protocol=d.read(d.OUT/'protocol.json')
    if d.digest(d.OUT/'protocol.json')!=d.read(d.OUT/'protocol.sha256.json')['sha256']:raise RuntimeError('Frozen protocol was modified')
    if d.read(d.ROOT/'config.json')!=protocol['material']:raise RuntimeError('Material changed after freezing')
    C,rho=validate_material(protocol['material']);s=dict(protocol['layout'],config=protocol['material'],C=C,rho=rho)
    s['plate']=Plate(.001,C[np.ix_([0,1,5],[0,1,5])],rho,120)
    s['plate'].at_frequency=lru_cache(maxsize=256)(s['plate'].at_frequency)
    return s

def settings(absorber=False):
    s=dict(base_settings())
    if absorber:s.update(absorb_length=s['absorb_length']*1.5,length=s['sponge']+s['absorb_length']*1.5)
    return s

def polys_for(model):
    if model['kind']=='reference':return []
    polys,_=load_polygons(project_path(model['geometry_path']))
    if d.geometry_hash(polys)!=model['geometry_hash']:raise RuntimeError('Frozen geometry changed')
    return polys

def old_folder(model,level,absorber):
    if absorber:return None
    if model['kind']=='reference':return d.OLD/('reference_L48' if level==0 else 'checks/reference_L48_fine') if level<=1 else None
    if model['kind']=='group' and model.get('seed') in range(101,113):
        if level==0:return d.dev_path(model['seed'])
        if level==1 and model['seed']==101:return d.OLD/'checks/d80_L48_fine'
    return None

def assemble(launches,positions):
    incoming=np.column_stack([x['pre_forward'] for x in launches]);cond=float(np.linalg.cond(incoming))
    if not np.isfinite(cond) or cond>100:raise RuntimeError('Ill-conditioned dual incident matrix')
    matrices={key:np.linalg.solve(incoming.T,np.column_stack([x[tag] for x in launches]).T).T for key,tag in [('T','post_forward'),('R','pre_backward'),('right','post_backward')]}
    p=np.asarray(launches[0]['powers']);energy=np.sum(p[:,None]/p[None,:]*(abs(matrices['T'])**2+abs(matrices['R'])**2),axis=0)
    return dict(**matrices,powers=p,energy=energy,condition=cond,residual=max(max(v['residuals']) for v in launches),right_max=float(abs(matrices['right']).max()))

def launch_from_raw(coeff,res,powers,wavelengths,positions):
    k=2*np.pi/np.asarray(wavelengths);v=dict(coefficients=coeff,residuals=res,powers=powers,wavelength_m=wavelengths)
    for port,sl in [('pre',slice(0,3)),('post',slice(3,6))]:
        v[port+'_forward']=np.mean(coeff[sl][:,[0,2]]*np.exp(1j*positions[sl,None]*k),axis=0)
        v[port+'_backward']=np.mean(coeff[sl][:,[1,3]]*np.exp(-1j*positions[sl,None]*k),axis=0)
    return v

def build(client,s,model,level):
    polys=polys_for(model);mesh=.00015/(2**level)
    minimum=min(mesh/30,min((np.sqrt(4*p.area/np.pi)*1e-6/(30*2**level) for p in polys),default=mesh/30))
    m=rt.create_2d(client,'detectability',s['config'],polys,size_m=(s['length'],s['h']),pore_offset_m=(s['start'],0),
        mesh_m=mesh,mesh_min_m=minimum,guided=True,sponge_start_m=s['sponge'],repair_tolerance_m=1e-9)
    geom=m.java.component('comp1').geom('geom1')
    count=int(geom.getNFiniteVoids())
    if count!=len(polys):client.remove(m);raise RuntimeError('Native CAD cavity count mismatch')
    measure=geom.measureFinal();measure.selection().geom('geom1',2);measure.selection().all()
    area=gc.matched_rectangle_area(client,s,1e-9)-float(measure.getArea())
    if len(polys)==1:area=gc.direct_pore_area(client,s,polys[0],1e-9)
    expected=sum(p.area for p in polys)*1e-12
    if abs(area-expected)/(.048*.001)>1e-4 or (expected and abs(area/expected-1)>.005):
        client.remove(m);raise RuntimeError('Native CAD area mismatch')
    return m,dict(pore_count=count,actual_area_m2=area,expected_area_m2=expected,
        absolute_area_fraction_error=abs(area-expected)/(.048*.001),boundary_selections=getattr(m,'_boundary_selection_checks',{}))

def frequency_source(model,s,frequency,mode):
    # The transient sponge is fixed at the pulse center. Keep exactly that
    # operator throughout the FD band while refreshing the SAFE launch shape.
    model.java.param().set('f0',f'{d.F0}[Hz]')
    model.java.study('std1').feature('freq').set('plist',str(frequency))
    gc.set_source(model,s['plate'].at_frequency(frequency,mode),s['h'])

def fd_case(client,model,level=1,absorber=False,frequencies=(200000.,),stage='P1'):
    frequencies=np.asarray(frequencies,float);band=len(frequencies)>1
    folder=d.OUT/'solves'/('band' if band else 'center')/model['case_id']/f'level_{level}'/('absorber' if absorber else 'base')
    folder.mkdir(parents=True,exist_ok=True);s=settings(absorber);positions=np.array(s['pre_positions']+s['post_positions'])
    sig=dict(geometry_hash=model['geometry_hash'],material=s['config'],layout=[s['start'],s['sponge'],s['length']],
        port_positions_m=positions.tolist(),solver_schema='uniform_ports_2source_matrix_v1',mesh=.00015/2**level,
        micro_mesh_divisor=30*2**level,absolute_repair_tolerance_m=1e-9,frequencies=frequencies.tolist(),primary_samples=401,additional_samples_for_native_fields=[801,1601],
        source='Frequency-specific SAFE source, conjugated exp(+iwt)',source_code_sha256=d.digest(d.ROOT/'guided_collective.py'))
    d.freeze(folder/'signature.json',sig)
    if (folder/'complete.json').exists():return load_fd(folder)
    old=old_folder(model,level,absorber);oldrecord=d.read(old/'results.json') if old else None
    if old:
        osig=oldrecord['signature']
        if osig['geometry_hash']!=sig['geometry_hash'] or osig['material']!=sig['material'] or osig['layout']!=sig['layout']:raise RuntimeError('Legacy result is not compatible')
        with np.load(old/'signals.npz') as z:np.testing.assert_allclose(z['positions_m'],positions,rtol=0,atol=1e-12)
    native=None;geometry=None
    timings=d.read(folder/'timings.json') if (folder/'timings.json').exists() else []
    sources=[d.read(p) for p in sorted(folder.glob('*.source.json'))]
    try:
        for f in frequencies:
            for mode in d.MODES:
                file=folder/f'{int(f)}_{mode}.npz'
                if file.exists():
                    validate_npz(file,required=('coeff_401','res_401','powers','wavelengths'))
                    continue
                center=d.OUT/'solves/center'/model['case_id']/f'level_{level}'/('absorber' if absorber else 'base')
                if band and f==d.F0 and (center/'complete.json').exists():
                    expected=dict(sig,frequencies=[d.F0])
                    if d.read(center/'signature.json')==expected:
                        cached=load_fd(center)
                        shutil.copy2(center/file.name,file)
                        d.write(folder/f'{int(f)}_{mode}.source.json',dict(frequency=f,mode=mode,reuse='Identical completed P1 center-frequency launch',source=str(center/file.name),sha256=d.digest(file)))
                        continue
                # Baseline cached raw coefficients are enough for a mesh pair;
                # fine solves carry denser sampling, including both incident modes.
                old_has=oldrecord and any(x['frequency_hz']==f for x in oldrecord['results'])
                if old_has and level==0:
                    row=next(x for x in oldrecord['results'] if x['frequency_hz']==f)
                    with np.load(old/'signals.npz') as z:
                        atomic_savez(file,coeff_401=z[f'{int(f)}_{mode}_coefficients'],res_401=z[f'{int(f)}_{mode}_residuals'],
                            powers=row['modal_power_flux'],wavelengths=row['wavelength_m'])
                    provenance=dict(frequency=f,mode=mode,reuse='raw dual-launch data',source=str(old),fingerprint=oldrecord['fingerprint'])
                    sources.append(provenance);d.write(folder/f'{int(f)}_{mode}.source.json',provenance)
                    continue
                wait_resources(f'{stage}_level_{level}')
                state('solving',stage_name=stage,case=model['case_id'],level=level,absorber=absorber,frequency_hz=f,mode=mode)
                existing_model=None
                if old_has and mode=='A0' and f==oldrecord['signature']['frequencies'][-1] and (old/'model.mph').exists():
                    existing_model=client.load(str((old/'model.mph').resolve()))
                    np.testing.assert_allclose(existing_model.java.param().evaluate('f0'),f,rtol=1e-12)
                    data=gc.modal_port_data(existing_model,s,f)
                    with np.load(old/'signals.npz') as z:
                        a=z[f'{int(f)}_A0_coefficients'];b=data['coefficients']
                        if np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-100)>1e-7:raise RuntimeError('Saved A0 field does not reproduce original probes')
                    current=existing_model;elapsed=0.
                    if int(existing_model.java.component('comp1').geom('geom1').getNFiniteVoids())!=len(polys_for(model)):raise RuntimeError('Saved native geometry cavity count mismatch')
                    provenance=dict(frequency=f,mode=mode,reuse='verified last-frequency A0 native field',source=str(old/'model.mph'),sha256=d.digest(old/'model.mph'))
                    sources.append(provenance);d.write(folder/f'{int(f)}_{mode}.source.json',provenance)
                else:
                    if native is None:
                        native,geometry=build(client,s,model,level)
                        native.java.component('comp1').physics('solid').feature('lemm1').feature('sponge').set('alpha_dM',
                            f'if(x>{s["sponge"]}[m],3*2*pi*f0*((x-{s["sponge"]}[m])/{s["absorb_length"]}[m])^3,0[1/s])')
                        gc.add_frequency_study(native,f)
                    frequency_source(native,s,f,mode)
                    started=time.perf_counter();native.solve();elapsed=time.perf_counter()-started;current=native
                    data=gc.modal_port_data(current,s,f)
                    d.write(folder/f'{int(f)}_{mode}.source.json',dict(frequency=f,mode=mode,reuse='new native solve',signature_sha256=d.signature(sig)))
                    timings.append(dict(case=model['case_id'],level=level,absorber=absorber,frequency_hz=f,mode=mode,seconds=elapsed))
                    d.write(folder/'timings.json',timings)
                    if not (d.OUT/'first_native_timing.json').exists():d.write(d.OUT/'first_native_timing.json',timings[-1])
                try:
                    arrays=dict(powers=data['powers'],wavelengths=data['wavelength_m'],coeff_401=data['coefficients'],res_401=data['residuals'])
                    for count in [801,1601]:
                        other=gc.modal_port_data(current,s,f,samples=count);arrays[f'coeff_{count}']=other['coefficients'];arrays[f'res_{count}']=other['residuals']
                    atomic_savez(file,**arrays)
                    d.write(folder/'progress.json',dict(last_frequency=f,last_mode=mode,last_solve_seconds=elapsed))
                finally:
                    if existing_model is not None:client.remove(existing_model)
    finally:
        if native is not None:client.remove(native)
    d.write(folder/'complete.json',dict(signature_sha256=d.signature(sig),geometry=geometry,legacy_sources=sources,finished_at=d.stamp(),
        absorber_anchor_frequency_hz=d.F0,
        launch_sha256={p.name:d.digest(p) for p in folder.glob('*.npz')}))
    return load_fd(folder)

def load_fd(folder):
    sig=d.read(folder/'signature.json');positions=np.array(sig['port_positions_m']);frequencies=sig['frequencies'];result=dict(folder=folder,signature=sig,frequencies=np.array(frequencies))
    complete=d.read(folder/'complete.json')
    if complete['signature_sha256']!=d.signature(sig):raise RuntimeError('Completed signature mismatch')
    for name,expected in complete.get('launch_sha256',{}).items():
        if d.digest(folder/name)!=expected:raise RuntimeError('Saved launch data changed')
    for count in [401,801,1601]:
        values=[]
        for f in frequencies:
            launches=[]
            for mode in d.MODES:
                validate_npz(folder/f'{int(f)}_{mode}.npz',required=('coeff_401','res_401','powers','wavelengths'))
                with np.load(folder/f'{int(f)}_{mode}.npz',allow_pickle=False) as z:
                    if f'coeff_{count}' not in z:break
                    launches.append(launch_from_raw(z[f'coeff_{count}'],z[f'res_{count}'],z['powers'],z['wavelengths'],positions))
            if len(launches)!=2:break
            values.append(assemble(launches,positions))
        if len(values)==len(frequencies):result[count]={k:np.array([v[k] for v in values]) for k in values[0]}
    return result

def normalized(a,r,count=401):return np.diagonal(a[count]['T'],axis1=1,axis2=2)/np.diagonal(r[count]['T'],axis1=1,axis2=2)
def compare(a,b,group):
    mag=np.max(abs(abs(a)-abs(b)),axis=0);phase=np.max(abs(np.angle(a/b)),axis=0)*180/np.pi
    relative=np.max(abs(abs(a)-abs(b))/np.maximum(np.maximum(abs(1-abs(a)),abs(1-abs(b))),1e-12),axis=0)
    return dict(magnitude_delta=mag.tolist(),phase_deg=phase.tolist(),relative_loss_change=relative.tolist(),
        passed=bool(max(mag)<=.0002 and max(phase)<=1 and (not group or max(relative)<=.05)))
def quality(a,reference=False):
    x=a[401];return bool(np.max(x['residual'])<=.01 and np.max(abs(x['energy']-1))<=.002 and (not reference or np.max(x['right_max'])<=.001))

def model_rows(model,fine,ref,u,coverage,category):
    idx=int(np.argmin(abs(fine['frequencies']-d.F0)))
    def convert(x):
        data=x[401]
        return dict(folder=x['folder'],signature=x['signature'],record=dict(fingerprint=d.signature(x['signature'])),positions=np.array(x['signature']['port_positions_m']),
            **{k:data[k][idx] for k in ['T','R','right','powers','energy','condition','residual','right_max']})
    return d.metrics(convert(fine),convert(ref),model['case_id'],category,model.get('seed'),u,coverage)

def run_p1(client):
    models=d.read(d.OUT/'p1_models.json');reference=models[0];records=[]
    d.write(d.OUT/'p1_workload.json',dict(geometries=len(models),initial_variants=3*len(models),initial_modal_solves_upper_bound=6*len(models),
        refinements_max=2,first_measured_timing=None,note='Upper bound before reuse; each solve preserves all three port sampling densities.'))
    for model in models[1:]:
        accepted=False
        for level in [1,2,3]:
            r0=fd_case(client,reference,level-1);r1=fd_case(client,reference,level);ra=fd_case(client,reference,level,True)
            coarse=fd_case(client,model,level-1);fine=fd_case(client,model,level);extended=fd_case(client,model,level,True)
            a,b,c=normalized(coarse,r0),normalized(fine,r1),normalized(extended,ra)
            dense=normalized(fine,r1,801)
            controls=dict(mesh=compare(a,b,model['kind']=='group'),absorber=compare(b,c,model['kind']=='group'),ports=compare(b,dense,model['kind']=='group'))
            u=np.sum([v['magnitude_delta'] for v in controls.values()],axis=0);loss=1-abs(b[0]);score=abs(loss)
            numerical=all(v['passed'] for v in controls.values()) and max(u)<=.0002 and all(quality(x,True) for x in [r0,r1,ra]) and all(quality(x) for x in [coarse,fine,extended])
            margin=bool(np.all(u<abs(loss-d.TAU)))
            single=bool(np.all(score+u<d.TAU)) if model['kind']=='single' else True
            check=dict(case_id=model['case_id'],kind=model['kind'],level=level,controls=controls,u_num=u.tolist(),loss=loss.tolist(),score=score.tolist(),
                numerical_passed=bool(numerical),threshold_margin_resolved=margin,single_control_passed=single if model['kind']=='single' else None,
                passed=bool(numerical and margin and single),contour_scope='fixed upstream polygon geometry',source=str(fine['folder']))
            d.write(d.OUT/'p1_checks'/model['case_id']/f'level_{level}.json',check)
            d.write(d.OUT/'metrics'/f'p1_{model["case_id"]}.json',model_rows(model,fine,r1,u,'Own paired mesh + absorber + 401/801 sampling; finite sensitivity estimate','matched_single' if model['kind']=='single' else 'representative_group'))
            if check['passed']:accepted=True;records.append(check);break
        if not accepted:
            gate=dict(status='needs_adjustment',passed=False,completed=records,failed=check,reason='P1 failed after initial controls and at most two additional levels; dependent stages stopped.')
            d.write(d.OUT/'p1_validation.json',gate);raise RuntimeError(gate['reason'])
        d.write(d.OUT/'p1_validation.json',dict(status='running',passed=None,completed=records))
    d.write(d.OUT/'p1_validation.json',dict(status='complete',passed=True,completed=records,all_single_controls_passed=True))
    return records

def main():
    d.OUT.mkdir(exist_ok=True);lock=(d.OUT/'.worker.lock').open('a+b');lock.seek(0)
    try:msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
    except OSError:raise RuntimeError('A detectability worker already owns this task')
    client=None
    try:
        verify_inputs()
        apply_process_priority()
        wait_resources(stage='P1_level_1',initial=True)
        verify_inputs()
        state('starting_P1',cores=2)
        client=rt.start_session(d.OUT,2)
        records=run_p1(client)
        verify_inputs()
        from detectability_followon import run
        run(client,records)
        state('complete')
    except Exception as e:
        current=d.read(d.OUT/'worker_state.json') if (d.OUT/'worker_state.json').exists() else {}
        phase=str(current.get('stage_name',current.get('pending_stage','P1')))[:2].lower()
        if phase in ['p1','p2','p3']:
            file=d.OUT/f'{phase}_validation.json';gate=d.read(file) if file.exists() else {}
            if gate.get('passed') is not False:
                gate.update(status='needs_adjustment',passed=False,failure_type='execution_error',reason=str(e));d.write(file,gate)
        d.write(d.OUT/'worker_error.json',dict(error=str(e),traceback=traceback.format_exc(),at=d.stamp()))
        state('needs_adjustment',error=str(e))
        raise
    finally:
        from detectability_report import render
        render();lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()

if __name__=='__main__':main()
