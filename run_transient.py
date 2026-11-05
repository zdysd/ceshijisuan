"""Matched-pore Solid Mechanics transient runs. Use the ASCII COMSOL Python.

Stages: pilot, spatial, verify, ensemble, report, audit. Validation gates control ensemble runs.
Existing frequency-domain models/results are never overwritten.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import mph
from cfrp_pores.portability import runtime_root,atomic_savez
import site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import argparse,json,hashlib,time,sys
from functools import lru_cache
import numpy as np
from cfrp_pores.comsol import create_2d,add_load,add_modal_displacement,add_frequency_study,interpolate,compact_transient,retain_transient_frames
from cfrp_pores.dispersion import Plate,validate_material
from cfrp_pores.transient import source_band,analyze_ports,modal_bank,tail_fraction,field_frame_indices,minimum_inplane_speed
from guided_collective import grouped_pores
from guided_statistics import attenuation_interval

SCHEMA='transient_matched_ports_v1'
MODES={'bulk':['P'],'guided':['S0','A0']}


def write_json(path,value):
    path=Path(path);temporary=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    temporary.replace(path)


def configure_scratch(root):
    # COMSOL's native JVM option handling on this Windows host misdecodes
    # non-ASCII directory names. Keep runtime paths ASCII; results remain local.
    scratch=runtime_root()
    temporary=scratch/f'temp_{os.getpid()}';recovery=scratch/f'recovery_{os.getpid()}'
    preferences=scratch/f'preferences_{os.getpid()}'
    for path in [temporary,recovery,preferences]:path.mkdir(parents=True,exist_ok=True)
    os.environ['TEMP']=str(temporary);os.environ['TMP']=str(temporary)
    import tempfile
    tempfile.tempdir=str(temporary)
    os.environ['MPLCONFIGDIR']=str(scratch/'matplotlib')
    return scratch,temporary,recovery,preferences


def start_session(root,cores):
    """Use process-local E-drive scratch and isolated COMSOL preferences."""
    runtime=root/'_runtime';runtime.mkdir(parents=True,exist_ok=True)
    scratch,temporary,recovery,preferences=configure_scratch(root)
    # The recovery manager can cache its directory during JVM initialization.
    # Supply an isolated preferences file and JVM options BEFORE initStandalone.
    escaped_recovery=recovery.as_posix().encode('unicode_escape').decode('ascii')
    (preferences/'comsol.prefs').write_text('tempfiles.recovery.recoverydir='+escaped_recovery+'\n'
        'tempfiles.recovery.autosave=on\ntempfiles.recovery.checkforrecoveries=off\n',encoding='ascii')
    import jpype
    original_start=jpype.startJVM
    def start_with_paths(*positional,**keywords):
        return original_start(*positional,f'-Djava.io.tmpdir={temporary}',f'-Dcs.tmpdir={temporary}',
                              f'-Dcs.prefsdir={preferences}',f'-Dcs.recoverydir={recovery}',**keywords)
    jpype.startJVM=start_with_paths
    try:client=mph.Client(cores=cores,version='6.2')
    finally:jpype.startJVM=original_start
    client.java.setPreference('tempfiles.recovery.recoverydir',str(recovery))
    client.java.setPreference('tempfiles.recovery.autosave','on')
    actual_temp=str(jpype.JClass('java.lang.System').getProperty('java.io.tmpdir'))
    actual_recovery=str(client.java.getPreference('tempfiles.recovery.recoverydir'))
    if not all(Path(v).resolve().is_relative_to(scratch.resolve()) for v in [actual_temp,actual_recovery]):
        raise RuntimeError('COMSOL scratch paths are not on the requested output drive.')
    probe=client.create('scratch_path_probe')
    actual_model_recovery=str(client.java.getRecoveryDir(probe.java.tag()))
    # Check the model's actual directory, not only the preference readback.
    checked_model_path=Path(actual_model_recovery).resolve()
    checked_recovery_path=recovery.resolve()
    valid_recovery=checked_model_path.is_relative_to(checked_recovery_path)
    client.remove(probe)
    if not valid_recovery:
        raise RuntimeError(f'Actual model recovery path {checked_model_path!s} is outside {checked_recovery_path!s}')
    client.java.showProgress(str(scratch/f'solver_{os.getpid()}.log'))
    write_json(runtime/f'session_{os.getpid()}.json',dict(pid=os.getpid(),temporary_directory=actual_temp,
               recovery_directory=actual_recovery,actual_model_recovery=actual_model_recovery,
               isolated_preferences_directory=str(preferences),preferences_saved_to_user_profile=False,
               automatic_recovery_enabled=bool(client.java.useRecovery())))
    return client


def layout(kind,length_mm,f0,C,rho,variant='base',profile='legacy'):
    low,high=source_band(f0);L=length_mm*1e-3
    if kind=='bulk':
        start=.002;gap=.0005;absorb=.002
        speed=np.sqrt(C[0,0]/rho)
        shear_speed=np.sqrt(C[5,5]/rho)
        end_time=2*(start+L+3*gap)/min(speed,shear_speed)+5/f0
        # Initial quadratic mesh: at least five elements per shortest oblique
        # wavelength. The mandatory half-size check determines actual accuracy.
        mesh=min(minimum_inplane_speed(C[np.ix_([0,1,5],[0,1,5])],rho)/high/5,4e-5)
    else:
        start=.020;gap=.002;absorb=.040
        end_time=140e-6+max(0,L-.008)/1000*3
        mesh=.00015
    if profile in ['refined','refined_v2','refined_v3','refined_v4','refined_v5','refined_v6']:
        mesh=.00015 if kind=='guided' else .00006
        absorb=.080 if kind=='guided' else .004
        if kind=='bulk':start=.007
        if profile in ['refined_v2','refined_v3','refined_v4','refined_v5','refined_v6'] and kind=='bulk':
            mesh=.00005;absorb=.008;start=.011
            if variant=='absorber':start+=.004
            if profile in ['refined_v3','refined_v4','refined_v5','refined_v6']:mesh=.000075
            if profile in ['refined_v5','refined_v6']:
                absorb=.012 if profile=='refined_v5' else .018
                start=absorb*(1.5 if variant=='absorber' else 1)+.003
        end_time=(140e-6+max(0,L-.008)/1000*3 if kind=='guided' else 50e-6+max(0,L-.008)/1000*5)
        if profile=='refined_v6' and kind=='bulk':end_time*=2
    post=start+L+np.array([1,2,3])*gap
    sponge=float(post[-1]+2*gap)
    if variant=='absorber':absorb*=1.5
    value=dict(kind=kind,f0=f0,fmin=low,fmax=high,roi_length=L,start=start,
                width=.001,sponge=sponge,length=sponge+absorb,absorb_length=absorb,
                positions=np.r_[start-np.array([3,2,1])*gap,post].tolist(),
                end_time=end_time,dt=(.5 if variant=='time' else 1)/(60*high),
                output_dt=1/(20*high),mesh=mesh*(.5 if variant=='mesh' else 1),variant=variant)
    if profile in ['refined','refined_v2','refined_v3','refined_v4','refined_v5','refined_v6']:
        value.update(profile=profile if kind=='bulk' else 'refined',displacement_order=3,curvature_factor=.15,port_samples=129,safe_elements=120,
                     dt=(.5 if variant=='time' else 1)/((60 if kind=='guided' else 200 if profile in ['refined_v4','refined_v5','refined_v6'] else 160)*high),damping_factor=.5 if kind=='bulk' else 1.)
        if kind=='guided':value['block_duration']=30e-6
        if kind=='bulk':value.update(source_x=.0045,left_sponge_end=.004,block_duration=5e-6)
        if profile in ['refined_v2','refined_v3','refined_v4','refined_v5','refined_v6'] and kind=='bulk':
            left=.012 if variant=='absorber' else .008
            if profile in ['refined_v5','refined_v6']:left=absorb
            value.update(source_x=left+.0005,left_sponge_end=left,adaptive_window=True,minimum_window=20e-6)
            if profile in ['refined_v3','refined_v4','refined_v5','refined_v6']:value.update(displacement_order=4,boundary_selection_revision=2)
            if profile in ['refined_v4','refined_v5','refined_v6']:value['periodic_mesh_match']=True
    if variant=='sampling':value['output_dt']=value['dt']
    return value


def builder_options(s):
    return dict(source_x_m=s.get('source_x'),displacement_order=s.get('displacement_order'),
                curvature_factor=s.get('curvature_factor',.35),periodic_mesh_match=s.get('periodic_mesh_match',False))


def analysis_frequencies(s):
    return np.linspace(s['fmin'],s['fmax'],17) if s.get('profile') else s['f0']*np.linspace(.7,1.3,17)


def cache_linear_matrices(model):
    """K, damping, mass and constraint Jacobian are time independent here.

    Loads and prescribed displacement values still depend on time and MUST
    be reassembled. This changes assembly cost, not the mechanical equations.
    """
    step=model.java.sol('sol1').feature('t1')
    tags=[str(t) for t in step.feature().tags() if str(step.feature(t).getType())=='Advanced']
    advanced=step.feature(tags[0]) if tags else step.create('constant_matrices','Advanced')
    advanced.set('keep','on')
    for key in ['K','D','E','N']:advanced.set(key,'on')
    advanced.set('L','off')
    advanced.set('M','off')
    advanced.set('cachepattern',True)


def source_and_study(model,s,plate,mode):
    j=model.java;p=j.component('comp1').physics('solid')
    j.param().set('f0',f'{s["f0"]}[Hz]')
    j.param().set('dt',f'{s["dt"]}[s]');j.param().set('tout',f'{s["output_dt"]}[s]')
    j.param().set('tend',f'{s["end_time"]}[s]')
    envelope='if(t<=5/f0,(1-cos(2*pi*f0*t/5))/2,0)'
    if mode=='P':add_load(model,[f'1e6[Pa]*({envelope})*sin(2*pi*f0*t)','0','0'])
    else:
        wave=plate.at_frequency(s['f0'],mode)
        wave=dict(wave,u=np.conj(wave['u']),v=np.conj(wave['v']))
        add_modal_displacement(model,wave,s['width'])
        # Re[(a+ib) exp(i wt)] = a*cos(wt)-b*sin(wt), preserving quadrature.
        p.feature('source').set('U0',[f'1e-8[m]*({envelope})*({a}r(y)*cos(2*pi*f0*t)-{a}i(y)*sin(2*pi*f0*t))' for a in ['u','v']]+['0'])
    if 'source_x' in s:p.feature('source').selection().named('excitation')
    if s['kind']=='guided' or 'left_sponge_end' in s:
        factor=s.get('damping_factor',3.)
        right=f'if(x>{s["sponge"]}[m],{factor}*2*pi*f0*((x-{s["sponge"]}[m])/{s["absorb_length"]}[m])^3,0[1/s])'
        left=(f'+if(x<{s["left_sponge_end"]}[m],{factor}*2*pi*f0*(({s["left_sponge_end"]}[m]-x)/{s["left_sponge_end"]}[m])^3,0[1/s])'
              if 'left_sponge_end' in s else '')
        p.feature('lemm1').feature('sponge').set('alpha_dM',
            right+left)
    j.study().create('std1');j.study('std1').create('time','Transient')
    j.study('std1').feature('time').set('tlist','range(0,tout,tend)')
    j.study('std1').createAutoSequences('all')
    step=j.sol('sol1').feature('t1');step.set('timemethod','genalpha')
    step.set('tstepsgenalpha','manual');step.set('timestepgenalpha','dt');step.set('rhoinf',1.)
    step.set('rtol',1e-5)
    if 'd1' not in list(step.feature().tags()):step.create('d1','Direct')
    step.feature('d1').set('linsolver','pardiso')
    if 'fc1' in list(step.feature().tags()):step.feature('fc1').set('linsolver','d1')
    if s.get('profile'):cache_linear_matrices(model)


def extract_history(model,s):
    y=np.linspace(1e-8,s['width']-1e-8,s.get('port_samples',65));positions=np.array(s['positions'])
    xgrid,ygrid=np.meshgrid(positions,y,indexing='ij')
    fields=interpolate(model,[xgrid.ravel(),ygrid.ravel()],['u','v','solid.sx','solid.sxy']).real
    t=np.asarray(model.evaluate('t')).ravel()
    if not len(t) or t[-1]<s['end_time']-s['output_dt']*(1+1e-8):
        raise RuntimeError('The solver did not reach the requested end of the observation window.')
    values=fields.transpose(1,2,0).reshape(len(t),len(positions),len(y),4).transpose(0,1,3,2)
    if not np.isfinite(values).all():raise RuntimeError('Receiver intersects a void or a nonfinite field.')
    return dict(times=t,values=values,y=y,positions=positions)


def extract_movie(model,s,folder,mode,times):
    # All frames are real solved transient fields; no phase animation of a harmonic solution.
    x=np.linspace(0,s['sponge'],360);y=np.linspace(1e-8,s['width']-1e-8,17)
    xx,yy=np.meshgrid(x,y,indexing='xy');selected=field_frame_indices(times,s['fmax'])
    # Select solution numbers before interpolation to avoid materializing every full field.
    import jpype
    j=model.java;n=j.result().numerical().create('movie_sample','Interp')
    n.set('data','dset1');n.set('expr',['u','v'])
    n.set('coord',jpype.JArray(jpype.JDouble,2)([xx.ravel().tolist(),yy.ravel().tolist()]))
    n.set('solnum',jpype.JArray(jpype.JInt)((selected+1).tolist()))
    data=np.asarray(n.getData(),float).reshape(2,len(selected),len(y),len(x))
    destination=folder/f'{mode}_fields.npz';reuse=False
    if destination.exists():
        with np.load(destination) as previous:
            reuse=(np.array_equal(previous['times'],times[selected]) and previous['u'].shape==data[0].shape
                   and np.allclose(previous['u'],data[0],rtol=1e-10,atol=1e-18,equal_nan=True)
                   and np.allclose(previous['v'],data[1],rtol=1e-10,atol=1e-18,equal_nan=True))
    if not reuse:atomic_savez(destination,times=times[selected],x=x,y=y,u=data[0],v=data[1])
    else:print('REUSE numerically matching movie grid',mode,flush=True)
    j.result().numerical().remove('movie_sample')


def frequency_checkpoint(model,s,folder,mode,frequencies=None):
    """Same geometry, mesh, source mode and absorbing layer; independent harmonic solve."""
    j=model.java;p=j.component('comp1').physics('solid')
    if mode=='P':p.feature('source').set('FperArea',['1e6[Pa]','0','0'])
    else:p.feature('source').set('U0',['1e-8[m]*(ur(y)+i*ui(y))','1e-8[m]*(vr(y)+i*vi(y))','0'])
    j.study().remove('std1');j.sol().remove('sol1')
    for tag in list(j.result().dataset().tags()):j.result().dataset().remove(tag)
    frequencies=np.asarray([s['f0']] if frequencies is None else frequencies,float)
    add_frequency_study(model,s['f0'])
    j.study('std1').feature('freq').set('plist',' '.join(f'{f:.17g}' for f in frequencies))
    model.solve()
    y=np.linspace(1e-8,s['width']-1e-8,s.get('port_samples',65));positions=np.array(s['positions'])
    xx,yy=np.meshgrid(positions,y,indexing='ij')
    fields=interpolate(model,[xx.ravel(),yy.ravel()],['u','v','solid.sx','solid.sxy'])
    values=fields.transpose(1,2,0).reshape(len(frequencies),len(positions),len(y),4).transpose(0,1,3,2)
    center=int(np.argmin(abs(frequencies-s['f0'])))
    if not np.isclose(frequencies[center],s['f0']):raise ValueError('Harmonic checkpoint must include the center frequency.')
    atomic_savez(folder/f'{mode}_harmonic.npz',values=values[center],y=y,positions=positions)
    if len(frequencies)>1:
        atomic_savez(folder/f'{mode}_harmonic_band.npz',values=values,y=y,positions=positions,frequencies=frequencies)


def ensure_recovery_stopped(source):
    import ctypes
    from ctypes import wintypes
    owner=(source/'recoveryid.txt').read_text().split(':')
    pid=int(owner[1])
    kernel=ctypes.windll.kernel32;kernel.OpenProcess.restype=ctypes.c_void_p
    kernel.CloseHandle.argtypes=[ctypes.c_void_p]
    kernel.GetExitCodeProcess.argtypes=[ctypes.c_void_p,ctypes.POINTER(wintypes.DWORD)]
    kernel.GetProcessTimes.argtypes=[ctypes.c_void_p]+[ctypes.POINTER(wintypes.FILETIME)]*4
    handle=kernel.OpenProcess(0x1000,0,pid)
    if handle:
        try:
            status=wintypes.DWORD()
            if not kernel.GetExitCodeProcess(handle,ctypes.byref(status)):raise RuntimeError('Cannot verify recovery owner state.')
            if status.value!=259:return  # An exited process can still have an open Windows handle.
            created,exited,kernel_time,user_time=[wintypes.FILETIME() for _ in range(4)]
            if not kernel.GetProcessTimes(handle,*[ctypes.byref(v) for v in [created,exited,kernel_time,user_time]]):
                raise RuntimeError('Cannot verify recovery owner creation time.')
            creation_seconds=((created.dwHighDateTime<<32)|created.dwLowDateTime)/1e7-11644473600
            if creation_seconds>(source/'recoveryid.txt').stat().st_mtime+1:return  # PID was reused by a newer process.
            raise RuntimeError(f'Case is still owned by running process {pid}; refusing duplicate solve.')
        finally:kernel.CloseHandle(handle)


def recover_mode(client,record,mode):
    """Copy a stopped process's real recovery model before continuing it."""
    if record.get('active_mode')!=mode:return None
    source=Path(record.get('active_recovery_directory',''))
    if not (source/'recovery.xml').is_file():return None
    scratch=runtime_root().resolve()
    if not source.resolve().is_relative_to(scratch):raise ValueError('Recovery path is outside this simulation scratch directory.')
    import shutil
    ensure_recovery_stopped(source)
    destination=scratch/f'restored_input_{os.getpid()}_{time.time_ns()}'
    print('COPY RECOVERY',source,flush=True)
    shutil.copytree(source,destination)
    model=mph.Model(client.java.loadRecovery('restored_'+mode,str(destination)))
    s=record['signature']['layout']
    if not np.isclose(model.java.param().evaluate('f0'),s['f0']):raise RuntimeError('Recovery frequency does not match the case.')
    return model


def solve_in_blocks(model,s,folder,mode,frames):
    """Export each short absolute-time block, then retain its final u and velocity."""
    import jpype
    blocks=folder/'blocks'/mode;blocks.mkdir(parents=True,exist_ok=True)
    progress=blocks/'progress.json'
    state=json.loads(progress.read_text(encoding='utf-8-sig')) if progress.exists() else dict(completed=0,checks=[])
    interval=s['output_dt'];total=int(np.floor(s['end_time']/interval+1e-8))
    step=max(1,int(s['block_duration']/interval));ends=list(range(step,total,step))+[total]
    off_index=int(np.ceil((5/s['f0'])/interval))+2
    if s.get('adaptive_window') and off_index<ends[0]:ends=[off_index]+ends
    if s['kind']=='guided':
        # Restart the prescribed-displacement launch only after it is off.
        # During excitation its constrained velocity requires extra consistent
        # initialization; we do not assume that getUDot alone provides it.
        ends=[n for n in ends if n==total or n*interval>5/s['f0']+2*interval]
    if state['completed']:
        from cfrp_pores.transient import continuation_block_ends
        saved=[]
        for i in range(state['completed']):
            with np.load(blocks/f'history_{i:03d}.npz') as h:
                last=float(h['times'][-1]);position=int(round(last/interval))
                if not np.isclose(last,position*interval,rtol=1e-12,atol=1e-17):
                    raise RuntimeError('Saved block is not on the requested output grid.')
                saved.append(position)
        ends=continuation_block_ends(ends,saved)
    j=model.java;sol=j.sol('sol1')
    C,rho=validate_material(json.loads((ROOT/'config.json').read_text()))
    stress_scale=2*np.pi*s['f0']*np.sqrt(rho*C[0,0]);weights=np.array([1.,1.,stress_scale,stress_scale])[None,:,None]
    for index in range(state['completed'],len(ends)):
        if state.get('window_converged'):break
        first=0 if index==0 else ends[index-1];last=ends[index]
        t0=first*interval;t1=last*interval;expected_velocity=None;prior=None
        if index:
            prior=dict(np.load(blocks/f'history_{index-1:03d}.npz'))
            np.testing.assert_allclose(np.asarray(model.evaluate('t')).ravel()[-1],t0,rtol=1e-12,atol=1e-17)
            expected_velocity=np.asarray(sol.getUDot(jpype.JInt(1)))
            init=f'block_init_{index}';sol.copySolution(init)
            v=j.sol('sol1').feature('v1');v.set('control','user');v.set('initmethod','sol')
            v.set('initsol',init);v.set('solnum','last')
            # The imported u and velocity already satisfy the unchanged
            # constraints. A backward-Euler consistency correction would
            # perturb the restart state, particularly for a prescribed source.
            sol.feature('t1').set('consistent','off')
            if index>1:j.sol().remove(f'block_init_{index-1}')
        expression=f'range({t0:.17g},tout,{t1:.17g})'
        j.study('std1').feature('time').set('tlist',expression);sol.feature('t1').set('tlist',expression)
        print('TIME BLOCK',mode,index+1,'/',len(ends),f'{t0*1e6:.3f}..{t1*1e6:.3f} us',flush=True)
        sol.runAll();history=extract_history(model,dict(s,end_time=t1))
        if prior is not None:
            seam=float(np.linalg.norm((history['values'][0]-prior['values'][-1])/weights)/max(np.linalg.norm(prior['values'][-1]/weights),1e-30))
            velocity=float(np.linalg.norm(np.asarray(sol.getUDot(jpype.JInt(1)))-expected_velocity)/max(np.linalg.norm(expected_velocity),1e-30))
            if max(seam,velocity)>1e-5:
                write_json(blocks/'failed_seam.json',dict(block=index,seam_relative_error=seam,velocity_relative_error=velocity))
                raise RuntimeError(f'Block restart continuity failed: field={seam:.3g}, velocity={velocity:.3g}.')
            state['checks'].append(dict(block=index,seam_relative_error=seam,velocity_relative_error=velocity))
        atomic_savez(blocks/f'history_{index:03d}.npz',**history)
        if frames:
            extract_movie(model,s,blocks,f'{mode}_{index:03d}',history['times'])
        # Commit a new checkpoint before publishing its pointer. A power or
        # quota interruption cannot overwrite the preceding committed state.
        checkpoint=blocks/f'restart_{index:03d}.mph'
        compact_transient(model,checkpoint,history['times'])
        state['completed']=index+1
        state['checkpoint_file']=checkpoint.name
        if s.get('adaptive_window') and t1>=s['minimum_window']:
            pieces=[dict(np.load(blocks/f'history_{i:03d}.npz')) for i in range(index+1)]
            combined=dict(times=np.concatenate([h['times'][int(i>0):] for i,h in enumerate(pieces)]),
                values=np.concatenate([h['values'][int(i>0):] for i,h in enumerate(pieces)]),y=history['y'],positions=history['positions'])
            frequencies=analysis_frequencies(s)
            full=analyze_ports([combined],frequencies,'bulk',(C,rho))
            trimmed=analyze_ports([combined],frequencies,'bulk',(C,rho),end_fraction=.9)
            delta=float(np.max(abs(full['T']-trimmed['T'])));tail=tail_fraction(combined)
            earlier=state.get('last_window_T');change=None
            if earlier is not None:
                previous=np.array(earlier['real'])+1j*np.array(earlier['imag'])
                change=float(np.max(abs(full['T']-previous)))
            stable=bool(delta<1e-4 and tail<1e-5 and change is not None and change<1e-4)
            state['stable_windows']=state.get('stable_windows',0)+1 if stable else 0
            state['window_converged']=state['stable_windows']>=2
            state['last_window_T']=dict(real=full['T'].real.tolist(),imag=full['T'].imag.tolist())
            state.setdefault('window_checks',[]).append(dict(end_time_s=t1,tail_signal_fraction=tail,
                truncated_window_T_delta=delta,previous_block_T_delta=change,passed=stable))
            print('WINDOW CHECK',state['window_checks'][-1],flush=True)
        write_json(progress,state)
    histories=[dict(np.load(blocks/f'history_{i:03d}.npz')) for i in range(state['completed'])]
    joined=dict(times=np.concatenate([h['times'][int(i>0):] for i,h in enumerate(histories)]),
                values=np.concatenate([h['values'][int(i>0):] for i,h in enumerate(histories)]),
                y=histories[0]['y'],positions=histories[0]['positions'])
    if frames:
        pieces=[dict(np.load(blocks/f'{mode}_{i:03d}_fields.npz')) for i in range(state['completed'])]
        movie={k:np.concatenate([v[k][int(i>0):] for i,v in enumerate(pieces)]) for k in ['times','u','v']}
        movie.update(x=pieces[0]['x'],y=pieces[0]['y'])
        atomic_savez(folder/f'{mode}_fields.npz',**movie)
    write_json(folder/f'{mode}_block_validation.json',dict(passed=True,blocks=state['completed'],checks=state['checks'],
        window_converged=state.get('window_converged'),window_checks=state.get('window_checks',[])))
    return joined


def run_case(client,root,kind,case,length_mm,seed,variant,args):
    cfg=json.loads((ROOT/'config.json').read_text());C,rho=validate_material(cfg)
    f0=args.bulk_frequency if kind=='bulk' else args.guided_frequency
    s=layout(kind,length_mm,f0,C,rho,variant,getattr(args,'profile','legacy'))
    if kind=='guided' and getattr(args,'native_guided_output',False):s['output_dt']=s['dt']
    if kind=='bulk' and case=='reference' and getattr(args,'profile','legacy')=='refined_v6':
        # The former fixed 12 us record failed its own tail/window checks.
        # Continue the reference until its independent window tests stabilize.
        arrival=(s['positions'][-1]-s['source_x'])/np.sqrt(C[0,0]/rho)+5/f0
        s.update(end_time=max(50e-6,4*arrival),adaptive_window=True,minimum_window=12e-6)
    if kind=='bulk':
        step_factor=getattr(args,'bulk_step_factor',1.)
        window_factor=getattr(args,'bulk_window_factor',1.)
        if step_factor>1:s.update(dt=s['dt']/step_factor,time_step_refinement=step_factor)
        if window_factor>1:s.update(end_time=s['end_time']*window_factor,window_extension_factor=window_factor)
    s['end_time']*=args.duration_factor
    if kind=='bulk' and getattr(args,'profile','legacy')=='refined_v6':
        from cfrp_pores.transient import bulk_window_limit
        s['end_time']=min(getattr(args,'bulk_window_end_time',s['end_time']),bulk_window_limit(case))
    if variant=='width' and kind=='bulk':s['width']=.002
    plate=Plate(s['width'],C[np.ix_([0,1,5],[0,1,5])],rho,s.get('safe_elements',80))
    plate.at_frequency=lru_cache(maxsize=128)(plate.at_frequency)
    if kind=='guided':
        s['mesh']=min(s['mesh'],min(plate.at_frequency(s['fmax'],m)['wavelength_m']/6 for m in MODES[kind]),s['width']/6)
        if s['fmax']>=plate.first_higher_cutoff_hz():raise ValueError('Pulse spectrum extends beyond supported A0/S0 modes.')
    pores,metadata=grouped_pores(args.diameter,length_mm,seed,args.porosity)
    if case=='reference':pores=[]
    elif case=='single':pores=[max(pores,key=lambda p:p.area)]
    if variant=='width':
        from shapely.affinity import translate
        extra=(list(pores) if getattr(args,'profile','legacy')!='legacy' else
               grouped_pores(args.diameter,length_mm,seed+100000,args.porosity)[0])
        if case=='group':pores=pores+[translate(p,yoff=1000) for p in extra]
        elif case=='single':raise ValueError('Width check applies to reference/group, not single pores.')
    folder=root/kind/variant/f'L{length_mm:g}'/f'seed_{seed}'/case;folder.mkdir(parents=True,exist_ok=True)
    signature=dict(schema=SCHEMA,layout=s,material=cfg,case=case,seed=seed,
                   geometry_hash=hashlib.sha256(b''.join(p.wkb for p in pores)).hexdigest(),frames=args.frames)
    fingerprint=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
    recordpath=folder/'run.json';previous=None
    if recordpath.exists():
        record=json.loads(recordpath.read_text(encoding='utf-8'))
        if record['fingerprint']!=fingerprint:
            old=dict(record['signature']);new=dict(signature)
            old_frames=old.pop('frames',False);new.pop('frames',None)
            presentation_only=old==new and (old_frames or not args.frames or case=='reference')
            if not presentation_only:raise ValueError('Cached settings differ; use a new output directory.')
            if args.frames and not old_frames:print('REUSE reference histories; original reference has no movie frames:',folder,flush=True)
        if record.get('status')=='complete':return record
        candidate=Path(record.get('active_recovery_directory',''))
        if (candidate/'recovery.xml').is_file():ensure_recovery_stopped(candidate)
        previous=record
    write_json(folder/'geometry.json',dict(unit='um',size_um=[length_mm*1000,s['width']*1e6],
               polygons=[dict(exterior=np.asarray(p.exterior.coords).tolist(),interiors=[np.asarray(r.coords).tolist() for r in p.interiors]) for p in pores]))
    record=dict(fingerprint=fingerprint,signature=signature,status='running',runs=[],pore_metadata=metadata)
    write_json(recordpath,record)
    for mode in MODES[kind]:
        raw=folder/f'{mode}_history.npz'
        checkpoint_frequencies=analysis_frequencies(s) if s.get('profile') else None
        if raw.exists() and (folder/f'{mode}_solve.json').exists():
            info=json.loads((folder/f'{mode}_solve.json').read_text(encoding='utf-8-sig'))
            if not (folder/f'{mode}_harmonic.npz').exists() or (s.get('profile') and not (folder/f'{mode}_harmonic_band.npz').exists()):
                resumed=client.load(str((folder/f'{mode}.mph').resolve()))
                try:frequency_checkpoint(resumed,s,folder,mode,checkpoint_frequencies)
                finally:client.remove(resumed)
                info['frequency_checkpoint_complete']=True;write_json(folder/f'{mode}_solve.json',info)
            record['runs'].append(info);continue
        model=None
        try:
            block_checkpoint=folder/'blocks'/mode/'restart.mph'
            if s.get('block_duration'):
                progress_path=block_checkpoint.parent/'progress.json'
                if progress_path.exists():
                    saved_progress=json.loads(progress_path.read_text(encoding='utf-8-sig'))
                    checkpoint_name=saved_progress.get('checkpoint_file','restart.mph')
                    if Path(checkpoint_name).name!=checkpoint_name:raise RuntimeError('Invalid checkpoint filename.')
                    block_checkpoint=block_checkpoint.with_name(checkpoint_name)
                    if saved_progress['completed'] and not block_checkpoint.exists():raise RuntimeError('Committed block checkpoint is missing.')
                model=client.load(str(block_checkpoint.resolve())) if block_checkpoint.exists() else None
            else:model=recover_mode(client,previous,mode) if previous else None
            recovered=model is not None
            if recovered:model.java.param().set('tend',f'{s["end_time"]}[s]')
            print('RESUME' if recovered else 'BUILD',kind,case,variant,'L',length_mm,mode,flush=True)
            minimum=min(s['mesh']/30,min((np.sqrt(4*p.area/np.pi)*1e-6/(60 if variant=='mesh' else 30) for p in pores),default=s['mesh']/30))
            if not recovered:
                model=create_2d(client,'transient_'+mode,cfg,pores,size_m=(s['length'],s['width']),
                     pore_offset_m=(s['start'],0),mesh_m=s['mesh'],mesh_min_m=minimum,
                     guided=kind=='guided',lateral_periodic=kind=='bulk',
                     sponge_start_m=s['sponge'] if kind=='guided' or 'left_sponge_end' in s else None,repair_tolerance_m=1e-9,
                     **builder_options(s))
                write_json(folder/f'{mode}_boundary_selections.json',getattr(model,'_boundary_selection_checks',{}))
            geometry=model.java.component('comp1').geom('geom1')
            if int(geometry.getNFiniteVoids())!=len(pores):raise RuntimeError('CAD merged or lost pores.')
            measure=geometry.measureFinal();measure.selection().geom('geom1',2);measure.selection().all()
            requested=sum(p.area for p in pores)*1e-12
            measured=s['length']*s['width']-float(measure.getArea())
            if abs(measured-requested)/(s['roi_length']*s['width'])>cfg['geometry_porosity_absolute_tolerance']:
                raise RuntimeError('Final CAD porosity differs from requested geometry.')
            record['requested_porosity']=requested/(s['roi_length']*s['width'])
            record['measured_porosity']=measured/(s['roi_length']*s['width'])
            if not recovered:source_and_study(model,s,plate,mode)
            record['active_mode']=mode
            record['active_recovery_directory']=str(client.java.getRecoveryDir(model.java.tag()))
            write_json(recordpath,record)
            print('SOLVE',mode,'elements',int(model.java.component('comp1').mesh('mesh1').getNumElem()),
                  'internal steps',int(np.ceil(s['end_time']/s['dt'])),flush=True)
            start=time.perf_counter()
            if s.get('block_duration'):
                history=solve_in_blocks(model,s,folder,mode,args.frames)
            elif recovered:
                before=np.asarray(model.evaluate('t')).ravel()
                target=np.floor(s['end_time']/s['output_dt']+1e-8)*s['output_dt']
                if before[-1]<target-s['output_dt']/2:model.java.sol('sol1').continueRun()
                after=np.asarray(model.evaluate('t')).ravel()
                if abs(after[0])>1e-15 or after[-1]<target-s['output_dt']/2:raise RuntimeError('Recovery did not preserve the full requested time history.')
            else:model.solve()
            elapsed=time.perf_counter()-start
            print('SOLVED',kind,case,variant,mode,round(elapsed,1),'seconds',flush=True)
            write_json(folder/f'{mode}_solve_progress.json',dict(stage='solved',solve_seconds=elapsed,
                       recovery_directory=record['active_recovery_directory']))
            if not s.get('block_duration'):history=extract_history(model,s)
            atomic_savez(raw,**history)
            if args.frames and not s.get('block_duration'):extract_movie(model,s,folder,mode,history['times'])
            # Keep representative actual FE states for COMSOL playback, independently
            # of the portable movie grid and the complete port time histories.
            if args.frames and not s.get('block_duration'):
                # Movie grids retain the full display cadence. Limit the editable
                # MPH snapshots to avoid holding hundreds of full FE vectors twice.
                count=min(96,len(history['times']))
                info=retain_transient_frames(model,folder/f'{mode}.mph',history['times'],count=count)
            else:
                saved_times=np.asarray(model.evaluate('t')).ravel()
                info=compact_transient(model,folder/f'{mode}.mph',saved_times)
                info['receiver_history_frames']=len(history['times'])
            check=extract_history(model,dict(s,end_time=float(history['times'][-1])))
            indices=info.get('retained_time_indices',[-1])
            np.testing.assert_allclose(check['values'],history['values'][indices],rtol=1e-8,atol=1e-12)
            info.update(mode=mode,solve_seconds=elapsed,dofs=len(model.java.sol('sol1').getU()),
                        recovered_from_checkpoint=recovered,
                        mesh_elements=int(model.java.component('comp1').mesh('mesh1').getNumElem()),
                        movie_frames=len(np.load(folder/f'{mode}_fields.npz')['times']) if args.frames else 0,
                        tail_energy_fraction=tail_fraction(history))
            write_json(folder/f'{mode}_solve.json',info)
            frequency_checkpoint(model,s,folder,mode,checkpoint_frequencies)
            info['frequency_checkpoint_complete']=True
            write_json(folder/f'{mode}_solve.json',info);record['runs'].append(info)
            write_json(recordpath,record)
        except Exception as error:
            record.update(status='failed',error=str(error));write_json(recordpath,record);raise
        finally:
            if model is not None:client.remove(model)
    record['status']='complete';write_json(recordpath,record)
    return record


def process_case(folder,bank_cache):
    record=json.loads((folder/'run.json').read_text(encoding='utf-8'));s=record['signature']['layout']
    kind=s['kind'];cfg=record['signature']['material'];C,rho=validate_material(cfg)
    histories=[dict(np.load(folder/f'{mode}_history.npz')) for mode in MODES[kind]]
    # Resolve the central useful band; the tails remain in stored raw histories.
    frequencies=analysis_frequencies(s)
    bank=None
    if kind=='guided':
        key=(s['f0'],s['width'],s.get('safe_elements',80))
        if key not in bank_cache:
            plate=Plate(s['width'],C[np.ix_([0,1,5],[0,1,5])],rho,s.get('safe_elements',80))
            bank_cache[key]=modal_bank(plate,frequencies,histories[0]['y'])
        bank=bank_cache[key]
    result=analyze_ports(histories,frequencies,kind,(C,rho),bank)
    trimmed=analyze_ports(histories,frequencies,kind,(C,rho),bank,end_fraction=.9)
    result['window_T_delta']=abs(result['T']-trimmed['T'])
    if all((folder/f'{m}_harmonic_band.npz').exists() for m in MODES[kind]):
        band=[dict(np.load(folder/f'{m}_harmonic_band.npz')) for m in MODES[kind]]
        for data in band:np.testing.assert_allclose(data['frequencies'],frequencies,rtol=1e-12)
        harmonic_band=analyze_ports(histories,frequencies,kind,(C,rho),bank,spectra_override=[h['values'] for h in band])
        result['harmonic_band_T']=harmonic_band['T'];result['harmonic_band_R']=harmonic_band['R']
        result['frequency_band_magnitude_delta']=abs(abs(result['T'])-abs(harmonic_band['T']))
    if all((folder/f'{m}_harmonic.npz').exists() for m in MODES[kind]):
        spectra=[np.load(folder/f'{m}_harmonic.npz')['values'][None] for m in MODES[kind]]
        harmonic=analyze_ports(histories,[s['f0']],kind,(C,rho),None if bank is None else [bank[8]],spectra_override=spectra)
        result['harmonic_T']=harmonic['T'][0]
        result['harmonic_R']=harmonic['R'][0]
        result['frequency_T_delta']=abs(result['T'][8]-harmonic['T'][0])
        result['frequency_T_magnitude_delta']=abs(abs(result['T'][8])-abs(harmonic['T'][0]))
        result['frequency_T_diagonal_phase_error']=np.angle(np.diag(result['T'][8])/np.diag(harmonic['T'][0]))
    atomic_savez(folder/'spectra.npz',**result)
    row=dict(kind=kind,case=record['signature']['case'],variant=s['variant'],seed=record['signature']['seed'],
             length_m=s['roi_length'],frequency_hz=s['f0'],folder=str(folder.resolve()),
             center_index=8,tail_energy_fraction=max(r['tail_energy_fraction'] for r in record['runs']),
             residual_max=float(result['residual'].max()),right_incoming_max=float(abs(result['right']).max()),
             window_delta_max=float(result['window_T_delta'].max()),
             nominal_pore_diameter_um=record['pore_metadata']['nominal_diameter_um'],
             actual_porosity=record.get('measured_porosity',0. if record['signature']['case']=='reference' else None),
             wavelength_m=result['wavelengths'][8].tolist())
    row['recorded_end_time_s']=float(histories[0]['times'][-1])
    diameter_key='diameter_max_um' if row['case']=='single' else 'diameter_mean_um'
    diameter=record['pore_metadata'][diameter_key] if row['case']!='reference' else None
    row['representative_pore_diameter_um']=diameter
    row['pore_diameter_to_wavelength']=(diameter*1e-6/result['wavelengths'][8]).tolist() if diameter is not None else []
    row['frequency_T_delta_max']=float(result['frequency_T_delta'].max()) if 'frequency_T_delta' in result else None
    if kind=='bulk':
        shortest=minimum_inplane_speed(C[np.ix_([0,1,5],[0,1,5])],rho)/s['fmax']
        row['shortest_oblique_wavelength_at_fmax_m']=shortest
        row['elements_per_shortest_oblique_wavelength']=shortest/s['mesh']
    row['frequency_T_magnitude_delta_max']=float(result['frequency_T_magnitude_delta'].max()) if 'frequency_T_delta' in result else None
    row['frequency_diagonal_phase_error_deg']=float(np.max(abs(result['frequency_T_diagonal_phase_error']))*180/np.pi) if 'frequency_T_delta' in result else None
    # Surface displacement peak delay across complete pore-free upstream ports.
    from scipy.signal import hilbert
    arrivals=[]
    for hi,(history,mode) in enumerate(zip(histories,MODES[kind])):
        component=1 if mode=='A0' else 0
        signal=history['values'][:,:3,component,-1]
        envelope=abs(hilbert(signal,axis=0));peaks=[]
        for p in range(3):
            index=int(np.argmax(envelope[:,p]));offset=0.
            if 0<index<len(envelope)-1:
                a,b,c=envelope[index-1:index+2,p]
                if a-2*b+c!=0:offset=.5*(a-c)/(a-2*b+c)
            peaks.append(history['times'][index]+offset*(history['times'][1]-history['times'][0]))
        peaks=np.asarray(peaks)
        expected=np.sqrt(C[0,0]/rho) if kind=='bulk' else bank[8][3][hi]['group_velocity_m_s']
        measured=(history['positions'][2]-history['positions'][0])/(peaks[2]-peaks[0]) if peaks[2]>peaks[0] else None
        arrivals.append(dict(mode=mode,expected_group_velocity_m_s=float(expected),measured_group_velocity_m_s=None if measured is None else float(measured),
                             relative_error=None if measured is None else float(abs(measured/expected-1))))
    row['arrival_checks']=arrivals
    return row,result,histories


def extend_bulk_window(client,root,args):
    """Continue from the saved displacement AND velocity without replaying the pulse."""
    import jpype,shutil
    if args.duration_factor<=1 or args.case not in ['single','group']:
        raise ValueError('extend requires --case single|group and --duration-factor greater than 1.')
    base=root/'bulk/base'/f'L{args.length_mm:g}'/f'seed_{args.seed}'/args.case
    original=json.loads((base/'run.json').read_text(encoding='utf-8'))
    if original['status']!='complete':raise RuntimeError('Complete the base case before extending its observation window.')
    verification=root/'restart_check.json'
    if not verification.exists() or not json.loads(verification.read_text())['passed']:
        raise RuntimeError('The monolithic-versus-restarted integration check must pass before window extension.')
    s=dict(original['signature']['layout']);s['variant']=f'window{args.duration_factor:g}';s['end_time']*=args.duration_factor
    folder=root/'bulk'/s['variant']/f'L{args.length_mm:g}'/f'seed_{args.seed}'/args.case;folder.mkdir(parents=True,exist_ok=True)
    signature=dict(original['signature'],layout=s,frames=False,base_fingerprint=original['fingerprint'])
    fingerprint=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
    previous=None;path=folder/'run.json'
    if path.exists():
        previous=json.loads(path.read_text(encoding='utf-8'))
        if previous['fingerprint']!=fingerprint:raise ValueError('Extension settings differ from cached case.')
        if previous['status']=='complete':return previous
    prefix=dict(np.load(base/'P_history.npz'));t0=float(prefix['times'][-1])
    record=dict(fingerprint=fingerprint,signature=signature,status='running',runs=[],pore_metadata=original['pore_metadata'],
                requested_porosity=original.get('requested_porosity'),measured_porosity=original.get('measured_porosity'),
                continuation_start_s=t0,base_case=str(base.resolve()))
    model=recover_mode(client,previous,'P') if previous else None;recovered=model is not None
    if model is None:
        model=client.load(str((base/'P.mph').resolve()))
        model.java.sol('sol1').copySolution('sol_initial')
    j=model.java;sol=j.sol('sol1');initial=j.sol('sol_initial')
    initial_index=jpype.JInt(len(initial.getPVals()))
    initial_v=initial.getUDot(initial_index)
    try:
        record.update(active_mode='P',active_recovery_directory=str(client.java.getRecoveryDir(j.tag())))
        write_json(path,record);shutil.copy2(base/'geometry.json',folder/'geometry.json')
        if not recovered:
            j.param().set('tend',f'{s["end_time"]}[s]')
            expression=f'range({t0:.17g},tout,tend)'
            j.study('std1').feature('time').set('tlist',expression);sol.feature('t1').set('tlist',expression)
            variables=sol.feature('v1');variables.set('control','user');variables.set('initmethod','sol')
            variables.set('initsol','sol_initial');variables.set('solnum','last')
        print('EXTEND bulk',args.case,'from',t0,'to',s['end_time'],'seconds',flush=True)
        started=time.perf_counter()
        if recovered:sol.continueRun()
        else:sol.runAll()
        elapsed=time.perf_counter()-started;segment=extract_history(model,s)
        if not np.isclose(segment['times'][0],t0,rtol=1e-10,atol=1e-16):raise RuntimeError('Continuation reset the absolute time origin.')
        C,rho=validate_material(original['signature']['material']);scale=2*np.pi*s['f0']*np.sqrt(rho*C[0,0])
        weights=np.array([1.,1.,scale,scale])[None,:,None]
        seam=float(np.linalg.norm((segment['values'][0]-prefix['values'][-1])/weights)/max(np.linalg.norm(prefix['values'][-1]/weights),1e-30))
        velocity=float(np.linalg.norm(np.asarray(sol.getUDot(jpype.JInt(1)))-np.asarray(initial_v))/max(np.linalg.norm(np.asarray(initial_v)),1e-30))
        if max(seam,velocity)>1e-5:raise RuntimeError('Displacement/stress/velocity continuity failed at the restart seam.')
        combined=dict(times=np.r_[prefix['times'],segment['times'][1:]],
                      values=np.concatenate([prefix['values'],segment['values'][1:]],axis=0),y=prefix['y'],positions=prefix['positions'])
        atomic_savez(folder/'P_history.npz',**combined)
        # Preserve COMSOL's native copied solution store unchanged. It has no
        # attached study and must not be rebuilt with createSolution().
        info=compact_transient(model,folder/'P.mph',segment['times'])
        np.testing.assert_allclose(extract_history(model,s)['values'][-1],combined['values'][-1],rtol=1e-8,atol=1e-12)
        info.update(mode='P',solve_seconds=elapsed,receiver_history_frames=len(combined['times']),
                    movie_frames=0,tail_energy_fraction=tail_fraction(combined),
                    initial_state_relative_error=seam,initial_velocity_relative_error=velocity,
                    frequency_checkpoint_complete=True,harmonic_reused_from=str((base/'P_harmonic.npz').resolve()))
        shutil.copy2(base/'P_harmonic.npz',folder/'P_harmonic.npz');write_json(folder/'P_solve.json',info)
        record.update(status='complete',runs=[info]);write_json(path,record)
    except Exception as error:
        record.update(status='failed',error=str(error));write_json(path,record);raise
    finally:client.remove(model)
    return record


def spatial_checks(client,root,kind,args,band=False):
    """Center-frequency spatial screening, separated from the transient dt check.

    These checks use the identical linear spatial operator and absorbing model.
    They do not claim broadband transient mesh convergence.
    """
    cfg=json.loads((ROOT/'config.json').read_text());C,rho=validate_material(cfg)
    f0=args.bulk_frequency if kind=='bulk' else args.guided_frequency
    refined=getattr(args,'profile','legacy')!='legacy'
    variants=(['base','mesh','absorber']+(['width'] if kind=='bulk' else []) if refined else
              ['mesh','absorber'] if kind=='guided' else ['mesh','width'])
    for variant in variants:
        for case in ['reference','group']:
            folder=root/'spatial_checks'/kind/variant/case;folder.mkdir(parents=True,exist_ok=True)
            s=layout(kind,args.length_mm,f0,C,rho,variant,getattr(args,'profile','legacy'))
            pores,_=grouped_pores(args.diameter,args.length_mm,args.seed,args.porosity)
            if case=='reference':pores=[]
            if variant=='width':
                from shapely.affinity import translate
                s['width']=.002
                if case=='group':
                    extra=(list(pores) if refined else grouped_pores(args.diameter,args.length_mm,args.seed+100000,args.porosity)[0])
                    pores += [translate(p,yoff=1000) for p in extra]
            signature=dict(schema=SCHEMA,layout=s,material=cfg,seed=args.seed,case=case,
                           geometry_hash=hashlib.sha256(b''.join(p.wkb for p in pores)).hexdigest())
            existing=folder/'settings.json'
            def spatial_signature(value):
                if not refined:return value
                value=dict(value);value['layout']={k:v for k,v in value['layout'].items()
                    if k not in ['dt','output_dt','end_time','block_duration','adaptive_window','minimum_window']}
                return value
            if existing.exists() and spatial_signature(json.loads(existing.read_text(encoding='utf-8')))!=spatial_signature(signature):
                raise ValueError('Spatial-check settings differ; use a new output directory.')
            write_json(existing,signature)
            result_name='spectra_band.npz' if band else 'spectra.npz'
            if (folder/result_name).exists():continue
            plate=Plate(s['width'],C[np.ix_([0,1,5],[0,1,5])],rho,s.get('safe_elements',80))
            frequencies=analysis_frequencies(s) if band else np.array([s['f0']])
            for mode in MODES[kind]:
                suffix='harmonic_band' if band else 'harmonic'
                if (folder/f'{mode}_{suffix}.npz').exists():continue
                print('SPECTRAL CHECK' if band else 'SPATIAL CHECK',kind,variant,case,mode,flush=True)
                minimum=min(s['mesh']/30,min((np.sqrt(4*p.area/np.pi)*1e-6/(60 if variant=='mesh' else 30) for p in pores),default=s['mesh']/30))
                model=create_2d(client,'spatial_'+mode,cfg,pores,size_m=(s['length'],s['width']),
                       pore_offset_m=(s['start'],0),mesh_m=s['mesh'],mesh_min_m=minimum,
                       guided=kind=='guided',lateral_periodic=kind=='bulk',
                       sponge_start_m=s['sponge'] if kind=='guided' or 'left_sponge_end' in s else None,repair_tolerance_m=1e-9,
                       **builder_options(s))
                try:
                    write_json(folder/f'{mode}_boundary_selections.json',getattr(model,'_boundary_selection_checks',{}))
                    if int(model.java.component('comp1').geom('geom1').getNFiniteVoids())!=len(pores):raise RuntimeError('Spatial check lost pores.')
                    source_and_study(model,s,plate,mode);frequency_checkpoint(model,s,folder,mode,frequencies)
                finally:client.remove(model)
            harmonics=[dict(np.load(folder/f'{mode}_{suffix}.npz')) for mode in MODES[kind]]
            bank=modal_bank(plate,frequencies,harmonics[0]['y']) if kind=='guided' else None
            result=analyze_ports(harmonics,frequencies,kind,(C,rho),bank,spectra_override=[h['values'] if band else h['values'][None] for h in harmonics])
            atomic_savez(folder/result_name,**result)


def render_case(folder,histories,modes):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation,PillowWriter,FFMpegWriter
    fig,axes=plt.subplots(len(histories),2,figsize=(11,3*len(histories)),squeeze=False)
    for i,(h,mode) in enumerate(zip(histories,modes)):
        for component in range(2):
            ax=axes[i,component]
            for index,label in [(1,'incident port'),(4,'transmitted port')]:
                ax.plot(h['times']*1e6,h['values'][:,index,component,-1]*1e9,label=label)
            ax.set(xlabel='Time (us)',ylabel=f'{["u","v"][component]} (nm)',title=mode);ax.legend()
    fig.tight_layout();fig.savefig(folder/'waveforms.png',dpi=140);plt.close(fig)
    for mode in modes:
        path=folder/f'{mode}_fields.npz'
        if not path.exists():continue
        data=dict(np.load(path))
        geometry=json.loads((folder/'geometry.json').read_text(encoding='utf-8'))
        settings=json.loads((folder/'run.json').read_text(encoding='utf-8'))['signature']['layout']
        cadence=float(np.max(np.diff(data['times']))*settings['fmax'])
        write_json(folder/f'{mode}_animation_sampling.json',dict(frames=len(data['times']),
            samples_per_highest_effective_period=1/cadence,above_temporal_nyquist=cadence<.5,
            meets_recommended_eight_samples=cadence<=.125+1e-12,
            time_start_s=float(data['times'][0]),time_end_s=float(data['times'][-1])))
        video=bool(settings.get('profile'))
        movie=folder/f'{mode}_displacement.{"mp4" if video else "gif"}';snapshot=folder/f'{mode}_displacement_snapshot.png'
        if movie.exists() and snapshot.exists() and min(movie.stat().st_mtime,snapshot.stat().st_mtime)>=path.stat().st_mtime:continue
        rendered_source=folder/'render_source'/path.name
        if movie.exists() and snapshot.exists() and rendered_source.exists():
            with np.load(path) as current,np.load(rendered_source) as previous:
                same=all(np.array_equal(current[k],previous[k]) for k in ['times','x','y'])
                same=same and all(np.allclose(current[k],previous[k],rtol=1e-10,atol=1e-18,equal_nan=True) for k in ['u','v'])
            if same:
                write_json(folder/f'{mode}_animation_validation.json',dict(reused=True,numerically_matching_transient_fields=True,
                    relative_tolerance=1e-10,absolute_displacement_tolerance_m=1e-18))
                continue
        if cadence>=.5:
            write_json(folder/f'{mode}_animation_status.json',dict(rendered=False,reason='Saved field cadence is below temporal Nyquist; do not animate aliased carrier waves.'))
            continue
        fig,axes=plt.subplots(2,1,figsize=(10,4));images=[]
        for ax,component in zip(axes,['u','v']):
            a=data[component]*1e9;scale=max(float(np.nanmax(abs(a))),1e-12)
            im=ax.imshow(a[0],extent=[data['x'][0]*1e3,data['x'][-1]*1e3,0,data['y'][-1]*1e3],
                         origin='lower',aspect='auto',cmap='RdBu_r',vmin=-scale,vmax=scale)
            ax.set(xlabel='x (mm)',ylabel='y (mm)');fig.colorbar(im,ax=ax,label=component+' (nm)');images.append(im)
            for polygon in geometry['polygons']:
                points=np.asarray(polygon['exterior'])/1000
                ax.plot(points[:,0]+settings['start']*1000,points[:,1],color='black',linewidth=.4)
        title=fig.suptitle('')
        def update(index):
            for im,c in zip(images,['u','v']):im.set_data(data[c][index]*1e9)
            title.set_text(f'{mode}: {data["times"][index]*1e6:.2f} us; illustrative material; thickness magnified')
            return images+[title]
        if video:
            media=Path('E:/comsol_runtime/media_tools')
            if str(media) not in sys.path:sys.path.insert(0,str(media))
            import imageio_ffmpeg
            plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
            writer=FFMpegWriter(fps=30,codec='libx264',extra_args=['-threads','1','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart'])
        else:writer=PillowWriter(fps=30)
        fig.tight_layout();FuncAnimation(fig,update,frames=len(data['times']),interval=33).save(movie,writer=writer)
        if video:
            reader=imageio_ffmpeg.read_frames(str(movie));next(reader);count=sum(1 for _ in reader)
            if count!=len(data['times']):raise RuntimeError('Encoded animation lost numerical time frames.')
            write_json(folder/f'{mode}_video_validation.json',dict(frames=count,full_decode_verified=True,source='Actual COMSOL transient displacement grid'))
        middle=int(np.argmin(abs(data['x']-(settings['start']+settings['roi_length']/2))))
        peak=int(np.argmax(np.nanmax(data['u'][:,:,middle]**2+data['v'][:,:,middle]**2,axis=1)))
        update(peak);fig.savefig(folder/f'{mode}_displacement_snapshot.png',dpi=140);plt.close(fig)


def report(root,render=False):
    rows=[];results={};cache={}
    for path in sorted(root.glob('*/*/L*/seed_*/*/run.json')):
        if json.loads(path.read_text(encoding='utf-8')).get('status')!='complete':continue
        print('ANALYZE',path.parent,flush=True)
        row,r,h=process_case(path.parent,cache);rows.append(row)
        key=(row['kind'],row['variant'],row['length_m'],row['seed'],row['case']);results[key]=r
        if render and row['variant']=='base':render_case(path.parent,h,MODES[row['kind']])
    for row in rows:
        reference_variant='base' if row['variant'].startswith('window') else row['variant']
        ref=results.get((row['kind'],reference_variant,row['length_m'],row['seed'],'reference'))
        if ref is None:continue
        row['reference_variant_used']=reference_variant
        r=results[(row['kind'],row['variant'],row['length_m'],row['seed'],row['case'])]
        if not np.array_equal(r['frequencies'],ref['frequencies']):raise ValueError('Case and reference have different frequency grids.')
        normalized=np.diagonal(r['T'],axis1=1,axis2=2)/np.diagonal(ref['T'],axis1=1,axis2=2)
        atomic_savez(Path(row['folder'])/'normalized.npz',frequencies=r['frequencies'],T=normalized,
                            R=r['R'],delta_R=r['R']-ref['R'],transmission_matrix=r['T'],power_balance=r['power_balance'])
        row.update(center_T_real=normalized[8].real.tolist(),center_T_imag=normalized[8].imag.tolist(),
                   center_amplitude_drop=(1-abs(normalized[8])).tolist(),
                   center_reflection_amplitude=abs(r['R'][8]).tolist(),
                   center_transmission_amplitude=abs(r['T'][8]).tolist(),
                   center_power_balance=r['power_balance'][8].tolist())
        weights=r['power_flux'][8,:,None]/r['power_flux'][8,None,:]
        row['center_transmission_power']=(weights*abs(r['T'][8])**2).tolist()
        row['center_reflection_power']=(weights*abs(r['R'][8])**2).tolist()
    write_json(root/'summary.json',rows)
    window_checks=[]
    for row in rows:
        if not row['variant'].startswith('window'):continue
        prior=next((v for v in rows if v['kind']==row['kind'] and v['variant']=='base' and v['case']==row['case'] and v['seed']==row['seed'] and v['length_m']==row['length_m']),None)
        if prior is None or 'center_T_real' not in row or 'center_T_real' not in prior:continue
        a=np.array(prior['center_T_real'])+1j*np.array(prior['center_T_imag'])
        b=np.array(row['center_T_real'])+1j*np.array(row['center_T_imag'])
        window_checks.append(dict(kind=row['kind'],case=row['case'],base_end_time_s=prior['recorded_end_time_s'],
            extended_end_time_s=row['recorded_end_time_s'],center_complex_T_change=abs(a-b).tolist(),
            center_magnitude_T_change=abs(abs(a)-abs(b)).tolist(),extended_tail_energy_fraction=row['tail_energy_fraction'],
            extended_window_delta_max=row['window_delta_max'],
            passed=bool(max(abs(a-b))<2e-4 and row['tail_energy_fraction']<1e-4 and row['window_delta_max']<2e-4)))
    write_json(root/'window_validation.json',window_checks)
    lines=['# 固体力学瞬态对照结果','',
           '示例补全 CFRP 材料；二维孔腔面外延伸；不代表真实铺层或实测吸收。',
           '孔形来自高斯随机场轮廓库，尺寸归一化后独立重排；不声称保留原随机场位置自相关。单孔与孔群的总孔隙率不同。',
           '接收历史含 u、v、σxx、σxy；导波按频率联合分离 S0/A0 前后行波并使用双入射矩阵归一化。',
           '动画来自真实瞬态求解帧，逐例采样信息见 *_animation_sampling.json。默认保存目标为最高有效频率每周期至少 8 帧；早期缓存可能较稀疏。低于 Nyquist 的帧序列拒绝生成动画；满足 Nyquist 但低于每周期 8 帧的回放仅作定性观察。分段求解的 MPH 保留最终位移和速度作为续算状态，完整接收历史及动画采样场另存 NPZ；未分段的代表模型最多保留 96 个有限元快照。','',
           '尾部比例是各接收截面最后 10% 时段的位移信号平方积分占比，取最大值；它不是总机械能或材料吸收率。',
           '|类型|案例|变体|孔段 mm|模态|中心频率 MHz|位移幅值下降 %|尾部信号比例|',
           '|---|---|---|---:|---|---:|---:|---:|']
    for row in rows:
        for mode,drop in zip(MODES[row['kind']],row.get('center_amplitude_drop',[])):
            lines.append(f'|{row["kind"]}|{row["case"]}|{row["variant"]}|{row["length_m"]*1e3:g}|{mode}|{row["frequency_hz"]/1e6:g}|{drop*100:.6g}|{row["tail_energy_fraction"]:.3g}|')
    lines+=['','仅有单个孔隙段长度或单个随机样本时不生成衰减系数及统计区间。',
            '体波仅提取截面平均 P 通道；散射到横向周期空间谐波及剪切波的能量不能当作吸收。',
            '不同频率的体波／导波不以总幅值下降直接排序灵敏度。数值验收见 validation.json。']
    for check in window_checks:
        lines.append(f'\n观察窗口由 {check["base_end_time_s"]*1e6:.4g} μs 延长至 {check["extended_end_time_s"]*1e6:.4g} μs；中心频率幅值变化为 {[round(x*100,5) for x in check["center_magnitude_T_change"]]} 个百分点；窗口验收：{check["passed"]}。')
    (root/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    if render:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,2,figsize=(12,7),squeeze=False)
        for column,kind in enumerate(['bulk','guided']):
            for row in rows:
                if row['kind']!=kind or row['variant']!='base' or row['case']=='reference' or 'center_T_real' not in row:continue
                d=np.load(Path(row['folder'])/'normalized.npz')
                raw=results[(kind,'base',row['length_m'],row['seed'],row['case'])]
                power=raw['power_flux'];weight=power[:,:,None]/power[:,None,:]
                for mi,mode in enumerate(MODES[kind]):
                    label=mode+' '+row['case']
                    axes[0,column].plot(d['frequencies']/1e6,100*(1-abs(d['T'][:,mi])),label=label)
                    reflected=np.sum(weight*abs(raw['R'])**2,axis=1)[:,mi]
                    axes[1,column].plot(d['frequencies']/1e6,100*reflected,label=label+' reflection')
                    if kind=='guided':
                        converted=weight[:,1-mi,mi]*abs(raw['T'][:,1-mi,mi])**2
                        axes[1,column].plot(d['frequencies']/1e6,100*converted,'--',label=label+' conversion')
            axes[0,column].set(title=kind,xlabel='Frequency (MHz)',ylabel='Normalized amplitude reduction (%)')
            axes[1,column].set(xlabel='Frequency (MHz)',ylabel='Incident modal power (%)')
            for ax in axes[:,column]:
                ax.grid(alpha=.25)
                if ax.lines:ax.legend(fontsize=7)
        fig.suptitle('Matched pore statistics; different frequencies and wavelengths')
        fig.tight_layout();fig.savefig(root/'transmission_reflection.png',dpi=160);plt.close(fig)
        fig,axes=plt.subplots(1,3,figsize=(12,4))
        for ax,(kind,mode,index) in zip(axes,[('bulk','P',0),('guided','S0',0),('guided','A0',1)]):
            selected=[v for v in rows if v['kind']==kind and v['variant']=='base' and v['case']=='group' and 'center_amplitude_drop' in v]
            if not selected:ax.set_visible(False);continue
            group=min(selected,key=lambda v:(v['length_m'],v['seed']))
            single=next((v for v in rows if v['kind']==kind and v['variant']=='base' and v['case']=='single' and v['seed']==group['seed'] and v['length_m']==group['length_m']),None)
            have_single=single is not None and 'center_amplitude_drop' in single
            values=[100*single['center_amplitude_drop'][index] if have_single else 0.,100*group['center_amplitude_drop'][index]]
            bars=ax.bar(['Single pore' if have_single else 'No pores','Pore group'],values,color=['#7c9cbf','#d78b42'])
            ax.bar_label(bars,labels=[f'{v:.3g}%' for v in values],padding=4)
            quality=(f'{group["recorded_end_time_s"]*1e6:.2f} us window; '
                     +('unconverged' if group['tail_energy_fraction']>=1e-4 else 'see validation')) if mode=='P' else 'One realization'
            ax.set(title=f'{mode}: {group["frequency_hz"]/1e6:g} MHz\n{quality}',ylabel='Transmission amplitude reduction (%)',
                   xlabel=f'Wavelength {group["wavelength_m"][index]*1000:.3g} mm; one realization')
            ax.margins(y=.22);ax.grid(axis='y',alpha=.2)
        fig.suptitle('Finite-window transmission comparison: illustrative material')
        fig.text(.5,.01,'Group area porosity: 1.5%. See the validation report for convergence and statistical limits.',ha='center',fontsize=9)
        fig.tight_layout(rect=(0,.04,1,.95));fig.savefig(root/'comparison_summary.png',dpi=160);plt.close(fig)
        for check in window_checks:
            candidates=[v for v in rows if v['kind']==check['kind'] and v['case']==check['case'] and v['variant'] in ['base','window2']]
            if check['kind']!='bulk' or len(candidates)!=2:continue
            fig,axes=plt.subplots(2,1,figsize=(10,6))
            for row in candidates:
                folder=Path(row['folder']);history=np.load(folder/'P_history.npz');normalized=np.load(folder/'normalized.npz')
                label=f'{row["recorded_end_time_s"]*1e6:.2f} us record'
                mean_u=np.trapezoid(history['values'][:,4,0],history['y'],axis=-1)/np.ptp(history['y'])
                axes[0].plot(history['times']*1e6,mean_u*1e9,label=label,alpha=.85)
                axes[1].plot(normalized['frequencies']/1e6,100*(1-abs(normalized['T'][:,0])),label=label)
            axes[0].axvline(check['base_end_time_s']*1e6,color='gray',linestyle='--')
            axes[0].set(xlabel='Time (us)',ylabel='Mean transmitted u (nm)')
            axes[1].set(xlabel='Frequency (MHz)',ylabel='Transmission amplitude reduction (%)')
            for ax in axes:ax.legend();ax.grid(alpha=.2)
            fig.suptitle('Bulk P: observation window is still insufficient; quantitative loss is unconverged')
            fig.tight_layout();fig.savefig(root/'observation_window_comparison.png',dpi=160);plt.close(fig)
    return rows,results


def validation(root,rows,results):
    checks=[]
    for row in rows:
        if row['variant']!='base':continue
        kind,L,seed,case=[row[k] for k in ['kind','length_m','seed','case']]
        base=results[(kind,'base',L,seed,case)]
        refined=json.loads((Path(row['folder'])/'run.json').read_text(encoding='utf-8'))['signature']['layout'].get('profile')
        required=['mesh','time','absorber'] if kind=='guided' else (['mesh','time','absorber','width'] if refined else ['mesh','time','width'])
        differences={}
        for variant in required:
            other=results.get((kind,variant,L,seed,case));ref=results.get((kind,variant,L,seed,'reference'))
            bref=results.get((kind,'base',L,seed,'reference'))
            method='transient'
            spatial=root/'spatial_checks'/kind/variant
            if other is None and (spatial/case/'spectra.npz').exists() and (spatial/'reference/spectra.npz').exists():
                settings=json.loads((spatial/case/'settings.json').read_text(encoding='utf-8'))
                if settings['seed']!=seed or not np.isclose(settings['layout']['roi_length'],L):continue
                if 'harmonic_T' not in base or bref is None or 'harmonic_T' not in bref:continue
                other=dict(np.load(spatial/case/'spectra.npz'));ref=dict(np.load(spatial/'reference/spectra.npz'))
                a0=base['harmonic_T'];ar=bref['harmonic_T'];b0=other['T'][0];br=ref['T'][0]
                method='matched center-frequency spatial screening'
            else:
                if other is None or ref is None or bref is None:continue
                a0=base['T'][8];ar=bref['T'][8];b0=other['T'][8];br=ref['T'][8]
            if case=='reference':a=np.diag(a0);b=np.diag(b0)
            else:a=np.diag(a0)/np.diag(ar);b=np.diag(b0)/np.diag(br)
            differences[variant]=dict(complex_T_delta=float(max(abs(a-b))),
                magnitude_delta=float(max(abs(abs(a)-abs(b)))),phase_error_deg=float(max(abs(np.angle(a/b)))*180/np.pi),
                method=method,
                by_mode={mode:dict(magnitude_delta=float(abs(abs(aa)-abs(bb))),phase_error_deg=float(abs(np.angle(aa/bb))*180/np.pi),
                    relative_loss_change=float(abs(abs(aa)-abs(bb))/max(abs(1-abs(aa)),abs(1-abs(bb)),1e-12)))
                    for mode,aa,bb in zip(MODES[kind],a,b)},
                relative_loss_change=float(max(abs(abs(a)-abs(b))/np.maximum(np.maximum(abs(1-abs(a)),abs(1-abs(b))),1e-12))))
            if refined:
                if variant=='time' and other is not None and len(other['T'])==len(base['T']):
                    aa=np.diagonal(base['T'],axis1=1,axis2=2);bb=np.diagonal(other['T'],axis1=1,axis2=2)
                    if case!='reference':
                        aa=aa/np.diagonal(bref['T'],axis1=1,axis2=2);bb=bb/np.diagonal(ref['T'],axis1=1,axis2=2)
                elif all((root/'spatial_checks'/kind/v/c/'spectra_band.npz').exists() for v in ['base',variant] for c in [case,'reference']):
                    def band_values(v,c):return np.load(root/'spatial_checks'/kind/v/c/'spectra_band.npz')['T']
                    aa=np.diagonal(band_values('base',case),axis1=1,axis2=2)
                    bb=np.diagonal(band_values(variant,case),axis1=1,axis2=2)
                    if case!='reference':
                        aa=aa/np.diagonal(band_values('base','reference'),axis1=1,axis2=2)
                        bb=bb/np.diagonal(band_values(variant,'reference'),axis1=1,axis2=2)
                else:continue
                differences[variant]['band']=dict(magnitude_delta=float(np.max(abs(abs(aa)-abs(bb)))),
                    phase_error_deg=float(np.max(abs(np.angle(aa/bb)))*180/np.pi),
                    relative_loss_change=float(np.max(abs(abs(aa)-abs(bb))/np.maximum(np.maximum(abs(1-abs(aa)),abs(1-abs(bb))),1e-12))))
        complete=all(v in differences for v in required)
        passed=complete and all(v['magnitude_delta']<=2e-4 and v['phase_error_deg']<=1. and (case=='reference' or v['relative_loss_change']<=.05) for v in differences.values())
        if refined:
            complete=complete and all('band' in d for d in differences.values())
            passed=passed and complete and all(d['band']['magnitude_delta']<=2e-4 and d['band']['phase_error_deg']<=1 and
                (case=='reference' or d['band']['relative_loss_change']<=.05) for d in differences.values())
        passed=passed and row['tail_energy_fraction']<1e-4 and row['window_delta_max']<2e-4
        passed=passed and row['frequency_T_magnitude_delta_max'] is not None and row['frequency_T_magnitude_delta_max']<2e-4
        passed=passed and row['frequency_diagonal_phase_error_deg'] is not None and row['frequency_diagonal_phase_error_deg']<1.
        frequency_band_consistency=None
        if refined:
            available='harmonic_band_T' in base and bref is not None and 'harmonic_band_T' in bref
            complete=complete and available;passed=passed and available
            if available:
                aa=np.diagonal(base['T'],axis1=1,axis2=2);bb=np.diagonal(base['harmonic_band_T'],axis1=1,axis2=2)
                if case!='reference':
                    aa=aa/np.diagonal(bref['T'],axis1=1,axis2=2)
                    bb=bb/np.diagonal(bref['harmonic_band_T'],axis1=1,axis2=2)
                magnitude=float(np.max(abs(abs(aa)-abs(bb))));phase=float(np.max(abs(np.angle(aa/bb)))*180/np.pi)
                relative=float(np.max(abs(abs(aa)-abs(bb))/np.maximum(np.maximum(abs(1-abs(aa)),abs(1-abs(bb))),1e-12)))
                frequency_band_consistency=dict(magnitude_delta=magnitude,phase_error_deg=phase,relative_loss_change=relative,
                    passed=bool(magnitude<=2e-4 and phase<=1 and (case=='reference' or relative<=.05)),
                    by_mode={mode:dict(magnitude_delta=float(np.max(abs(abs(aa[:,i])-abs(bb[:,i])))),
                        relative_loss_change=float(np.max(abs(abs(aa[:,i])-abs(bb[:,i]))/np.maximum(np.maximum(abs(1-abs(aa[:,i])),abs(1-abs(bb[:,i]))),1e-12))),
                        phase_error_deg=float(np.max(abs(np.angle(aa[:,i]/bb[:,i])))*180/np.pi))
                        for i,mode in enumerate(MODES[kind])})
                passed=passed and frequency_band_consistency['passed']
        if case=='reference':passed=passed and all(a['relative_error'] is not None and a['relative_error']<.05 for a in row['arrival_checks'])
        if kind=='guided':passed=passed and row['residual_max']<.01 and row['right_incoming_max']<.001 and float(np.max(abs(base['power_balance'][8]-1)))<.005
        checks.append(dict(kind=kind,case=case,length_m=L,seed=seed,complete=complete,passed=bool(passed),
                           differences=differences,tail_energy_fraction=row['tail_energy_fraction'],
                           window_delta_max=row['window_delta_max'],residual_max=row['residual_max'],
                           harmonic_magnitude_delta=row['frequency_T_magnitude_delta_max'],
                           harmonic_diagonal_phase_error_deg=row['frequency_diagonal_phase_error_deg'],
                           frequency_band_consistency=frequency_band_consistency,
                           arrival_checks=row['arrival_checks']))
    # All reference and group cases in both physics families must pass; missing cases cannot pass vacuously.
    needed={(k,c) for k in ['bulk','guided'] for c in ['reference','group']}
    eligible=[c for c in checks if (c['kind'],c['case']) in needed]
    ready=needed.issubset({(c['kind'],c['case']) for c in eligible}) and all(c['passed'] for c in eligible)
    value=dict(ready_for_ensemble=bool(ready),checks=checks,
               thresholds=dict(absolute_amplitude=.0002,phase_degrees=1.,relative_nonzero_loss=.05,
                               modal_residual=.01,right_incoming_amplitude=.001,tail_energy=.0001,arrival_relative_error=.05),
               note='No near-zero absolute tolerance is substituted for relative attenuation precision.')
    write_json(root/'validation.json',value);return value


def statistics(root,rows):
    output=[]
    for kind in ['bulk','guided']:
        selected=[r for r in rows if r['kind']==kind and r['variant']=='base' and r['case']=='group' and 'center_T_real' in r]
        lengths=sorted(set(r['length_m'] for r in selected));seeds=sorted(set(r['seed'] for r in selected))
        complete=[s for s in seeds if all(any(r['seed']==s and r['length_m']==L for r in selected) for L in lengths)]
        if len(lengths)<3 or len(complete)<3:continue
        for index,mode in enumerate(MODES[kind]):
            data=np.array([[next(complex(r['center_T_real'][index],r['center_T_imag'][index]) for r in selected if r['seed']==s and r['length_m']==L) for L in lengths] for s in complete])
            item=attenuation_interval(lengths,data);item.update(kind=kind,mode=mode,seeds=complete,material_constant_validated=False)
            if len(complete)>=6:
                half=attenuation_interval(lengths,data[:len(complete)//2])['fit']['alpha_np_m'];full=item['fit']['alpha_np_m']
                item['sample_count_relative_change']=abs(full-half)/max(abs(full),abs(half),1e-30)
                item['sample_stability_passed']=item['sample_count_relative_change']<=.05
            output.append(item)
    write_json(root/'attenuation.json',output)


def append_validation_report(root,gates):
    lines=['','## 验证状态','',
           '**基础检查通过，可开展长度／样本扩展。**' if gates['ready_for_ensemble'] else
           '**基础检查尚未全部通过；当前结果用于波场与透射对照，不作为收敛的材料衰减常数。**','',
           '|类型／模态|案例|检查|幅值差（百分点）|相位差（度）|损失相对变化|方法|',
           '|---|---|---|---:|---:|---:|---|']
    for check in gates['checks']:
        for name,d in check['differences'].items():
            for mode,values in d.get('by_mode',{'all':d}).items():
                relative='不适用（无孔参考）' if check['case']=='reference' else f'{100*values["relative_loss_change"]:.3g}%'
                lines.append(f'|{check["kind"]}/{mode}|{check["case"]}|{name}|{100*values["magnitude_delta"]:.6g}|{values["phase_error_deg"]:.5g}|{relative}|{d["method"]}|')
        if not check['complete']:lines.append(f'|{check["kind"]}|{check["case"]}|尚有检查未完成|—|—|—|见 validation.json|')
    lines+=['','中心频率幅值绝对差目标为 0.02 个百分点、相位差为 1°；孔隙损失的相对变化目标为 5%。',
            '上表列中心频率控制；refined 配置另要求 17 点激励频带的空间控制、时间步与瞬态／频域一致性全部通过，详见 validation.json。refined 体波宽度检查复制同一周期单元，只验证数值一致性，不证明随机代表体积收敛；legacy 配置的宽度检查增加独立孔隙样本。',
            '完整置信区间仅在三个以上孔隙段长度及三个以上完整随机种子族齐全时计算；尚未满足时 attenuation.json 为空。','',
            '## 结果可信范围','',
            '- 本报告只统计 run.json 标记完成的算例；正在计算和未运行的算例不记为完成。已有动画适合检查传播、散射和模态形状；播放速度经过放慢，厚度方向放大。',
            '- 可信范围限定为当前示例材料、几何和随机样本；中心频率的稳定性不能外推为整个频带或随机总体已收敛。A0 需关注孔边／厚度网格和吸收区，体波需关注长时散射、周期单元宽度和更小内部步长。',
            '- 表内未完成的检查仍需补齐；在基础检查通过前，不将当前结果升级为材料衰减系数或统计置信区间。',
            '- 同模态透射下降、反射和模态转换均已分列在 summary.json 与频谱图中。它们不能合并解释成材料吸收；体波还未计入全部剪切波和横向空间谐波通道。','',
            '## 当前孔群的数值敏感性','']
    rows=json.loads((root/'summary.json').read_text(encoding='utf-8'))
    for row in rows:
        if row['variant']!='base' or row['case']!='group':continue
        check=next((c for c in gates['checks'] if all(c[k]==row[k] for k in ['kind','case','length_m','seed'])),None)
        if check is None:continue
        for mode,drop in zip(MODES[row['kind']],row.get('center_amplitude_drop',[])):
            changes=[f'{name}: {100*d["by_mode"][mode]["relative_loss_change"]:.3g}%' for name,d in check['differences'].items() if mode in d.get('by_mode',{})]
            lines.append(f'- {row["kind"]}/{mode}，{row["frequency_hz"]/1e6:g} MHz：幅值下降 {100*drop:.6g}%；损失相对变化（{", ".join(changes)}）。')
    for check in gates['checks']:
        if not check['complete']:
            lines.append(f'- {check["kind"]}/{check["case"]} 尚未完成全部独立检查；不标为通过。')
    audit_path=root/'saved_model_audit_both.json'
    if audit_path.exists():
        audits=json.loads(audit_path.read_text(encoding='utf-8'))
        lines+=['',f'已重新打开 {len(audits)} 个保存模型，核对保留时刻及接收位移／应力；记录见 saved_model_audit_both.json。']
    lines+=['','## 图像与模型','']
    provenance=root/'guided/base/L8/seed_101/group/movie_field_provenance.json'
    if provenance.exists():
        items=json.loads(provenance.read_text(encoding='utf-8-sig'))['checks']
        movie_end=max(c['movie_end_time_s'] for c in items)*1e6
        receiver_end=max(c['receiver_validation_end_time_s'] for c in items)*1e6
        lines += [f'导波动画显示前 {movie_end:.4g} μs 的实际有限元位移场；接收验收延长到 {receiver_end:.4g} μs。动画帧来自相同几何、材料和内部时间步的原求解，已核对全部共有接收时刻，未伪造后续时间帧。',
                  f'[动画复用来源及核对记录]({provenance.resolve().as_posix()})','']
    for pattern,label in [('**/*_displacement.mp4','位移动画视频'),('**/*_displacement.gif','位移动画'),('**/waveforms.png','接收波形'),('comparison_summary.png','中心频率对比'),('observation_window_comparison.png','体波观察窗口对比'),('transmission_reflection.png','透射／反射频谱')]:
        for path in sorted(root.glob(pattern)):
            if '_runtime' in path.parts:continue
            lines.append(f'- [{label}：{path.relative_to(root).as_posix()}]({path.resolve().as_posix()})')
    lines+=['','本次临时与恢复文件使用 E:/comsol_runtime/transient_matched/；启动时验证实际模型恢复路径，清单见 _runtime/session_*.json。']
    with (root/'REPORT.md').open('a',encoding='utf-8') as stream:stream.write('\n'.join(lines)+'\n')


def audit_saved_models(client,root,kind='both'):
    """Reopen saved base models and compare every retained port state to raw data."""
    audits=[]
    for path in sorted(root.glob('*/*/L*/seed_*/*/run.json')):
        record=json.loads(path.read_text(encoding='utf-8'));s=record['signature']['layout']
        if record.get('status')!='complete' or kind not in ['both',s['kind']]:continue
        folder=path.parent
        for mode in MODES[s['kind']]:
            print('REOPEN CHECK',folder,mode,flush=True)
            raw=dict(np.load(folder/f'{mode}_history.npz'))
            info=json.loads((folder/f'{mode}_solve.json').read_text(encoding='utf-8-sig'))
            indices=info.get('retained_time_indices',[-1])
            model=client.load(str((folder/f'{mode}.mph').resolve()))
            try:
                actual=extract_history(model,dict(s,end_time=float(raw['times'][-1])))
                np.testing.assert_allclose(actual['times'],raw['times'][indices],rtol=1e-12,atol=1e-18)
                np.testing.assert_allclose(actual['values'],raw['values'][indices],rtol=1e-8,atol=1e-12)
                recovery=str(client.java.getRecoveryDir(model.java.tag()))
                if not Path(recovery).resolve().is_relative_to(runtime_root()):raise RuntimeError('Reopened model recovery escaped configured runtime directory.')
                audits.append(dict(model=str((folder/f'{mode}.mph').resolve()),frames=len(indices),
                                   saved_times_verified=True,all_retained_ports_verified=True,recovery_directory=recovery))
            finally:client.remove(model)
    write_json(root/f'saved_model_audit_{kind}.json',audits)
    return audits


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['pilot','spatial','spectral','verify','ensemble','report','audit','extend'])
    parser.add_argument('--output',type=Path,default=ROOT/'output/transient_matched')
    parser.add_argument('--kind',choices=['bulk','guided','both'],default='both')
    parser.add_argument('--case',choices=['reference','single','group','all'],default='all')
    parser.add_argument('--variant',choices=['base','mesh','time','absorber','width'],default='base')
    parser.add_argument('--profile',choices=['legacy','refined','refined_v2','refined_v3','refined_v4','refined_v5','refined_v6'],default='legacy')
    parser.add_argument('--length-mm',type=float,default=8.)
    parser.add_argument('--seed',type=int,default=101);parser.add_argument('--samples',type=int,default=6)
    parser.add_argument('--diameter',type=float,default=80.);parser.add_argument('--porosity',type=float,default=.015)
    parser.add_argument('--bulk-frequency',type=float,default=5e6);parser.add_argument('--guided-frequency',type=float,default=2e5)
    parser.add_argument('--duration-factor',type=float,default=1.,help='Multiply the observation duration; changed settings require a fresh output directory.')
    parser.add_argument('--native-guided-output',action='store_true',help='Save receivers every internal step for guided waves.')
    parser.add_argument('--bulk-step-factor',type=float,default=1.,help='Divide the bulk internal step by this factor; use a fresh case directory.')
    parser.add_argument('--bulk-window-factor',type=float,default=1.,help='Multiply the bulk observation limit; retain committed blocks when extending.')
    parser.add_argument('--cores',type=int,default=4);parser.add_argument('--frames',action=argparse.BooleanOptionalAction,default=True)
    parser.add_argument('--render',action='store_true');args=parser.parse_args()
    if args.length_mm<=0 or args.samples<3 or args.diameter<=0 or args.duration_factor<1 or args.bulk_step_factor<1 or args.bulk_window_factor<1 or not 0<args.porosity<.1:parser.error('Invalid geometry/sample/time parameters.')
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=True)
    configure_scratch(root)
    if args.stage=='audit':
        audit_saved_models(start_session(root,args.cores),root,args.kind);return
    if args.stage=='extend':
        if args.kind!='bulk':parser.error('Window extension currently requires --kind bulk.')
        extend_bulk_window(start_session(root,args.cores),root,args)
        rows,results=report(root,args.render);gates=validation(root,rows,results);statistics(root,rows);append_validation_report(root,gates)
        return
    if args.stage!='report':
        if args.stage=='ensemble':
            gate=root/'validation.json'
            if not gate.exists() or not json.loads(gate.read_text())['ready_for_ensemble']:
                raise RuntimeError('Basic numerical gates have not passed. Run pilot/verify and inspect validation.json first.')
        client=start_session(root,args.cores)
        for kind in (['bulk','guided'] if args.kind=='both' else [args.kind]):
            if args.stage in ['spatial','spectral','verify']:
                spatial_checks(client,root,kind,args,band=args.stage=='spectral')
                if args.stage in ['spatial','spectral']:continue
            variants=([args.variant] if args.stage!='verify' else ['time'])
            lengths=([args.length_mm] if args.stage!='ensemble' else [8.,24.,48.])
            seeds=([args.seed] if args.stage!='ensemble' else list(range(args.seed,args.seed+args.samples)))
            for variant in variants:
                if variant=='width' and kind!='bulk':raise ValueError('Width variant is for bulk only.')
                for length in lengths:
                    for seed in seeds:
                        cases=((['reference','group','single'] if args.profile!='legacy' else ['reference','single','group']) if args.stage=='pilot' else ['reference','group']) if args.case=='all' else [args.case]
                        for case in cases:run_case(client,root,kind,case,length,seed,variant,args)
    rows,results=report(root,args.render);gates=validation(root,rows,results);statistics(root,rows);append_validation_report(root,gates)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');sys.stderr.reconfigure(encoding='utf-8');main()
