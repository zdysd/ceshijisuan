"""Test one exact Newton correction for the linear bulk transient operator.

Production settings are untouched. Compare full receiver histories including
traction and, for the checkpoint test, the final FE displacement AND velocity.
"""
import argparse,json,shutil,time
import numpy as np
import run_transient as rt


def linear_iteration(model):
    import jpype
    fc=model.java.sol('sol1').feature('t1').feature('fc1')
    fc.set('dtech','const');fc.set('damp',1.)
    fc.set('ntermconst','iter');fc.set('niter',jpype.JInt(1))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full-checkpoint',action='store_true');args=parser.parse_args()
    root=rt.ROOT/'output/transient_refined/linear_bulk_check'/f'attempt_{time.time_ns()}'
    root.mkdir(parents=True);client=rt.start_session(root,2)
    cfg=json.loads((rt.ROOT/'config.json').read_text());C,rho=rt.validate_material(cfg)
    s=rt.layout('bulk',.2,5e6,C,rho,profile='refined_v6')
    s.update(width=.0002,start=.0007,roi_length=.0002,sponge=.0015,length=.002,
             absorb_length=.0005,source_x=.0003,left_sponge_end=.0002,
             positions=[.0004,.0005,.0006,.001,.0011,.0012],mesh=.0001,
             end_time=1.5e-6,port_samples=17,adaptive_window=False)
    source=None
    if args.full_checkpoint:
        folder=rt.ROOT/'output/transient_converged/bulk/base/L8/seed_101/group'
        record=json.loads((folder/'run.json').read_text());s=record['signature']['layout']
        progress=json.loads((folder/'blocks/P/progress.json').read_text())
        source=folder/'blocks/P'/progress.get('checkpoint_file','restart.mph')
        before=(source.stat().st_size,source.stat().st_mtime_ns)
        shutil.copy2(source,root/'source_checkpoint.mph')
        if before!=(source.stat().st_size,source.stat().st_mtime_ns):raise RuntimeError('Source checkpoint changed during copy; retry later.')
    histories=[];states=[];velocities=[];timings=[]
    for optimized in [False,True]:
        if args.full_checkpoint:
            model=client.load(str((root/'source_checkpoint.mph').resolve()))
            j=model.java;sol=j.sol('sol1');t0=float(np.asarray(model.evaluate('t')).ravel()[-1])
            initial=rt.extract_history(model,dict(s,end_time=t0))['values'][-1]
            import jpype
            initial_velocity=np.asarray(sol.getUDot(jpype.JInt(1))).copy()
            sol.copySolution('linear_test_initial')
            v=sol.feature('v1');v.set('control','user');v.set('initmethod','sol');v.set('initsol','linear_test_initial');v.set('solnum','last')
            sol.feature('t1').set('consistent','off')
            end=t0+10*s['output_dt']
            expression=f'range({t0:.17g},tout,{end:.17g})'
            j.study('std1').feature('time').set('tlist',expression);sol.feature('t1').set('tlist',expression)
        else:
            model=rt.create_2d(client,'linear_bulk_test',cfg,size_m=(s['length'],s['width']),mesh_m=s['mesh'],
                 lateral_periodic=True,guided=False,sponge_start_m=s['sponge'],**rt.builder_options(s))
            rt.source_and_study(model,s,None,'P');end=s['end_time']
        try:
            if optimized:linear_iteration(model)
            print('SOLVE',dict(full_checkpoint=args.full_checkpoint,one_iteration=optimized),flush=True)
            started=time.perf_counter()
            if args.full_checkpoint:model.java.sol('sol1').runAll()
            else:model.solve()
            timings.append(time.perf_counter()-started)
            h=rt.extract_history(model,dict(s,end_time=end));histories.append(h)
            import jpype
            sol=model.java.sol('sol1');states.append(np.asarray(sol.getU(jpype.JInt(len(h['times'])))).copy())
            velocities.append(np.asarray(sol.getUDot(jpype.JInt(len(h['times'])))).copy())
            if args.full_checkpoint:
                np.testing.assert_allclose(h['values'][0],initial,rtol=1e-8,atol=1e-15)
                np.testing.assert_allclose(np.asarray(sol.getUDot(jpype.JInt(1))),initial_velocity,rtol=1e-8,atol=1e-15)
            np.savez_compressed(root/('one_iteration.npz' if optimized else 'default.npz'),**h)
        finally:client.remove(model)
    np.testing.assert_allclose(histories[0]['times'],histories[1]['times'],rtol=1e-12,atol=1e-17)
    scale=np.array([1,1,2*np.pi*s['f0']*np.sqrt(rho*C[0,0]),2*np.pi*s['f0']*np.sqrt(rho*C[0,0])])[None,None,:,None]
    a=histories[0]['values']/scale;b=histories[1]['values']/scale
    history_error=float(np.linalg.norm(a-b)/np.linalg.norm(a))
    maximum_error=float(abs(a-b).max()/abs(a).max())
    state_error=float(np.linalg.norm(states[1]-states[0])/max(np.linalg.norm(states[0]),1e-100))
    velocity_error=float(np.linalg.norm(velocities[1]-velocities[0])/max(np.linalg.norm(velocities[0]),1e-100))
    result=dict(full_checkpoint=args.full_checkpoint,history_relative_error=history_error,maximum_scaled_error=maximum_error,
        final_state_relative_error=state_error,final_velocity_relative_error=velocity_error,tolerance=1e-10,
        passed=max(history_error,maximum_error,state_error,velocity_error)<1e-10,
        default_and_one_iteration_seconds=timings,speed_ratio=timings[0]/timings[1],source_checkpoint=str(source),
        run_folder=str(root),reference='COMSOL 6.2 FullyCoupled: dtech=const, ntermconst=iter, niter=1; only for this linear operator.')
    rt.write_json(root/'validation.json',result)
    rt.write_json(root.parent/('full_checkpoint_validation.json' if args.full_checkpoint else 'pulse_validation.json'),result)
    print(json.dumps(result),flush=True)
    if not result['passed']:raise RuntimeError('One-iteration linear solve changed the computed fields.')


if __name__=='__main__':main()
