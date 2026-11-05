"""Compare the refined block restart against one uninterrupted COMSOL solve."""
import run_transient as rt
import argparse,json,time,numpy as np
from pathlib import Path


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--kind',choices=['bulk','guided'],default='bulk')
    parser.add_argument('--order',type=int,choices=[3,4],default=3)
    parser.add_argument('--extension-check',action='store_true')
    args=parser.parse_args()
    base=rt.ROOT/'output/transient_refined'/('restart_smoke' if args.kind=='bulk' else 'restart_smoke_guided')
    root=base/f'attempt_{time.time_ns()}'
    client=rt.start_session(base,2)
    cfg=json.loads((rt.ROOT/'config.json').read_text());C,rho=rt.validate_material(cfg)
    s=rt.layout('bulk',.2,5e6,C,rho,profile='refined')
    s.update(width=.0002,start=.0007,roi_length=.0002,sponge=.0015,length=.002,
             absorb_length=.0005,source_x=.0003,left_sponge_end=.0002,
             positions=[.0004,.0005,.0006,.001,.0011,.0012],mesh=.0001,
             end_time=1.5e-6,block_duration=.5e-6,port_samples=17)
    mode='P';plate=None
    if args.kind=='guided':
        s=rt.layout('guided',1.,2e5,C,rho,profile='refined')
        s.update(width=.0002,start=.001,roi_length=.001,sponge=.003,length=.004,
                 absorb_length=.001,positions=[.0002,.0004,.0006,.0022,.0024,.0026],mesh=.0002,
                 end_time=30e-6,block_duration=10e-6,port_samples=17)
        mode='A0';plate=rt.Plate(s['width'],C[np.ix_([0,1,5],[0,1,5])],rho,20)
    s['displacement_order']=args.order
    if args.kind=='bulk' and args.order==4:
        s.update(periodic_mesh_match=True,dt=1/(200*s['fmax']))
    outputs=[];timings=[]
    for blocked in [False,True]:
        m=rt.create_2d(client,'restart_test',cfg,size_m=(s['length'],s['width']),
              mesh_m=s['mesh'],lateral_periodic=args.kind=='bulk',guided=args.kind=='guided',
              sponge_start_m=s['sponge'],**rt.builder_options(s))
        try:
            # Compare default automatic assembly to explicit constant-matrix
            # reuse AND block restarts in one complete-history check.
            rt.source_and_study(m,dict(s,profile=None) if not blocked else s,plate,mode)
            started=time.perf_counter()
            if blocked:
                if args.extension_check:
                    rt.solve_in_blocks(m,dict(s,end_time=.87*s['end_time']),root,mode,False)
                    progress=json.loads((root/'blocks'/mode/'progress.json').read_text())
                    client.remove(m)
                    m=client.load(str((root/'blocks'/mode/progress['checkpoint_file']).resolve()))
                history=rt.solve_in_blocks(m,s,root,mode,False)
            else:m.solve();history=rt.extract_history(m,s)
            outputs.append(history);timings.append(time.perf_counter()-started)
        finally:client.remove(m)
    np.testing.assert_allclose(outputs[0]['times'],outputs[1]['times'],atol=1e-17,rtol=1e-12)
    scale=np.array([1,1,2*np.pi*s['f0']*np.sqrt(rho*C[0,0]),2*np.pi*s['f0']*np.sqrt(rho*C[0,0])])[None,None,:,None]
    error=float(np.linalg.norm((outputs[0]['values']-outputs[1]['values'])/scale)/np.linalg.norm(outputs[0]['values']/scale))
    result=dict(kind=args.kind,weighted_complete_history_relative_error=error,tolerance=1e-5,passed=error<1e-5,
                displacement_order=args.order,extended_short_final_block=args.extension_check,solve_seconds_default_and_cached_blocks=timings)
    filename='block_restart_validation.json' if args.kind=='bulk' else 'block_restart_validation_guided.json'
    if args.order!=3:filename=f'block_restart_validation_{args.kind}_p{args.order}.json'
    if args.extension_check:filename=f'window_extension_validation_{args.kind}_p{args.order}.json'
    rt.write_json(base.parent/filename,result)
    print(json.dumps(result),flush=True)
    if not result['passed']:raise RuntimeError('Block integration changed the receiver solution.')


if __name__=='__main__':main()
