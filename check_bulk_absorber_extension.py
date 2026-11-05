"""Test the longer absorber at the previously worst frequency before a full sweep."""
import run_transient as rt
import numpy as np,json,argparse


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--step',type=int,choices=[1,2],default=1);args=parser.parse_args()
    root=rt.ROOT/'output/transient_refined'/('absorber18_probe' if args.step==1 else 'absorber27_probe');client=rt.start_session(root,4)
    cfg=json.loads((rt.ROOT/'config.json').read_text());C,rho=rt.validate_material(cfg)
    frequencies=np.array([3925781.25,5e6]);s=rt.layout('bulk',8.,5e6,C,rho,'absorber','refined_v5' if args.step==1 else 'refined_v6')
    results={}
    for case in ['reference','group']:
        folder=root/case;folder.mkdir(parents=True,exist_ok=True)
        pores,_=rt.grouped_pores(80.,8.,101,.015)
        if case=='reference':pores=[]
        if not (folder/'P_harmonic_band.npz').exists():
            minimum=min(s['mesh']/30,min((np.sqrt(4*p.area/np.pi)*1e-6/30 for p in pores),default=s['mesh']/30))
            model=rt.create_2d(client,'absorber_probe',cfg,pores,size_m=(s['length'],s['width']),pore_offset_m=(s['start'],0),
                mesh_m=s['mesh'],mesh_min_m=minimum,lateral_periodic=True,sponge_start_m=s['sponge'],repair_tolerance_m=1e-9,**rt.builder_options(s))
            try:
                rt.write_json(folder/'boundary_selections.json',model._boundary_selection_checks)
                rt.source_and_study(model,s,None,'P');rt.frequency_checkpoint(model,s,folder,'P',frequencies)
            finally:client.remove(model)
        h=dict(np.load(folder/'P_harmonic_band.npz'))
        results[case]=rt.analyze_ports([h],frequencies,'bulk',(C,rho),spectra_override=[h['values']])['T'][:,0,0]
        rt.write_json(folder/'settings.json',dict(layout=s,material=cfg,case=case,frequencies=frequencies.tolist()))
    baseline=rt.ROOT/('output/transient_converged/spatial_checks/bulk/absorber' if args.step==1 else 'output/transient_refined/absorber18_probe')
    f=np.load(baseline/'group/P_harmonic_band.npz')['frequencies'];indices=[int(np.argmin(abs(f-value))) for value in frequencies]
    np.testing.assert_allclose(f[indices],frequencies)
    baseline_results={}
    for case in ['reference','group']:
        h=dict(np.load(baseline/case/'P_harmonic_band.npz'))
        r=rt.analyze_ports([h],frequencies,'bulk',(C,rho),spectra_override=[h['values'][indices]])
        baseline_results[case]=r['T'][:,0,0]
    a=baseline_results['group']/baseline_results['reference']
    b=results['group']/results['reference'];delta=abs(abs(a)-abs(b));phase=abs(np.angle(a/b))*180/np.pi
    relative=delta/np.maximum(np.maximum(abs(1-abs(a)),abs(1-abs(b))),1e-12)
    result=dict(scope='Two-frequency diagnostic only, not full-band acceptance.',frequencies_hz=frequencies.tolist(),
        base_absorber_mm=12. if args.step==1 else 18.,test_absorber_mm=18. if args.step==1 else 27.,magnitude_delta_percentage_points=(100*delta).tolist(),
        relative_loss_change_percent=(100*relative).tolist(),phase_difference_deg=phase.tolist(),
        passed=bool(max(delta)<=2e-4 and max(relative)<=.05 and max(phase)<=1))
    rt.write_json(root/'comparison.json',result);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
