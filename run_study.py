"""COMSOL screening, coherent ensembles and numerical convergence checks.

Run with the installed COMSOL Python runtime, like run_comsol.py. Screening is
guarded by numerical-quality checks; outputs are not declared validated merely
because a solver completed. Completed receiver files can be resumed.
"""
import mph
import site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))
import argparse
import json
import hashlib
import sys
import numpy as np
from cfrp_pores.comsol import start_client,create_2d,add_load,add_frequency_study,save_and_solve,receiver_amplitudes
from cfrp_pores.pores import PoreSpec,export_sample,load_polygons,simplify_geometry,write_csv
from cfrp_pores.attenuation import fit_attenuation,fit_quality,uncertainty,coherent_mean,converged_alpha
from cfrp_pores.dispersion import validate_material
sys.stdout.reconfigure(encoding='utf-8')
_ROOT=Path(__file__).parent
IMPLEMENTATION_HASH=hashlib.sha256(b''.join(p.read_bytes() for p in [_ROOT/'run_study.py',*sorted((_ROOT/'cfrp_pores').glob('*.py'))])).hexdigest()


def run_one(client,root,cfg,porosity,seed,frequency=5e6,mesh=4e-5,length=.001,width=.001,tolerance=.1):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    C,rho=validate_material(cfg)
    mesh=min(mesh,np.sqrt(C[5,5]/rho)/(6*frequency))
    result_path=root/'measurement.json'
    if result_path.exists():
        data=json.loads(result_path.read_text())
        expected=dict(porosity=porosity,seed=seed,frequency_hz=frequency,mesh_m=mesh,length_m=length,width_m=width,requested_contour_tolerance_um=tolerance,
                      implementation_hash=IMPLEMENTATION_HASH,material_hash=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest())
        if all(data.get(k)==v for k,v in expected.items()):return data
        raise ValueError(f'Resume configuration mismatch in {result_path}; use a new output directory.')
    # Matched solid edges permit a transverse periodic plane-wave cell. Small
    # pore-free edge buffers are reported and eliminated in the width study.
    edge_buffer=40e-6
    interior_width=width-2*edge_buffer
    interior_porosity=porosity*width/interior_width
    spec=PoreSpec((length*1e6,interior_width*1e6),cfg['pixel_um'],tuple(cfg['correlation_lengths_um']),interior_porosity,seed,cfg['orientation_deg'])
    if porosity:
        export_sample(root/'geometry',spec)
        polygons,_=load_polygons(root/'geometry/geometry.json')
        original=polygons
        trial_tolerance=tolerance
        while True:
            polygons,simplification=simplify_geometry(original,trial_tolerance)
            actual=sum(p.area for p in polygons)/(length*width*1e12)
            if abs(actual-porosity)<=cfg['geometry_porosity_absolute_tolerance'] and abs(simplification['relative_void_area_change'])<=.02:break
            if trial_tolerance==0:raise RuntimeError('Actual geometric porosity misses tolerance.')
            trial_tolerance=trial_tolerance/2 if trial_tolerance>.00625 else 0.
    else:polygons=[];simplification=None
    receivers=np.linspace(.0003+.15*length,.0003+.95*length,9)
    model=None
    try:
        print(f'RUN P={porosity*100:g}% seed={seed} mesh={mesh:g} L={length:g} W={width:g}',flush=True)
        model=create_2d(client,'screen',cfg,polygons,size_m=(length+.001,width),pore_offset_m=(.0003,edge_buffer),mesh_m=mesh,lateral_periodic=True)
        add_load(model,['1e6[Pa]','0','0']);add_frequency_study(model,frequency)
        data=save_and_solve(model,root/'model.mph')
        values,missing=receiver_amplitudes(model,receivers,width,801,bulk_impedance=np.sqrt(rho*C[0,0]))
        data.update(porosity=porosity,seed=seed,frequency_hz=frequency,mesh_m=mesh,length_m=length,width_m=width,
                    requested_contour_tolerance_um=tolerance,contour_tolerance_um=simplification['tolerance_um'] if simplification else tolerance,simplification=simplification,
                    positions_m=receivers.tolist(),amplitude_real=values.real.tolist(),amplitude_imag=values.imag.tolist(),
                    missing_receiver_fraction=missing,receiver_definition='reference-impedance forward displacement/stress projection',status='solved',
                    lateral_boundary='intentional_periodic_cell_with_width_convergence',edge_buffer_m=edge_buffer,
                    pore_roi_porosity=interior_porosity,
                    actual_vector_porosity=sum(p.area for p in polygons)/(length*width*1e12),implementation_hash=IMPLEMENTATION_HASH,
                    material_hash=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest())
        result_path.write_text(json.dumps(data,indent=2),encoding='utf-8')
        return data
    finally:
        if model is not None:client.remove(model)


def amplitudes(record):return np.array(record['amplitude_real'])+1j*np.array(record['amplitude_imag'])


def validation(client,root,cfg):
    results=[];root=Path(root)
    # Same physical random field and mask in the three discretization comparisons.
    for label,mesh,tol,length,width in [
        ('base',4e-5,.1,.001,.001),('mesh_fine',2e-5,.1,.001,.001),
        ('contour_full',4e-5,0.,.001,.001),('length_long',4e-5,.1,.002,.001),
        ('width_wide',4e-5,.1,.001,.002)]:
        reference=run_one(client,root/label/'reference',cfg,0,101,mesh=mesh,length=length,width=width,tolerance=tol)
        pore=run_one(client,root/label/'pores',cfg,.015,101,mesh=mesh,length=length,width=width,tolerance=tol)
        fitted=fit_attenuation(pore['positions_m'],amplitudes(pore),amplitudes(reference))
        row=dict(case=label,**fitted,**fit_quality(fitted),solve_seconds=pore['solve_seconds'],mesh_elements=pore['mesh_elements'],dofs=pore['dofs'],
                 actual_porosity=pore['actual_vector_porosity'],requested_porosity=pore['porosity'])
        results.append(row);write_csv(root/'validation.csv',results)
    base=results[0];gates={}
    for row in results[1:]:
        gates[row['case']]=converged_alpha(base['alpha_np_m'],row['alpha_np_m'],base['distance_span_m'])
    gates['single_sample_exponential_fit']=dict(passed=base['exponential_fit_accepted'])
    gates['geometry_porosity']=dict(passed=all(abs(r['actual_porosity']-r['requested_porosity'])<=cfg['geometry_porosity_absolute_tolerance'] for r in results))
    gates['material_is_measured']=dict(passed=cfg['material_status']=='measured_CFRP')
    # Material accuracy is distinct from numerical implementation accuracy.
    numerical_pass=all(gates[k]['passed'] for k in ['mesh_fine','contour_full','length_long','width_wide','single_sample_exponential_fit','geometry_porosity'])
    (root/'gates.json').write_text(json.dumps(dict(gates=gates,numerical_pass=numerical_pass,
       relative_domain_coefficient_passed=all(gates[k]['relative_coefficient_passed'] for k in ['length_long','width_wide']),
       explanation='Domain changes also change finite random samples. If they fail, ensemble convergence is required before inferring material attenuation.'),indent=2))
    return numerical_pass


def ensemble(client,root,cfg,porosity=.015,samples=10,frequency=5e6):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    ref=run_one(client,root/'reference',cfg,0,101,frequency)
    individual=[];signals=[];running=[]
    for seed in range(101,101+samples):
        record=run_one(client,root/f'seed_{seed}',cfg,porosity,seed,frequency)
        z=amplitudes(record);signals.append(z)
        fit=fit_attenuation(record['positions_m'],z,amplitudes(ref))
        individual.append(dict(seed=seed,porosity=porosity,**fit,**fit_quality(fit)))
        mean_fit=fit_attenuation(record['positions_m'],coherent_mean(signals),amplitudes(ref))
        running.append(dict(samples=len(signals),**mean_fit,**fit_quality(mean_fit)))
        write_csv(root/'individual.csv',individual);write_csv(root/'coherent_running.csv',running)
    summary=dict(porosity=porosity,frequency_hz=frequency,material_status=cfg['material_status'],
                 individual_statistics=uncertainty([r['alpha_np_m'] for r in individual]),
                 coherent_attenuation=running[-1],
                 warning='CI of individual fitted slopes is not the CI of the coherent attenuation estimator.')
    # Bootstrap entire realizations, preserving every receiver's phase.
    rng=np.random.default_rng(1234);bootstrap=[]
    for _ in range(500):
        chosen=rng.integers(samples,size=samples)
        bootstrap.append(fit_attenuation(record['positions_m'],np.array(signals)[chosen].mean(axis=0),amplitudes(ref))['alpha_np_m'])
    summary['coherent_bootstrap_95_percent_ci_np_m']=np.percentile(bootstrap,[2.5,97.5]).tolist()
    if samples>=10:
        summary['sample_count_gate']=converged_alpha(running[samples//2-1]['alpha_np_m'],running[-1]['alpha_np_m'],running[-1]['distance_span_m'])
    (root/'summary.json').write_text(json.dumps(summary,indent=2))
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['validate','ensemble','screen'])
    p.add_argument('--output',default='output/study');p.add_argument('--samples',type=int,default=10)
    p.add_argument('--porosity',type=float,default=.015);p.add_argument('--frequency',type=float,default=5e6)
    p.add_argument('--config',default='config.json')
    args=p.parse_args();cfg=json.loads(Path(args.config).read_text());client=start_client()
    if args.command=='validate':validation(client,args.output,cfg)
    elif args.command=='ensemble':
        if args.samples<2: p.error('An ensemble needs at least two independent samples.')
        ensemble(client,args.output,cfg,args.porosity,args.samples,args.frequency)
    else:
        gates_path=Path(args.output)/'validation/gates.json'
        passed=validation(client,Path(args.output)/'validation',cfg)
        if not passed:
            raise SystemExit('Numerical gates failed. Use ensemble/domain validation before launching a full porosity sweep.')
        for phi in cfg['porosity_fractions']:
            if phi:ensemble(client,Path(args.output)/f'phi_{phi:g}',cfg,phi,args.samples,args.frequency)
