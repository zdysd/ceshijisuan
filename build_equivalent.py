"""Verified scalar COMSOL propagation model for ONE calibrated P-wave direction.

This is a modal/coherent propagation surrogate, not a full CFRP material tensor.
It does not automatically transfer bulk coefficients to A0/S0 or other frequencies.
"""
import mph
import site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))
import argparse,json,sys
import jpype
import numpy as np
from cfrp_pores.comsol import interpolate
from cfrp_pores.dispersion import validate_material
from cfrp_pores.attenuation import fit_attenuation,fit_quality
sys.stdout.reconfigure(encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('ensemble')
    parser.add_argument('--validation',default='output/revalidation_corrected_domain/forward_gates.json')
    parser.add_argument('--output',default='output/equivalent_bulk_5MHz')
    args=parser.parse_args();root=Path(args.ensemble);out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    summary=json.loads((root/'forward_summary.json').read_text());gates=json.loads(Path(args.validation).read_text())
    if gates.get('frequency_hz')!=summary['frequency_hz']:raise SystemExit('Validation and calibration frequencies do not match.')
    if not gates['numerical_pass'] or not summary['quality']['exponential_fit_accepted'] or not summary['sample_count_gate']['passed']:
        raise SystemExit('Equivalent model blocked: domain/mesh/ensemble gates have not passed.')
    if summary['fit']['alpha_np_m']<=0 or summary['bootstrap_ci_95_np_m'][0]<=0:
        raise SystemExit('Equivalent positive attenuation is not statistically identified.')
    cfg=json.loads(Path('config.json').read_text());C,rho=validate_material(cfg)
    f=summary['frequency_hz'];k0=2*np.pi*f/np.sqrt(C[0,0]/rho)
    kr=k0-summary['fit']['phase_slope_rad_m'];alpha=summary['fit']['alpha_np_m'];length=.001
    parameters=dict(frequency_hz=f,k_real_rad_m=kr,alpha_np_m=alpha,alpha_ci_95_np_m=summary['bootstrap_ci_95_np_m'],
                    phase_velocity_m_s=2*np.pi*f/kr,phase_convention='exp(+i*w*t), spatial exp(-(alpha+i*k_real)*x)',
                    material_status=cfg['material_status'],mode='P, global x',full_material_tensor=False,
                    demonstration_length_m=length,extrapolation_to_large_specimens_validated=False,
                    validation_source=str(Path(args.validation)),
                    note='Numerically calibrated propagation surrogate; dimensional/ensemble errors and material assumptions remain. Do not use as a general laminate constitutive law.')
    (out/'parameters.json').write_text(json.dumps(parameters,indent=2))
    client=mph.Client(cores=2,version='6.2');m=client.create('coherent_P_equivalent');j=m.java
    j.param().set('L',f'{length}[m]');j.param().set('kr',f'{kr}[1/m]');j.param().set('alpha',f'{alpha}[1/m]');j.param().set('f0',f'{f}[Hz]')
    c=j.component().create('comp1',True);g=c.geom().create('geom1',1)
    interval=g.create('i1','Interval');interval.set('p1','0');interval.set('p2','L');g.run()
    phys=c.physics().create('c','CoefficientFormPDE','geom1')
    equation=phys.feature('cfeq1');equation.set('c','1');equation.set('a','-(kr-i*alpha)^2');equation.set('f','0');equation.set('ea','0');equation.set('da','0')
    left=phys.create('left','DirichletBoundary',0);left.selection().set(jpype.JArray(jpype.JInt)([1]));left.set('r','1')
    right=phys.create('right','FluxBoundary',0);right.selection().set(jpype.JArray(jpype.JInt)([2]));right.set('q','i*(kr-i*alpha)');right.set('g','0')
    mesh=c.mesh().create('mesh1');sz=mesh.create('size1','Size');sz.set('custom',True);sz.set('hmax',str(2*np.pi/kr/16));mesh.create('edge','Edge');mesh.run()
    j.study().create('std1');j.study('std1').create('freq','Frequency');j.study('std1').feature('freq').set('plist','f0')
    m.solve();points=np.linspace(0,length,121);solution=interpolate(m,[points],['u'])[0,-1]
    exact=np.exp(-(alpha+1j*kr)*points);fit=fit_attenuation(points,solution)
    error=float(np.max(abs(solution-exact)))
    result=dict(max_complex_amplitude_error=error,analytic_comparison_passed=error<.01,
                fitted_alpha_np_m=fit['alpha_np_m'],relative_alpha_error=abs(fit['alpha_np_m']/alpha-1),
                dofs=len(j.sol('sol1').getU()),mesh_elements=j.component('comp1').mesh('mesh1').getNumElem())
    if not result['analytic_comparison_passed'] or result['relative_alpha_error']>.05:raise RuntimeError('Scalar surrogate failed its analytic verification.')
    pg=j.result().create('wave','PlotGroup1D');pg.create('line','LineGraph');pg.feature('line').set('expr','abs(u)')
    m.save(str((out/'equivalent_P_5MHz.mph').resolve()))
    np.savez_compressed(out/'comparison.npz',x_m=points,comsol=solution,analytic=exact)
    (out/'verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
