"""Recheck solved mesh/contour/domain cases with forward-wave projection."""
import mph
import site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))
import argparse,json,sys
from extract_forward import extract
from cfrp_pores.attenuation import fit_attenuation,fit_quality,converged_alpha
from cfrp_pores.pores import write_csv
sys.stdout.reconfigure(encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder');args=p.parse_args()
    root=Path(args.folder);cfg=json.loads(Path('config.json').read_text());client=mph.Client(cores=2,version='6.2')
    rows=[]
    for case in ['base','mesh_fine','contour_full','length_long','width_wide']:
        reference=json.loads((root/case/'reference/measurement.json').read_text())
        pores=json.loads((root/case/'pores/measurement.json').read_text())
        fref,_=extract(client,reference['model'],cfg,reference);forward,_=extract(client,pores['model'],cfg,pores)
        fit=fit_attenuation(pores['positions_m'],forward,fref)
        rows.append(dict(case=case,**fit,**fit_quality(fit),dofs=pores['dofs'],mesh_elements=pores['mesh_elements'],
                         actual_porosity=pores['actual_vector_porosity'],requested_porosity=pores['porosity']))
        print(case,fit['alpha_np_m'],flush=True)
    write_csv(root/'forward_validation.csv',rows)
    base=rows[0]
    gates={r['case']:converged_alpha(base['alpha_np_m'],r['alpha_np_m'],base['distance_span_m']) for r in rows[1:]}
    gates['exponential_fit']=dict(passed=base['exponential_fit_accepted'])
    gates['geometry_porosity']=dict(passed=all(abs(r['actual_porosity']-r['requested_porosity'])<=cfg['geometry_porosity_absolute_tolerance'] for r in rows))
    (root/'forward_gates.json').write_text(json.dumps(dict(gates=gates,frequency_hz=reference['frequency_hz'],alpha_range_np_m=[min(r['alpha_np_m'] for r in rows),max(r['alpha_np_m'] for r in rows)],numerical_pass=all(g['passed'] for g in gates.values()),
        relative_domain_coefficient_passed=all(gates[k]['relative_coefficient_passed'] for k in ['length_long','width_wide']),
        warning='Near-zero absolute error gates do not establish 5% relative accuracy in a small nonzero coefficient.'),indent=2))
