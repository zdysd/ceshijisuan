"""Refine the matched single-pore meshes and recheck normalized transmission."""
import mph
import site
from pathlib import Path
ROOT = Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import json
import sys
import numpy as np
from verify_outputs import read, forward
from cfrp_pores.comsol import save_and_solve
from cfrp_pores.dispersion import validate_material
sys.stdout.reconfigure(encoding='utf-8')


def main():
    out = ROOT/'output/revalidation'
    cfg = read(ROOT/'config.json'); C,rho = validate_material(cfg)
    impedance=np.sqrt(rho*C[0,0]); client=mph.Client(cores=2,version='6.2')
    coarse=read(out/'incident_normalized_transmission.json')['results']
    results=[]
    for tag in ['5','10']:
        f=float(tag)*1e6; root=ROOT/f'output/final_ensemble_{tag}MHz'
        mesh_size=read(root/'reference/measurement.json')['mesh_m']/2
        samples={}
        for label,path in [('reference',root/'reference/model.mph'),
                           ('median_single',out/f'median_single_{tag}MHz.mph'),
                           ('largest_single',out/f'single_pore_{tag}MHz.mph')]:
            model=client.load(str(path))
            try:
                mesh=model.java.component('comp1').mesh('mesh1')
                mesh.feature('size1').set('hmax',str(mesh_size))
                mesh.feature('size1').set('hmin',str(mesh_size/30))
                mesh.run()
                save_and_solve(model,out/f'{label}_{tag}MHz_fine.mph')
                samples[label]=forward(model,[.00015,.00145,.00155,.00165],.001,f,impedance)
            finally:
                client.remove(model)
            print('FINE TRANSMISSION',tag,label,flush=True)
        reference=samples['reference'][1:]/samples['reference'][0]
        old=next(r for r in coarse if r['frequency_hz']==f)
        for label in ['median_single','largest_single']:
            transmission=(samples[label][1:]/samples[label][0])/reference
            change=float(np.mean(1-abs(transmission)))
            previous=old[label+'_amplitude_reduction']
            results.append(dict(frequency_hz=f,case=label,coarse_amplitude_reduction=previous,
                                fine_amplitude_reduction=change,
                                absolute_reduction_change=abs(change-previous),
                                relative_reduction_change=abs(change-previous)/max(abs(change),abs(previous)),
                                phase_shift_rad=float(np.angle(transmission.mean())),
                                post_section_variation=float(np.ptp(abs(transmission)))))
    # The corrected 5 MHz domain study also supplies matched full-group meshes.
    corrected=ROOT/'output/revalidation_corrected_domain'
    if (corrected/'mesh_fine/pores/model.mph').exists():
        values={}
        for case in ['base','mesh_fine']:
            for label in ['reference','pores']:
                model=client.load(str(corrected/case/label/'model.mph'))
                try:
                    signal=forward(model,[.00015,.00145,.00155,.00165],.001,5e6,impedance)
                    values[(case,label)]=signal[1:]/signal[0]
                finally:client.remove(model)
        coarse_T=values[('base','pores')]/values[('base','reference')]
        fine_T=values[('mesh_fine','pores')]/values[('mesh_fine','reference')]
        c=float(np.mean(1-abs(coarse_T))); f=float(np.mean(1-abs(fine_T)))
        results.append(dict(frequency_hz=5e6,case='group_seed101',coarse_amplitude_reduction=c,
                            fine_amplitude_reduction=f,absolute_reduction_change=abs(f-c),
                            relative_reduction_change=abs(f-c)/max(abs(c),abs(f))))
    (out/'transmission_mesh_check.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps(results,indent=2),flush=True)


if __name__=='__main__':main()
