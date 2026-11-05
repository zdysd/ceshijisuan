"""Re-evaluate solved frequency models using displacement + stress, no re-solve."""
import mph
import site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))
import argparse
import json
import sys
import numpy as np
from cfrp_pores.comsol import interpolate,integrate_section
from cfrp_pores.attenuation import fit_attenuation,fit_quality,forward_bulk_displacement,coherent_mean,converged_alpha
from cfrp_pores.pores import write_csv
from cfrp_pores.dispersion import validate_material
sys.stdout.reconfigure(encoding='utf-8')


def extract(client,path,config,record):
    C,rho=validate_material(config);Z=np.sqrt(C[0,0]*rho)
    model=client.load(str(Path(path).resolve()))
    y=np.linspace(record['width_m']*1e-5,record['width_m']*(1-1e-5),801)
    frequency=record['frequency_hz'];values=[];backward=[]
    for x in record['positions_m']:
        data=integrate_section(model,x,record['width_m'],['u','solid.sx','1'])
        average_u=data[0]/data[2].real
        average_stress=data[1]/record['width_m']
        forward=forward_bulk_displacement(average_u,average_stress,frequency,Z)
        values.append(forward);backward.append(average_u-forward)
    client.remove(model)
    return np.array(values),np.array(backward)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('folder');parser.add_argument('--config',default='config.json')
    args=parser.parse_args();root=Path(args.folder);config=json.loads(Path(args.config).read_text())
    client=mph.Client(cores=2,version='6.2')
    reference=json.loads((root/'reference/measurement.json').read_text())
    forward_ref,back_ref=extract(client,reference['model'],config,reference)
    rows=[];signals=[]
    for folder in sorted(root.glob('seed_*')):
        path=folder/'measurement.json'
        if not path.exists():continue
        record=json.loads(path.read_text());print('PROJECT',folder.name,flush=True)
        forward,back=extract(client,record['model'],config,record)
        signals.append(forward)
        fit=fit_attenuation(record['positions_m'],forward,forward_ref)
        rows.append(dict(seed=record['seed'],**fit,**fit_quality(fit),
                         backward_to_forward_amplitude_max=float(np.max(abs(back)/abs(forward)))))
        write_csv(root/'forward_individual.csv',rows)
    if not signals:raise RuntimeError('No completed pore models found.')
    averaged=coherent_mean(signals);fit=fit_attenuation(reference['positions_m'],averaged,forward_ref)
    rng=np.random.default_rng(1234);boot=[]
    for _ in range(500):
        ids=rng.integers(len(signals),size=len(signals))
        boot.append(fit_attenuation(reference['positions_m'],np.array(signals)[ids].mean(axis=0),forward_ref)['alpha_np_m'])
    summary=dict(samples=len(signals),frequency_hz=reference['frequency_hz'],
                 receiver_definition='intrinsic displacement + extrinsic normal stress; reference-impedance forward projection',
                 fit=fit,quality=fit_quality(fit),bootstrap_ci_95_np_m=np.percentile(boot,[2.5,97.5]).tolist(),
                 reference_backward_to_forward_max=float(np.max(abs(back_ref)/abs(forward_ref))))
    if len(signals)>=10:
        first=fit_attenuation(reference['positions_m'],np.array(signals[:len(signals)//2]).mean(axis=0),forward_ref)
        summary['sample_count_gate']=converged_alpha(first['alpha_np_m'],fit['alpha_np_m'],fit['distance_span_m'])
    (root/'forward_summary.json').write_text(json.dumps(summary,indent=2))
    write_csv(root/'forward_receivers.csv',[dict(x_m=x,real=a.real,imag=a.imag,reference_real=r.real,reference_imag=r.imag) for x,a,r in zip(reference['positions_m'],averaged,forward_ref)])
    print(json.dumps(summary,indent=2),flush=True)
