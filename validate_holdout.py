"""Unseen random geometries check the calibrated forward propagation surrogate."""
import mph
import site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))
import json,sys
import numpy as np
from run_study import run_one,amplitudes
sys.stdout.reconfigure(encoding='utf-8')

if __name__=='__main__':
    cfg=json.loads(Path('config.json').read_text())
    parameters=json.loads(Path('output/equivalent_bulk_5MHz/parameters.json').read_text())
    root=Path('output/held_out');root.mkdir(parents=True,exist_ok=True)
    client=mph.Client(cores=2,version='6.2')
    reference=run_one(client,root/'reference',cfg,0,211)
    rows=[];signals=[]
    kr=parameters['k_real_rad_m'];alpha=parameters['alpha_np_m']
    C=np.array(cfg['stiffness_GPa_voigt_xx_yy_zz_yz_xz_xy'])*1e9
    k0=2*np.pi*parameters['frequency_hz']/np.sqrt(C[0,0]/cfg['density_kg_m3'])
    x=np.array(reference['positions_m']);prediction=np.exp(-(alpha+1j*(kr-k0))*(x-x[0]))
    for seed in [211,212]:
        r=run_one(client,root/f'seed_{seed}',cfg,.015,seed)
        ratio=amplitudes(r)/amplitudes(reference);ratio/=ratio[0]
        signals.append(ratio)
        rows.append(dict(seed=seed,max_relative_magnitude_error=float(np.max(abs(abs(ratio)/abs(prediction)-1))),
                         max_complex_amplitude_error=float(np.max(abs(ratio-prediction)))))
    mean=np.array(signals).mean(axis=0)
    result=dict(unseen_seeds=[211,212],cases=rows,
                ensemble_max_magnitude_error=float(np.max(abs(abs(mean)/abs(prediction)-1))),
                ensemble_max_complex_amplitude_error=float(np.max(abs(mean-prediction))),
                tested_distance_span_m=float(np.ptp(x)),
                warning='Only same local dimensions and pore statistics; not validation of centimetre-scale extrapolation or other morphologies.')
    result['local_amplitude_prediction_passed']=result['ensemble_max_magnitude_error']<=.01
    (root/'verification.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(root/'comparison.npz',x_m=x,prediction=prediction,observed=mean)
    print(json.dumps(result,indent=2),flush=True)
