"""Convergence of the separate through-thickness reference eigenproblem."""
import json
from pathlib import Path
import numpy as np
from cfrp_pores.dispersion import Plate,validate_material

root=Path(__file__).resolve().parent
config=json.loads((root/'config.json').read_text(encoding='utf-8'))
C,rho=validate_material(config);block=C[np.ix_([0,1,5],[0,1,5])]
coarse=Plate(.001,block,rho,120);fine=Plate(.001,block,rho,240)
rows=[]
for f in [50000.,100000.,200000.]:
    for mode in ['S0','A0']:
        a=coarse.at_frequency(f,mode);b=fine.at_frequency(f,mode)
        av=np.concatenate([a['u'],a['v']]);bv=np.concatenate([np.interp(a['y_m'],b['y_m'],b[k]) for k in ['u','v']])
        scale=np.vdot(bv,av)/np.vdot(bv,bv)
        rows.append(dict(frequency_hz=f,mode=mode,wavelength_m=a['wavelength_m'],
            relative_wavelength_difference=abs(a['wavelength_m']/b['wavelength_m']-1),
            relative_profile_difference=float(np.linalg.norm(av-scale*bv)/np.linalg.norm(av))))
result=dict(coarse_elements=120,fine_elements=240,rows=rows,
    first_higher_cutoff_hz=fine.first_higher_cutoff_hz())
(root/'output/guided_collective/basis_convergence.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
