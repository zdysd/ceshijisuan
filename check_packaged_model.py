"""Solve a delivered copy in memory only; its file must remain solution-free."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import json
import numpy as np
from guided_collective import settings,modal_port_data,OUT
from package_comsol_models import DEST

client=mph.Client(cores=2,version='6.2');s=settings()
path=DEST/'main_models/checks/d20_maximum_fine/model.mph'
raw=np.load(OUT/'checks/d20_maximum_fine/signals.npz');rows=[]
model=client.load(str(path))
try:
    assert model.java.sol('sol1').isEmpty()
    for fi,mode_id,f,mode in [(3,1,200000.,'A0'),(1,0,50000.,'S0')]:
        model.java.param().set('freq_id',str(fi));model.java.param().set('mode_id',str(mode_id))
        model.solve();data=modal_port_data(model,s,f)
        original=raw[f'{int(f)}_{mode}_coefficients']
        error=float(np.max(abs(data['coefficients']-original))/np.max(abs(original)))
        rows.append(dict(frequency_hz=f,mode=mode,relative_port_coefficient_difference=error))
        assert error<1e-5,rows[-1]
finally:client.remove(model)
model=client.load(str(path))
try:assert model.java.sol('sol1').isEmpty()
finally:client.remove(model)
(DEST/'recompute_smoke_check.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
print(json.dumps(rows,indent=2),flush=True)
