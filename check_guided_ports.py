"""Recheck old guide models at uniform ports before/after the porous segment."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import json,sys
import numpy as np
from cfrp_pores.comsol import interpolate
from cfrp_pores.dispersion import Plate,validate_material,modal_direction_basis,separate_modal_directions
sys.stdout.reconfigure(encoding='utf-8')


def ports(model,plate,f,positions,samples=401):
    h=plate.thickness_m;y=np.linspace(h*1e-5,h*(1-1e-5),samples)
    basis,labels,scale=modal_direction_basis(y-h/2,plate,f)
    rows=[]
    for x in positions:
        data=interpolate(model,[np.full(len(y),x),y],['u','v','solid.sx','solid.sxy'])
        coeff,residual=separate_modal_directions(y-h/2,*[v[-1] for v in data],basis,labels,scale)
        rows.append(dict(x_m=float(x),residual=residual,coefficients=coeff))
    return rows


def main():
    cfg=json.loads((ROOT/'config.json').read_text());C,rho=validate_material(cfg)
    client=mph.Client(cores=2,version='6.2');out=ROOT/'output/guided_collective';out.mkdir(exist_ok=True)
    positions=[.004,.005,.006,.007,.050,.051,.052,.053,.054]
    arrays={};results=[]
    for kind in ['A0','S0']:
        reference=None
        for case in ['no_pore','random_pores']:
            model=client.load(str(ROOT/'output/comsol_guided_corrected'/f'{kind}_{case}.mph'))
            try:
                rows=ports(model,Plate(.001,C[np.ix_([0,1,5],[0,1,5])],rho,60),200000,positions)
                refined=ports(model,Plate(.001,C[np.ix_([0,1,5],[0,1,5])],rho,120),200000,positions,801)
            finally:client.remove(model)
            signal=np.array([r['coefficients'][kind+'_forward'] for r in rows])
            fine=np.array([r['coefficients'][kind+'_forward'] for r in refined])
            wave=Plate(.001,C[np.ix_([0,1,5],[0,1,5])],rho,120).at_frequency(200000,kind)
            corrected=signal*np.exp(1j*wave['k_rad_m']*np.array(positions))
            normalized=corrected[4:]/corrected[:4].mean()
            if reference is None:reference=normalized
            T=normalized/reference
            results.append(dict(mode=kind,case=case,frequency_hz=200000,
                                pre_residual_max=max(r['residual'] for r in rows[:4]),
                                post_residual_max=max(r['residual'] for r in rows[4:]),
                                refined_post_residual_max=max(r['residual'] for r in refined[4:]),
                                post_backward_forward_ratio_max=max(abs(r['coefficients'][kind+'_backward']/r['coefficients'][kind+'_forward']) for r in rows[4:]),
                                projection_relative_amplitude_difference=float(np.max(abs(abs(fine)/abs(signal)-1))),
                                transmission_magnitude_mean=float(abs(T).mean()),
                                amplitude_reduction=float(np.mean(1-abs(T))),
                                note='Uniform-port extraction with reference-wavenumber phase correction before averaging input amplitudes.'))
            arrays[f'{kind}_{case}_coefficients']=np.array([[r['coefficients'][label] for label in ['S0_forward','S0_backward','A0_forward','A0_backward']] for r in rows])
            print(json.dumps(results[-1]),flush=True)
    np.savez_compressed(out/'old_uniform_port_signals.npz',positions_m=positions,**arrays)
    (out/'old_uniform_port_check.json').write_text(json.dumps(results,indent=2),encoding='utf-8')


if __name__=='__main__':main()
