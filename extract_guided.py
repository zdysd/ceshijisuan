"""Joint modal displacement/traction projection of solved A0/S0 models."""
import mph
import site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))
import argparse,json,sys
import numpy as np
from cfrp_pores.comsol import interpolate
from cfrp_pores.dispersion import Plate,validate_material,modal_direction_basis,separate_modal_directions
from cfrp_pores.attenuation import fit_attenuation,fit_quality
from cfrp_pores.pores import write_csv
sys.stdout.reconfigure(encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder');args=p.parse_args()
    root=Path(args.folder);cfg=json.loads(Path('config.json').read_text());C,rho=validate_material(cfg)
    h=cfg['plate_thickness_m'];frequency=200000.
    plate=Plate(h,C[np.ix_([0,1,5],[0,1,5])],rho,60)
    y=np.linspace(h*1e-5,h*(1-1e-5),401)
    basis,labels,scale=modal_direction_basis(y-h/2,plate,frequency)
    client=mph.Client(cores=2,version='6.2');results=[]
    for kind in ['S0','A0']:
        ref=None
        for case in ['no_pore','random_pores']:
            model=client.load(str((root/(kind+'_'+case+'.mph')).resolve()))
            positions=np.linspace(.012,.044,5);rows=[];values=[]
            for x in positions:
                d=interpolate(model,[np.full(len(y),x),y],['u','v','solid.sx','solid.sxy'])
                coeff,residual=separate_modal_directions(y-h/2,*[v[-1] for v in d],basis,labels,scale)
                values.append(coeff[kind+'_forward'])
                rows.append(dict(x_m=x,nonmodal_residual=residual,**{name+'_abs':abs(value) for name,value in coeff.items()}))
            client.remove(model);values=np.array(values)
            if ref is None:ref=values
            fit=fit_attenuation(positions,values,ref)
            results.append(dict(mode=kind,case=case,frequency_hz=frequency,**fit,**fit_quality(fit),
                                nonmodal_residual_max=max(row['nonmodal_residual'] for row in rows)))
            write_csv(root/(kind+'_'+case+'_directions.csv'),rows)
            write_csv(root/(kind+'_'+case+'_forward_receivers.csv'),[dict(x_m=x,real=z.real,imag=z.imag) for x,z in zip(positions,values)])
    (root/'forward_results.json').write_text(json.dumps(results,indent=2))
    print(json.dumps(results,indent=2),flush=True)
