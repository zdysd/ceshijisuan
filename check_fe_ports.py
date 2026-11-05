"""Check element-integrated A0 port amplitudes on the final saved 200 kHz fields."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import argparse,json,sys
import numpy as np
from guided_collective import settings,OUT,dephase_average
from guided_length_study import length_settings
from guided_projection_fe import prepare_projection,project_port
from cfrp_pores.comsol import integrate_section
sys.stdout.reconfigure(encoding='utf-8')


def main():
    p=argparse.ArgumentParser();p.add_argument('--all-completed',action='store_true');args=p.parse_args()
    s=settings();client=mph.Client(cores=1,version='6.2');rows=[];ref={};point_ref={}
    folders=[OUT/'reference',*sorted(OUT.glob('reference_L*'))]
    folders+=sorted(OUT.glob('d*/single_*/model.mph'))
    folders+=sorted(OUT.glob('d*/L*/seed_*/model.mph'))
    if not args.all_completed:
        folders=[OUT/'reference',OUT/'d20/L8/seed_101',OUT/'d80/L8/seed_101',OUT/'d80/single_maximum']
    for item in folders:
        folder=item.parent if item.suffix=='.mph' else item
        path=folder/'model.mph'
        if not path.exists():continue
        length=8
        if folder.name.startswith('reference_L'):length=int(folder.name.replace('reference_L',''))
        elif folder.parent.name.startswith('L'):length=int(folder.parent.name[1:])
        local=length_settings(s,length)
        model=client.load(str(path))
        try:
            f=float(model.java.param().evaluate('f0'))
            if f!=200000.:raise ValueError('This check expects saved A0/200 kHz state.')
            proj=prepare_projection(model,s['plate'],f)
            ports=[project_port(model,x,proj) for x in local['pre_positions']+local['post_positions']]
            net_power=[float(integrate_section(model,x,.001,
                ['-0.5*real(solid.sx*conj(i*2*pi*f0*u)+solid.sxy*conj(i*2*pi*f0*v))'],8)[0].real)
                for x in local['pre_positions']+local['post_positions']]
        finally:client.remove(model)
        k=s['plate'].at_frequency(f,'A0')['k_rad_m']
        a=np.array([r['A0_forward'] for r in ports]);b=np.array([r['A0_backward'] for r in ports])
        coefficient=dephase_average(a[3:],local['post_positions'],k)/dephase_average(a[:3],local['pre_positions'],k)
        if folder.name.startswith('reference'):ref[length]=coefficient
        reference=ref.get(length,ref[8])
        T=coefficient/reference
        raw=np.load(folder/'signals.npz');point=raw['200000_A0_coefficients'][:,2]
        point_coefficient=dephase_average(point[3:],local['post_positions'],k)/dephase_average(point[:3],local['pre_positions'],k)
        if folder.name.startswith('reference'):point_ref[length]=point_coefficient
        point_T=point_coefficient/point_ref.get(length,point_ref[8])
        row=dict(case=str(folder.relative_to(OUT)),length_mm=length,frequency_hz=f,mode='A0',
                 T_real=float(T.real),T_imag=float(T.imag),amplitude_reduction=float(1-abs(T)),
                 right_back_forward_max=float(np.max(abs(b[3:]/a[3:]))),
                 least_squares_direct_reduction=float(1-abs(point_T)),
                 method_magnitude_difference=float(abs(abs(T)-abs(point_T))),
                 net_power_out_over_in=float(np.mean(net_power[3:])/np.mean(net_power[:3])),
                 extraction='Betti reciprocity with COMSOL finite-element cut-line integration')
        rows.append(row);print(json.dumps(row),flush=True)
    (OUT/'fe_port_check.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')


if __name__=='__main__':main()
