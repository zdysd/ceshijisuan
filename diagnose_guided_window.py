"""Compare receiver end tapers with the independently solved harmonic response.

This diagnostic never changes validation or production data.
"""
import json
from pathlib import Path
import run_transient as rt
import numpy as np


def main():
    root=rt.ROOT/'output/transient_converged'
    base=root/'guided/base/L8/seed_101'
    record=json.loads((base/'group/run.json').read_text(encoding='utf-8'))
    s=record['signature']['layout'];C,rho=rt.validate_material(record['signature']['material'])
    histories={c:[dict(np.load(base/c/f'{m}_history.npz')) for m in ['S0','A0']] for c in ['reference','group']}
    frequencies=rt.analysis_frequencies(s)
    plate=rt.Plate(s['width'],C[np.ix_([0,1,5],[0,1,5])],rho,s['safe_elements'])
    bank=rt.modal_bank(plate,frequencies,histories['group'][0]['y'])
    harmonic={c:np.load(base/c/'spectra.npz')['harmonic_band_T'] for c in histories}
    target=np.diagonal(harmonic['group'],axis1=1,axis2=2)/np.diagonal(harmonic['reference'],axis1=1,axis2=2)
    rows=[]
    for end in [100.,120.,140.]:
        for taper in [0.,10.,20.,40.]:
            outputs={}
            for case,launches in histories.items():
                subset=[]
                for h in launches:
                    keep=h['times']<=end*1e-6
                    t=h['times'][keep];v=h['values'][keep].copy()
                    if taper:
                        w=.5*(1+np.cos(np.pi*np.clip((t-(t[-1]-taper*1e-6))/(taper*1e-6),0,1)))
                        v*=w[:,None,None,None]
                    subset.append(dict(h,times=t,values=v))
                outputs[case]=rt.analyze_ports(subset,frequencies,'guided',(C,rho),bank)['T']
            actual=np.diagonal(outputs['group'],axis1=1,axis2=2)/np.diagonal(outputs['reference'],axis1=1,axis2=2)
            delta=abs(abs(actual)-abs(target))
            relative=delta/np.maximum(np.maximum(abs(1-abs(actual)),abs(1-abs(target))),1e-12)
            row=dict(end_us=end,taper_us=taper,magnitude_delta_pp=float(delta.max()*100),relative_loss_percent=float(relative.max()*100),
                     center_loss_percent=((1-abs(actual[8]))*100).tolist(),normalized_T_real=actual.real.tolist(),normalized_T_imag=actual.imag.tolist())
            rows.append(row)
            print({k:v for k,v in row.items() if not k.startswith('normalized')},flush=True)
    rt.write_json(root/'guided_window_diagnostic.json',dict(frequencies=frequencies.tolist(),rows=rows))


if __name__=='__main__':main()
