"""Read-only, single-launch cadence check on completed sampling blocks.

This is an output-quadrature diagnostic, not a replacement for dual-mode gates.
"""
import json
import numpy as np
import run_transient as rt
from cfrp_pores.transient import fourier_history
from cfrp_pores.dispersion import separate_modal_directions


def main():
    folder=rt.ROOT/'output/transient_refined/guided_sampling/guided/sampling/L8/seed_101/reference'
    record=json.loads((folder/'run.json').read_text());s=record['signature']['layout']
    C,rho=rt.validate_material(record['signature']['material'])
    progress=json.loads((folder/'blocks/S0/progress.json').read_text())
    pieces=[dict(np.load(folder/'blocks/S0'/f'history_{i:03d}.npz')) for i in range(progress['completed'])]
    h=dict(times=np.concatenate([p['times'][int(i>0):] for i,p in enumerate(pieces)]),
           values=np.concatenate([p['values'][int(i>0):] for i,p in enumerate(pieces)]),y=pieces[0]['y'],positions=pieces[0]['positions'])
    # Identical endpoints make the cadence comparison independent of record length.
    count=((len(h['times'])-1)//3)*3+1
    h.update(times=h['times'][:count],values=h['values'][:count])
    f=np.array([s['fmin'],s['f0'],s['fmax']])
    bank=rt.modal_bank(rt.Plate(s['width'],C[np.ix_([0,1,5],[0,1,5])],rho,s['safe_elements']),f,h['y'])
    signals=[]
    for stride in [1,3]:
        spectrum=fourier_history(h['times'][::stride],h['values'][::stride],f)
        values=[]
        for fi,(basis,labels,scale,modes,power) in enumerate(bank):
            forward=[]
            for port in spectrum[fi]:
                coefficients,_=separate_modal_directions(h['y']-s['width']/2,*port,basis,labels,scale)
                forward.append(coefficients[labels[0]])
            forward=np.array(forward)*np.exp(1j*np.array(h['positions'])*modes[0]['k_rad_m'])
            values.append(forward[3:].mean()/forward[:3].mean())
        signals.append(np.array(values))
    result=dict(record_end_us=float(h['times'][-1]*1e6),mode='S0',case='reference',frequencies_hz=f.tolist(),
                magnitude_difference_percentage_points=(100*abs(abs(signals[0])-abs(signals[1]))).tolist(),
                complex_difference=abs(signals[0]-signals[1]).tolist(),
                note='Same completed FE trajectory, identical endpoints; single-launch projection only; not acceptance.')
    rt.write_json(folder/'sampling_checkpoint_diagnostic.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
