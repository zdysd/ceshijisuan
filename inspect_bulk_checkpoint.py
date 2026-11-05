"""Read only committed P histories. Partial spectra never count as acceptance."""
import json
from datetime import datetime,timezone,timedelta
from pathlib import Path
import numpy as np
from cfrp_pores.transient import analyze_ports,tail_fraction
from cfrp_pores.dispersion import validate_material


def main():
    root=Path(__file__).resolve().parent/'output/transient_converged'
    folder=root/'bulk/base/L8/seed_101/group';record=json.loads((folder/'run.json').read_text(encoding='utf-8-sig'))
    s=record['signature']['layout'];C,rho=validate_material(record['signature']['material'])
    state=json.loads((folder/'blocks/P/progress.json').read_text(encoding='utf-8-sig'))
    pieces=[dict(np.load(folder/'blocks/P'/f'history_{i:03d}.npz')) for i in range(state['completed'])]
    h=dict(times=np.concatenate([p['times'][int(i>0):] for i,p in enumerate(pieces)]),
           values=np.concatenate([p['values'][int(i>0):] for i,p in enumerate(pieces)]),y=pieces[0]['y'],positions=pieces[0]['positions'])
    reference_path=root/'bulk/base/L8/seed_101/reference/P_history.npz'
    if not reference_path.exists():reference_path=reference_path.with_name('P_history_12us.npz')
    reference=dict(np.load(reference_path))
    f=np.linspace(s['fmin'],s['fmax'],17)
    a=analyze_ports([h],f,'bulk',(C,rho));r=analyze_ports([reference],f,'bulk',(C,rho))
    trimmed=analyze_ports([h],f,'bulk',(C,rho),end_fraction=.9)
    group_fd=np.load(root/'spatial_checks/bulk/base/group/spectra_band.npz')['T'][:,0,0]
    reference_fd=np.load(root/'spatial_checks/bulk/base/reference/spectra_band.npz')['T'][:,0,0]
    actual=a['T'][:,0,0]/r['T'][:,0,0];target=group_fd/reference_fd
    result=dict(updated_at=datetime.now(timezone(timedelta(hours=8))).isoformat(),provisional=True,
        completed_blocks=state['completed'],end_time_us=float(h['times'][-1]*1e6),tail_signal_fraction=tail_fraction(h),
        reference_history_source=str(reference_path),
        window_complex_T_delta=float(abs(a['T']-trimmed['T']).max()),
        maximum_magnitude_delta_percentage_points=float(100*abs(abs(actual)-abs(target)).max()),
        frequencies_hz=f.tolist(),transient_amplitude_loss_percent=(100*(1-abs(actual))).tolist(),
        harmonic_amplitude_loss_percent=(100*(1-abs(target))).tolist(),
        note='Incomplete observation window, not acceptance and not an attenuation coefficient. Read-only committed checkpoints; running solver untouched.')
    output=root.parent/'transient_verified/P_checkpoint_diagnostic.json'
    temporary=output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2),encoding='utf-8');temporary.replace(output)
    print(json.dumps({k:v for k,v in result.items() if not isinstance(v,list)},indent=2))


if __name__=='__main__':main()
