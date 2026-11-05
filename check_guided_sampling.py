"""Independent output-cadence check with unchanged internal integration steps."""
import argparse,json
from pathlib import Path
import run_transient as rt
from run_refined_pipeline import case_args
import numpy as np


def cadence_controls(rows):
    fine=np.array(rows[0]['normalized_T_real'])+1j*np.array(rows[0]['normalized_T_imag'])
    coarse=np.array(rows[1]['normalized_T_real'])+1j*np.array(rows[1]['normalized_T_imag'])
    checks=[]
    for i,mode in enumerate(['S0','A0']):
        a,b=fine[:,i],coarse[:,i];delta=abs(abs(a)-abs(b))
        magnitude=float(delta.max());relative=float(np.max(delta/np.maximum(np.maximum(abs(1-abs(a)),abs(1-abs(b))),1e-12)))
        phase=float(np.max(abs(np.angle(a/b)))*180/np.pi)
        checks.append(dict(mode=mode,magnitude_delta=magnitude,relative_loss_change=relative,phase_error_deg=phase,
                           passed=magnitude<=2e-4 and relative<=.05 and phase<=1))
    return checks


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--cores',type=int,default=2)
    parser.add_argument('--case',choices=['reference','group','both'],default='both')
    parser.add_argument('--analyze-only',action='store_true')
    args=parser.parse_args()
    root=rt.ROOT/'output/transient_refined/guided_sampling'
    config=case_args(args.cores,'refined_v6')
    rt.configure_scratch(root)
    if not args.analyze_only:
        client=rt.start_session(root,args.cores)
        for case in (['reference','group'] if args.case=='both' else [args.case]):
            rt.run_case(client,root,'guided',case,8.,101,'sampling',config)
    base=root/'guided/sampling/L8/seed_101'
    if not all((base/c/'run.json').exists() and json.loads((base/c/'run.json').read_text())['status']=='complete' for c in ['reference','group']):return
    records={c:rt.process_case(base/c,{}) for c in ['reference','group']}
    s=json.loads((base/'group/run.json').read_text())['signature']['layout']
    C,rho=rt.validate_material(json.loads((rt.ROOT/'config.json').read_text()))
    frequencies=rt.analysis_frequencies(s)
    bank=rt.modal_bank(rt.Plate(s['width'],C[np.ix_([0,1,5],[0,1,5])],rho,s['safe_elements']),frequencies,records['group'][2][0]['y'])
    results={};rows=[]
    harmonic={c:records[c][1]['harmonic_band_T'] for c in records}
    target=np.diagonal(harmonic['group'],axis1=1,axis2=2)/np.diagonal(harmonic['reference'],axis1=1,axis2=2)
    # All strides use the same endpoint (LCM of 1, 2, 3). Otherwise changing
    # cadence also changes the observation window by one or more samples.
    count=(min(len(h['times'])-1 for c in records for h in records[c][2])//6)*6+1
    common_times=records['reference'][2][0]['times'][:count]
    for c in records:
        for h in records[c][2]:np.testing.assert_allclose(h['times'][:count],common_times,rtol=1e-12,atol=1e-17)
    for stride in [1,2,3]:
        for case in records:
            histories=[dict(h,times=h['times'][:count:stride],values=h['values'][:count:stride]) for h in records[case][2]]
            results[case]=rt.analyze_ports(histories,frequencies,'guided',(C,rho),bank)['T']
        actual=np.diagonal(results['group'],axis1=1,axis2=2)/np.diagonal(results['reference'],axis1=1,axis2=2)
        delta=abs(abs(actual)-abs(target))
        relative=delta/np.maximum(np.maximum(abs(1-abs(actual)),abs(1-abs(target))),1e-12)
        row=dict(output_stride=stride,magnitude_delta_pp=float(delta.max()*100),relative_loss_percent=float(relative.max()*100),
                 center_loss_percent=((1-abs(actual[8]))*100).tolist(),normalized_T_real=actual.real.tolist(),normalized_T_imag=actual.imag.tolist())
        rows.append(row);print({k:v for k,v in row.items() if not k.startswith('normalized')},flush=True)
    rt.write_json(root/'sampling_validation.json',dict(frequencies=frequencies.tolist(),rows=rows,
                  comparison_end_time_s=float(common_times[-1]),identical_endpoints_verified=True,
                  one_vs_two_step_output_checks=cadence_controls(rows),
                  note='Same solved time histories and identical endpoints, decimated offline; no internal-step, geometry or observation-window change.'))


if __name__=='__main__':main()
