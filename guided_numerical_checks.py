"""Matched finer meshes and geometric tolerances for detectability/alpha gates."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import json,sys
import numpy as np
from guided_collective import OUT,settings,run_case,grouped_pores,read
from guided_length_study import length_settings
from guided_statistics import normalized_transmission
from cfrp_pores.attenuation import fit_attenuation
sys.stdout.reconfigure(encoding='utf-8')


def main():
    s=settings();client=mph.Client(cores=6,version='6.2');out=OUT/'checks';rows=[]
    fine_reference=run_case(client,out/'reference_fine',s,mesh=.000075,repair_tolerance_m=1e-9)
    original_reference=read(OUT/'reference/results.json')
    for d in [20,40,80]:
        pores,meta=grouped_pores(d,8,101);largest=max(pores,key=lambda p:p.area)
        fine=run_case(client,out/f'd{d}_maximum_fine',s,[largest],mesh=.000075,repair_tolerance_m=1e-9,micro_divisor=60)
        original=read(OUT/f'd{d}/single_maximum/results.json')
        for f in s['frequencies']:
            coarse_T=normalized_transmission(original,original_reference,f)
            fine_T=normalized_transmission(fine,fine_reference,f)
            for index,mode in enumerate(['S0','A0']):
                error=float(abs(abs(fine_T[index])-abs(coarse_T[index])))
                rows.append(dict(case=f'd{d}_maximum',mode=mode,frequency_hz=f,
                                 coarse_magnitude_change=float(abs(1-abs(coarse_T[index]))),
                                 fine_magnitude_change=float(abs(1-abs(fine_T[index]))),
                                 mesh_magnitude_error=error,
                                 numerical_error_below_primary_threshold=error<.001,
                                 primary_undetectability_with_mesh_margin=bool(abs(1-abs(fine_T[index]))+error<.001)))
    length_transmission={f:[] for f in [100000.,200000.]}
    for length in [8,24,48]:
        local=length_settings(s,length)
        ref=fine_reference if length==8 else run_case(client,out/f'reference_L{length}_fine',local,
            frequencies=[100000.,200000.],mesh=.000075,repair_tolerance_m=1e-9)
        coarse_ref=original_reference if length==8 else read(OUT/f'reference_L{length}/results.json')
        pores,meta=grouped_pores(80,length,101)
        fine=run_case(client,out/f'd80_L{length}_fine',local,pores,frequencies=[100000.,200000.],
                      mesh=.000075,repair_tolerance_m=1e-9,micro_divisor=60)
        coarse=read(OUT/f'd80/L{length}/seed_101/results.json')
        for f in [100000.,200000.]:
            fine_T=normalized_transmission(fine,ref,f);coarse_T=normalized_transmission(coarse,coarse_ref,f)
            length_transmission[f].append((coarse_T,fine_T))
            for index,mode in enumerate(['S0','A0']):
                rows.append(dict(case=f'd80_L{length}',mode=mode,frequency_hz=f,
                                 coarse_amplitude_reduction=float(1-abs(coarse_T[index])),
                                 fine_amplitude_reduction=float(1-abs(fine_T[index])),
                                 mesh_magnitude_error=float(abs(abs(fine_T[index])-abs(coarse_T[index])))))
    alpha=[]
    for f,pairs in length_transmission.items():
        for index,mode in enumerate(['S0','A0']):
            coarse=fit_attenuation(np.array([8,24,48])*.001,[a[index] for a,b in pairs])
            fine=fit_attenuation(np.array([8,24,48])*.001,[b[index] for a,b in pairs])
            change=abs(coarse['alpha_np_m']-fine['alpha_np_m'])/max(abs(coarse['alpha_np_m']),abs(fine['alpha_np_m']),1e-30)
            alpha.append(dict(frequency_hz=f,mode=mode,coarse_fit=coarse,fine_fit=fine,
                              relative_alpha_change=change,relative_alpha_mesh_passed=bool(change<=.05)))
    result=dict(transmission_checks=rows,alpha_checks=alpha,
                note='Fixed seed101 numerical comparisons; ensemble uncertainty is evaluated separately.')
    (OUT/'numerical_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
