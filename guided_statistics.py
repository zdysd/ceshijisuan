"""Coherent transmission statistics and length-based modal attenuation fits."""
import json
from pathlib import Path
import numpy as np
from cfrp_pores.attenuation import fit_attenuation


def matrix(record,frequency):
    r=next(r for r in record['results'] if r['frequency_hz']==frequency)
    return np.asarray(r['transmission_real'])+1j*np.asarray(r['transmission_imag'])


def normalized_transmission(record,reference,frequency):
    return np.diag(matrix(record,frequency))/np.diag(matrix(reference,frequency))


def coherent_interval(signals,replicates=1000,seed=2026):
    data=np.asarray(signals,complex)
    if data.ndim!=1 or len(data)<3:raise ValueError('At least three independent complex transmissions required.')
    rng=np.random.default_rng(seed)
    delta=1-abs(data.mean())
    boot=[];typical=[]
    for _ in range(replicates):
        sample=data[rng.integers(len(data),size=len(data))]
        boot.append(1-abs(sample.mean()));typical.append(float(np.mean(1-abs(sample))))
    return dict(samples=len(data),amplitude_reduction=float(delta),
                ci95_amplitude_reduction=np.percentile(boot,[2.5,97.5]).tolist(),
                mean_individual_amplitude_reduction=float(np.mean(1-abs(data))),
                ci95_mean_individual_reduction=np.percentile(typical,[2.5,97.5]).tolist(),
                fraction_individuals_above_0p1percent=float(np.mean(1-abs(data)>.001)),
                individual_magnitude_standard_deviation=float(abs(data).std(ddof=1)))


def attenuation_interval(lengths_m,signals,replicates=1000):
    """Resample whole seed families across lengths; retain intercept/interface loss."""
    lengths=np.asarray(lengths_m,float);data=np.asarray(signals,complex)
    if data.ndim!=2 or data.shape[1]!=len(lengths) or data.shape[0]<3:
        raise ValueError('Expected at least 3 seed families x >=3 distinct lengths.')
    fit=fit_attenuation(lengths,data.mean(axis=0))
    typical_fit=fit_attenuation(lengths,np.exp(np.log(abs(data)).mean(axis=0)))
    rng=np.random.default_rng(2026)
    boot=[];typical=[]
    for _ in range(replicates):
        sample=data[rng.integers(len(data),size=len(data))]
        boot.append(fit_attenuation(lengths,sample.mean(axis=0))['alpha_np_m'])
        typical.append(fit_attenuation(lengths,np.exp(np.log(abs(sample)).mean(axis=0)))['alpha_np_m'])
    ci=np.percentile(boot,[2.5,97.5])
    return dict(samples=data.shape[0],lengths_m=lengths.tolist(),fit=fit,ci95_alpha_np_m=ci.tolist(),
                typical_fit=typical_fit,ci95_typical_alpha_np_m=np.percentile(typical,[2.5,97.5]).tolist(),
                positive_attenuation_identified=bool(ci[0]>0),exponential_fit_passed=bool(fit['r_squared']>=.95),
                bootstrap_unit='Matched random-seed family across all tested lengths',
                local_alpha_np_m_by_adjacent_interval=[float(-np.log(abs(data[:,i+1].mean()/data[:,i].mean()))/(lengths[i+1]-lengths[i])) for i in range(len(lengths)-1)])


def summarize(root):
    root=Path(root);reference=json.loads((root/'reference/results.json').read_text())
    rows=[]
    for d in [20,40,80]:
        records=[]
        for folder in sorted((root/f'd{d}/L8').glob('seed_*')):
            if (folder/'results.json').exists():records.append(json.loads((folder/'results.json').read_text()))
        if not records:continue
        for f in [50000.,100000.,200000.]:
            transmission=np.array([normalized_transmission(r,reference,f) for r in records])
            for index,mode in enumerate(['S0','A0']):
                row=dict(diameter_group_um=d,frequency_hz=f,mode=mode,samples=len(records),
                         amplitude_reduction=float(1-abs(transmission[:,index].mean())))
                if len(records)>=3:row.update(coherent_interval(transmission[:,index]))
                single=root/f'd{d}/single_maximum/results.json'
                if single.exists():
                    r=json.loads(single.read_text());row['largest_single_magnitude_change']=float(abs(1-abs(normalized_transmission(r,reference,f)[index])))
                    sr=next(v for v in r['results'] if v['frequency_hz']==f)
                    rr=next(v for v in reference['results'] if v['frequency_hz']==f)
                    reflection=(np.asarray(sr['reflection_real'])+1j*np.asarray(sr['reflection_imag']))-(np.asarray(rr['reflection_real'])+1j*np.asarray(rr['reflection_imag']))
                    row['largest_single_reflected_mode_change_max']=float(np.max(abs(reflection[:,index])))
                rows.append(row)
    (root/'sweep_summary.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
    return rows


if __name__=='__main__':
    print(json.dumps(summarize('output/guided_collective'),indent=2))
