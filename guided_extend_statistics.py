"""Extend the most promising 200 kHz/80 um length comparison to 12 seed families."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import argparse,json,sys
import numpy as np
from guided_collective import OUT,settings,run_case,grouped_pores,export_geometry,read
from guided_length_study import length_settings
from guided_statistics import normalized_transmission,attenuation_interval,coherent_interval
sys.stdout.reconfigure(encoding='utf-8')


def family_path(length,seed):
    return OUT/f'd80/L{length}/seed_{seed}' if seed<=106 else OUT/f'extended_samples/d80/L{length}/seed_{seed}'


def summarize(target=12):
    results=[]
    for index,mode in enumerate(['S0','A0']):
        values=[];seeds=[]
        for seed in range(101,101+target):
            family=[]
            for length in [8,24,48]:
                path=family_path(length,seed)/'results.json'
                if not path.exists():break
                reference=OUT/('reference' if length==8 else f'reference_L{length}')/'results.json'
                family.append(normalized_transmission(read(path),read(reference),200000.)[index])
            if len(family)==3:values.append(family);seeds.append(seed)
        if len(values)<3:continue
        values=np.asarray(values)
        result=dict(frequency_hz=200000.,mode=mode,diameter_group_um=80,seeds=seeds,
                    **attenuation_interval(np.array([8,24,48])*.001,values))
        result['longest_segment_contrast']=coherent_interval(values[:,-1])
        if len(values)>=6:
            first=attenuation_interval(np.array([8,24,48])*.001,values[:len(values)//2],100)
            alpha=result['fit']['alpha_np_m'];previous=first['fit']['alpha_np_m']
            relative=abs(alpha-previous)/max(abs(alpha),abs(previous),1e-30)
            result['sample_count_relative_change']=relative
            result['sample_count_5percent_passed']=relative<=.05
        results.append(result)
    (OUT/'extended_attenuation.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    return results


def main():
    p=argparse.ArgumentParser();p.add_argument('--target',type=int,default=12);p.add_argument('--summarize-only',action='store_true');args=p.parse_args()
    if not args.summarize_only:
        s=settings();client=mph.Client(cores=2,version='6.2')
        for seed in range(107,101+args.target):
            for length in [8,24,48]:
                local=length_settings(s,length);folder=family_path(length,seed)
                pores,meta=grouped_pores(80,length,seed);export_geometry(folder/'geometry',pores,meta)
                run_case(client,folder,local,pores,frequencies=[200000.])
            summarize(args.target)
    print(json.dumps(summarize(args.target),indent=2),flush=True)


if __name__=='__main__':main()
