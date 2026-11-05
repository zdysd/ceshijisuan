"""Test attenuation versus porous length, with matched uniform-port references."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import argparse,json,sys
import numpy as np
from guided_collective import settings,OUT,run_case,grouped_pores,export_geometry,read
from guided_statistics import normalized_transmission,attenuation_interval
sys.stdout.reconfigure(encoding='utf-8')


def length_settings(s,length_mm):
    s=dict(s)
    if length_mm>16:
        end=s['start']+length_mm*.001
        s['post_positions']=[end+.008,end+.010,end+.012]
    if max(s['post_positions'])>=s['sponge']:
        raise ValueError('Receiving port enters the absorber; length is too large for this layout.')
    return s


def summarize(samples):
    rows=[]
    for f in [100000.,200000.]:
        for index,mode in enumerate(['S0','A0']):
            signals=[];seeds=[]
            for seed in range(101,101+samples):
                values=[]
                for length in [8,24,48]:
                    folder=OUT/f'd80/L{length}/seed_{seed}'
                    ref=OUT/('reference' if length==8 else f'reference_L{length}')
                    if not (folder/'results.json').exists():break
                    values.append(normalized_transmission(read(folder/'results.json'),read(ref/'results.json'),f)[index])
                if len(values)==3:signals.append(values);seeds.append(seed)
            if len(signals)>=3:
                result=attenuation_interval(np.array([8,24,48])*.001,signals)
                rows.append(dict(diameter_group_um=80,frequency_hz=f,mode=mode,seeds=seeds,**result))
    (OUT/'length_attenuation.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('--samples',type=int,default=6);p.add_argument('--cores',type=int,default=6)
    p.add_argument('--summarize-only',action='store_true');args=p.parse_args()
    if not args.summarize_only:
        s=settings();client=mph.Client(cores=args.cores,version='6.2')
        for length in [24,48]:
            local=length_settings(s,length)
            run_case(client,OUT/f'reference_L{length}',local,frequencies=[100000.,200000.])
        for seed in range(101,101+args.samples):
            for length in [24,48]:
                local=length_settings(s,length);folder=OUT/f'd80/L{length}/seed_{seed}'
                pores,meta=grouped_pores(80,length,seed);export_geometry(folder/'geometry',pores,meta)
                run_case(client,folder,local,pores,frequencies=[100000.,200000.])
                summarize(args.samples)
    print(json.dumps(summarize(args.samples),indent=2),flush=True)


if __name__=='__main__':main()
