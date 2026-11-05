"""Compare single and averaged receiving sections using existing modal data."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'output/guided_collective'
OUT=ROOT/'output/receiving_section_comparison';OUT.mkdir(exist_ok=True)
FREQUENCY=200000


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def transmission(folder,pre=(0,1,2),post=(3,4,5)):
    record=read(folder/'results.json')
    row=next(r for r in record['results'] if r['frequency_hz']==FREQUENCY)
    k=2*np.pi/np.asarray(row['wavelength_m'])
    with np.load(folder/'signals.npz') as data:
        position=data['positions_m'];incoming=[];outgoing=[]
        for mode in ['S0','A0']:
            c=data[f'{FREQUENCY}_{mode}_coefficients'][:,[0,2]]
            referenced=c*np.exp(1j*position[:,None]*k[None,:])
            incoming.append(referenced[list(pre)].mean(axis=0))
            outgoing.append(referenced[list(post)].mean(axis=0))
    result=np.column_stack(outgoing)@np.linalg.inv(np.column_stack(incoming))
    if len(pre)==3 and len(post)==3:
        original=np.array(row['transmission_real'])+1j*np.array(row['transmission_imag'])
        assert np.max(abs(result-original))<1e-10
    return result


def folder(length,seed):
    base=DATA if seed<=106 else DATA/'extended_samples'
    return base/f'd80/L{length}/seed_{seed}'


def main():
    lengths=np.array([8,24,48]);results={}
    for name,pre,post in [('three_sections',(0,1,2),(3,4,5)),('middle_section',(1,),(4,))]:
        families=[]
        for seed in range(101,113):
            values=[]
            for length in lengths:
                ref=DATA/('reference' if length==8 else f'reference_L{length}')
                values.append(np.diag(transmission(folder(length,seed),pre,post))/np.diag(transmission(ref,pre,post)))
            families.append(values)
        means=np.asarray(families).mean(axis=0)
        results[name]={mode:dict(longest_amplitude_reduction_percent=float((1-abs(means[-1,i]))*100),
                                alpha_np_m=float(-np.polyfit(lengths*.001,np.log(abs(means[:,i])),1)[0])) for i,mode in enumerate(['S0','A0'])}
    record={}
    for mode in ['S0','A0']:
        a=results['three_sections'][mode];b=results['middle_section'][mode]
        record[mode]=dict(three_sections=a,middle_section=b,
                         amplitude_difference_percentage_points=abs(a['longest_amplitude_reduction_percent']-b['longest_amplitude_reduction_percent']),
                         alpha_relative_difference=abs(a['alpha_np_m']-b['alpha_np_m'])/abs(a['alpha_np_m']))
    output=dict(frequency_hz=FREQUENCY,diameter_group_um=80,porosity=.015,samples=12,
                lengths_mm=lengths.tolist(),comparison=record,
                note='Existing modal coefficients only; identical phase correction, reference normalization and two independent launches.')
    (OUT/'single_vs_three_sections.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
