"""Separate receiver quadrature error from FE mesh error on the same solved field."""
import run_transient as rt
import run_refined_pipeline as pipeline
from cfrp_pores.comsol import integrate_section
import numpy as np,json


def main():
    root=rt.ROOT/'output/transient_refined/port_quadrature'
    client=rt.start_session(root,4)
    C,rho=rt.validate_material(json.loads((rt.ROOT/'config.json').read_text()))
    original=rt.frequency_checkpoint
    def checkpoint(model,s,folder,mode,frequencies=None):
        original(model,s,folder,mode,frequencies)
        results={}
        for order in [8,12]:
            mean=np.array([integrate_section(model,x,s['width'],['u','v','solid.sx','solid.sxy'],order)/s['width'] for x in s['positions']])
            h=dict(values=np.repeat(mean[:,:,None],2,axis=-1),y=np.array([0.,s['width']]),positions=np.array(s['positions']))
            r=rt.analyze_ports([h],[s['f0']],'bulk',(C,rho),spectra_override=[h['values'][None]])
            results[f'FE_quadrature_{order}']=r['T'][0,0,0]
        y=np.linspace(1e-8,s['width']-1e-8,2049);xx,yy=np.meshgrid(s['positions'],y,indexing='ij')
        fields=rt.interpolate(model,[xx.ravel(),yy.ravel()],['u','v','solid.sx','solid.sxy'])[:,0,:]
        values=fields.T.reshape(6,len(y),4).transpose(0,2,1)
        for stride in [16,4,1]:
            h=dict(values=values[:,:,::stride],y=y[::stride],positions=np.array(s['positions']))
            r=rt.analyze_ports([h],[s['f0']],'bulk',(C,rho),spectra_override=[h['values'][None]])
            results[f'trapezoid_{len(h["y"])}']=r['T'][0,0,0]
        rt.write_json(folder/'quadrature_comparison.json',{k:dict(real=float(v.real),imag=float(v.imag),magnitude=float(abs(v))) for k,v in results.items()})
        print('QUADRATURE',s['variant'],{k:float(abs(v)) for k,v in results.items()},flush=True)
    rt.frequency_checkpoint=checkpoint
    for variant in ['base','mesh']:
        pipeline.center_case(client,root,'bulk',8.,101,'group',variant,profile='refined_v3')


if __name__=='__main__':main()
