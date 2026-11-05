"""Inspect transverse Fourier channels at the failed bulk absorber frequencies."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def main():
    root=ROOT/'output/transient_converged/spatial_checks/bulk'
    cfg=json.loads((ROOT/'config.json').read_text());C=np.array(cfg['stiffness_GPa_voigt_xx_yy_zz_yz_xz_xy'])*1e9
    rho=cfg['density_kg_m3'];cp=np.sqrt(C[0,0]/rho);Z=rho*cp
    data={v:dict(np.load(root/v/'group/P_harmonic_band.npz')) for v in ['base','absorber']}
    normalized={v:np.load(root/v/'group/spectra_band.npz')['T'][:,0,0]/np.load(root/v/'reference/spectra_band.npz')['T'][:,0,0] for v in ['base','absorber']}
    frequencies=data['base']['frequencies'];rows=[]
    for fi,f in enumerate(frequencies):
        omega=2*np.pi*f;scale=omega*Z;row=dict(frequency_hz=float(f),magnitude_delta_percentage_points=float(100*abs(abs(normalized['base'][fi])-abs(normalized['absorber'][fi]))))
        for variant,h in data.items():
            y=h['y'];width=.001;ports=h['values'][fi];mean=np.trapezoid(ports,y,axis=-1)/np.ptp(y)
            incoming=.5*(mean[:3,0]-mean[:3,2]/(1j*scale))
            amplitude=np.mean(incoming*np.exp(1j*omega/cp*h['positions'][:3]))
            input_power=.5*Z*omega**2*abs(amplitude)**2;reflected_power=[];decays=[];velocities=[]
            for n in range(-5,6):
                ky=2*np.pi*n/width
                a=C[0,0]*C[5,5]
                b=(C[0,0]*C[1,1]+C[5,5]**2-(C[0,1]+C[5,5])**2)*ky**2-rho*omega**2*(C[0,0]+C[5,5])
                c=(C[5,5]*ky**2-rho*omega**2)*(C[1,1]*ky**2-rho*omega**2)
                ks=[sign*np.sqrt(complex(z)) for z in np.roots([a,b,c]) for sign in [1,-1]]
                basis=[];flux=[]
                for k in ks:
                    matrix=np.array([[C[0,0]*k*k+C[5,5]*ky*ky-rho*omega**2,-(C[0,1]+C[5,5])*k*ky],
                                     [-(C[0,1]+C[5,5])*k*ky,C[5,5]*k*k+C[1,1]*ky*ky-rho*omega**2]],complex)
                    _,_,vh=np.linalg.svd(matrix);u,v=vh.conj().T[:,-1]
                    sx=-1j*C[0,0]*k*u+1j*C[0,1]*ky*v;sy=1j*C[5,5]*(ky*u-k*v)
                    basis.append([u,v,sx/scale,sy/scale])
                    power=float(-.5*np.real(sx*np.conj(1j*omega*u)+sy*np.conj(1j*omega*v)))
                    propagating=abs(k.imag)<1e-6*max(abs(k),1)
                    flux.append(power if propagating else 0.)
                    if propagating and abs(power)>0:velocities.append(abs(power/(.5*rho*omega**2*(abs(u)**2+abs(v)**2))))
                    if not propagating and abs(k.imag)>1e-8:decays.append(1/abs(k.imag))
                basis=np.array(basis).T
                fourier=np.trapezoid(ports[3:]*np.exp(-1j*ky*y)[None,None,:],y,axis=-1)/width
                fourier[:,2:]/=scale
                coefficients=np.linalg.solve(basis,fourier.T)
                backward=np.sum(np.maximum(-np.array(flux),0)[:,None]*abs(coefficients)**2,axis=0)/input_power
                reflected_power.append(backward)
            row[variant]=dict(incoming_power_at_right_ports_percent=(100*np.sum(reflected_power,axis=0)).tolist(),
                longest_evanescent_decay_mm=float(max(decays,default=0)*1000),minimum_propagating_group_speed_x_m_s=float(min(velocities,default=0)))
        rows.append(row)
    output=ROOT/'output/transient_converged/bulk_absorber_diagnosis.json'
    output.write_text(json.dumps(rows,indent=2),encoding='utf-8')
    print(json.dumps(sorted(rows,key=lambda r:r['magnitude_delta_percentage_points'],reverse=True)[:5],indent=2))


if __name__=='__main__':main()
