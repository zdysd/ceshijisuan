"""Transient port spectra, direction separation and incident-normalized scattering.

Real histories use Fourier exp(-i omega t), consistent with COMSOL exp(+i omega t).
All gates retain absolute time. Frequencies outside the source band are excluded.
"""
from functools import lru_cache
import numpy as np
from .dispersion import modal_direction_basis, separate_modal_directions


def minimum_inplane_speed(stiffness,density):
    """Minimum phase speed from the 2D Christoffel tensor, not Voigt eigenvalues.

    The engineering-shear Voigt matrix has a different strain metric; its
    minimum eigenvalue need not give the slowest oblique shear wave.
    """
    from scipy.optimize import minimize_scalar
    C=np.asarray(stiffness,float)
    if C.shape!=(3,3) or density<=0:raise ValueError('Expected plane-strain 3x3 stiffness and positive density.')
    def eigenvalue(angle):
        x,y=np.cos(angle),np.sin(angle)
        a=C[0,0]*x*x+2*C[0,2]*x*y+C[2,2]*y*y
        d=C[2,2]*x*x+2*C[1,2]*x*y+C[1,1]*y*y
        b=C[0,2]*x*x+(C[0,1]+C[2,2])*x*y+C[1,2]*y*y
        return .5*(a+d-np.sqrt((a-d)**2+4*b*b))
    angles=np.linspace(0,np.pi,361);values=eigenvalue(angles);index=int(np.argmin(values))
    step=np.pi/360
    refined=minimize_scalar(eigenvalue,bounds=(angles[index]-step,angles[index]+step),method='bounded')
    minimum=min(float(values[index]),float(refined.fun))
    if minimum<=0:raise ValueError('Nonpositive acoustic tensor eigenvalue.')
    return float(np.sqrt(minimum/density))


def hann_pulse(times, frequency, cycles=5):
    t=np.asarray(times); duration=cycles/frequency
    return np.where((t>=0)&(t<=duration),.5*(1-np.cos(2*np.pi*t/duration))*np.sin(2*np.pi*frequency*t),0.)


def continuation_block_ends(planned,completed):
    """Keep saved block endpoints when extending a formerly final short block."""
    planned=list(planned);completed=list(completed)
    if not planned or any(b<=a for a,b in zip([0]+planned[:-1],planned)):
        raise ValueError('Planned endpoints must be positive and increasing.')
    if completed and (completed[-1]>planned[-1] or any(b<=a for a,b in zip([0]+completed[:-1],completed))):
        raise ValueError('Saved endpoints exceed the requested duration or are not increasing.')
    return completed+[n for n in planned if not completed or n>completed[-1]]


def bulk_window_limit(case):
    return 50e-6 if case=='reference' else 100e-6


def next_bulk_window(end_time,limit=50e-6):
    """Double an observation window without exceeding its hard limit."""
    if not np.isfinite(end_time) or end_time<=0:raise ValueError('Positive finite window required.')
    if not np.isfinite(limit) or limit<=0:raise ValueError('Positive finite limit required.')
    return None if end_time>=limit-1e-15 else min(2*end_time,limit)


def bulk_adjustment(gate,rows):
    """Choose a standard refinement only after complete P controls exist."""
    checks=[c for c in gate['checks'] if c['kind']=='bulk' and c['case'] in ['reference','group']]
    if len(checks)!=2 or not all(c['complete'] for c in checks):
        raise ValueError('Complete both P reference/group and half-step controls before choosing a refinement.')
    limits=gate['thresholds'];windows=[]
    for r in rows:
        if r['kind']=='bulk' and r['variant'] in ['base','time'] and r['case'] in ['reference','group']:
            if r['tail_energy_fraction']>=limits['tail_energy'] or r['window_delta_max']>=limits['absolute_amplitude']:
                windows.append([r['variant'],r['case']])
    if windows:return dict(action='extend_window',cases=windows)
    if all(c['passed'] for c in checks):return dict(action='accepted',cases=[])
    for c in checks:
        temporal=c['differences'].get('time',{}).get('band',{})
        frequency=c.get('frequency_band_consistency') or {}
        if frequency.get('passed') is False or c.get('harmonic_magnitude_delta',0)>=limits['absolute_amplitude'] or \
           c.get('harmonic_diagonal_phase_error_deg',0)>=limits['phase_degrees'] or \
           temporal.get('magnitude_delta',0)>limits['absolute_amplitude'] or \
           temporal.get('phase_error_deg',0)>limits['phase_degrees'] or \
           (c['case']=='group' and temporal.get('relative_loss_change',0)>limits['relative_nonzero_loss']):
            return dict(action='refine_time',cases=[['base','reference'],['base','group'],['time','reference'],['time','group']])
    raise RuntimeError('P acceptance failed for a reason requiring analysis; do not blindly refine or change thresholds.')


def field_frame_indices(times,fmax):
    """Keep at least eight frames per shortest resolved temporal period."""
    times=np.asarray(times,float)
    if len(times)<2 or fmax<=0 or np.any(np.diff(times)<=0):raise ValueError('Increasing times and positive bandwidth required.')
    stride=max(1,int(np.floor(1/(8*fmax*np.max(np.diff(times))))))
    count=min(len(times),max(96,int(np.ceil((len(times)-1)/stride))+1))
    return np.unique(np.linspace(0,len(times)-1,count,dtype=int))


def source_band(frequency, cycles=5, energy_fraction=.999):
    """Resolve 99.9% of energy of both quadratures of the finite pulse."""
    dt=1/(256*frequency); t=np.arange(0,cycles/frequency,dt)
    envelope=.5*(1-np.cos(2*np.pi*t*frequency/cycles))
    n=32768
    power=sum(abs(np.fft.rfft(envelope*fn(2*np.pi*frequency*t),n))**2 for fn in [np.sin,np.cos])
    f=np.fft.rfftfreq(n,dt); c=np.cumsum(power)/power.sum()
    return float(f[np.searchsorted(c,(1-energy_fraction)/2)]),float(f[np.searchsorted(c,(1+energy_fraction)/2)])


def fourier_history(times, values, frequencies, end_fraction=1.):
    t=np.asarray(times,float); v=np.asarray(values)
    if t.ndim!=1 or v.shape[0]!=len(t) or np.any(np.diff(t)<=0):
        raise ValueError('Expected increasing times and time-first histories.')
    if not 0<end_fraction<=1:raise ValueError('Invalid gate fraction.')
    keep=t<=t[0]+end_fraction*(t[-1]-t[0]); t=t[keep];v=v[keep]
    if len(t)<8:raise ValueError('Too few samples in time gate.')
    weights=np.empty(len(t));weights[0]=(t[1]-t[0])/2;weights[-1]=(t[-1]-t[-2])/2
    weights[1:-1]=(t[2:]-t[:-2])/2
    kernel=np.exp(-2j*np.pi*np.asarray(frequencies)[:,None]*t)*weights
    return np.tensordot(kernel,v,axes=(1,0))


def scattering(incoming,outgoing,reflected,right_incoming,max_condition=100.):
    """Columns are independent launches; rows are physical modes."""
    incoming=np.asarray(incoming)
    condition=float(np.linalg.cond(incoming))
    if not np.isfinite(condition) or condition>max_condition:
        raise ValueError('Incident modal matrix is ill-conditioned.')
    def solve(a):return np.linalg.solve(incoming.T,np.asarray(a).T).T
    return solve(outgoing),solve(reflected),solve(right_incoming),condition


def modal_bank(plate,frequencies,y):
    plate.at_frequency=lru_cache(maxsize=256)(plate.at_frequency)
    bank=[]
    for frequency in frequencies:
        if frequency>=plate.first_higher_cutoff_hz():raise ValueError('Pulse band includes higher Lamb modes.')
        basis,labels,scale=modal_direction_basis(y-plate.thickness_m/2,plate,frequency)
        modes=[plate.at_frequency(frequency,m) for m in ['S0','A0']]
        power=[]
        for col in [0,2]:
            u,v,sx,sy=basis[:,col].reshape(4,len(y));sx=sx*scale;sy=sy*scale
            power.append(float(np.trapezoid(-.5*np.real(sx*np.conj(2j*np.pi*frequency*u)+sy*np.conj(2j*np.pi*frequency*v)),y)))
        if min(power)<=0:raise ValueError('Invalid modal power direction.')
        bank.append((basis,labels,scale,modes,np.array(power)))
    return bank


def analyze_ports(histories,frequencies,kind,material,bank=None,end_fraction=1.,spectra_override=None):
    """history fields: times, values[t,port,component,y], y, positions."""
    C,rho=material; spectra=[]; residuals=[]
    if spectra_override is not None:spectra=spectra_override
    else:
        for history in histories:
            spectra.append(fourier_history(history['times'],history['values'],frequencies,end_fraction))
    positions=histories[0]['positions'];y=histories[0]['y'];npre=len(positions)//2
    matrices=[];reflections=[];rights=[];conditions=[];energies=[];wave_lengths=[];powers=[]
    for fi,f in enumerate(frequencies):
        launches=[];res=[]
        if kind=='guided':
            basis,labels,scale,modes,power=bank[fi]
            k=np.array([m['k_rad_m'] for m in modes]);wave_lengths.append([m['wavelength_m'] for m in modes])
        else:
            speed=np.sqrt(C[0,0]/rho);k=np.array([2*np.pi*f/speed]);power=np.ones(1)
            wave_lengths.append([speed/f])
        for spectrum in spectra:
            if kind=='guided':
                coeff=[]
                for port in spectrum[fi]:
                    a,r=separate_modal_directions(y- (y[0]+y[-1])/2,*port,basis,labels,scale)
                    coeff.append([a[label] for label in labels]);res.append(r)
                coeff=np.array(coeff);forward=coeff[:,[0,2]];backward=coeff[:,[1,3]]
            else:
                mean=np.trapezoid(spectrum[fi],y,axis=-1)/np.ptp(y)
                u=mean[:,0];stress=mean[:,2];Z=np.sqrt(rho*C[0,0])
                forward=(.5*(u-stress/(2j*np.pi*f*Z)))[:,None]
                backward=(.5*(u+stress/(2j*np.pi*f*Z)))[:,None]
            def avg(a,sl,sign):
                return np.mean(a[sl]*np.exp(sign*1j*positions[sl,None]*k),axis=0)
            launches.append((avg(forward,slice(0,npre),1),avg(forward,slice(npre,None),1),
                             avg(backward,slice(0,npre),-1),avg(backward,slice(npre,None),-1)))
        T,R,right,cond=scattering(*[np.column_stack([v[i] for v in launches]) for i in range(4)])
        matrices.append(T);reflections.append(R);rights.append(right);conditions.append(cond)
        powers.append(power)
        energies.append(np.sum((power[:,None]/power[None,:])*(abs(T)**2+abs(R)**2),axis=0))
        residuals.append(max(res,default=0.))
    return dict(frequencies=np.asarray(frequencies),T=np.array(matrices),R=np.array(reflections),
                right=np.array(rights),condition=np.array(conditions),residual=np.array(residuals),
                power_balance=np.array(energies),wavelengths=np.array(wave_lengths),power_flux=np.array(powers))


def tail_fraction(history):
    """Last-10% fraction of the displacement-squared signal integral, not elastic energy."""
    t=history['times'];v=history['values'][:,:,:2]
    power=np.sum(v*v,axis=(2,3));total=np.trapezoid(power,t,axis=0)
    sel=t>=t[0]+.9*(t[-1]-t[0])
    return float(np.max(np.trapezoid(power[sel],t[sel],axis=0)/np.maximum(total,1e-100)))
