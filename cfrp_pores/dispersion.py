"""Through-thickness SAFE discretization of 2D plane-strain elastic plates.

This resolves the plate thickness explicitly; it is not a shell theory. The
first S/A branches are obtained from the Hermitian eigenproblem at real k.
Full 3D anisotropic coupling and nonsymmetric laminates are deliberately rejected.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq


def validate_material(config):
    C=np.asarray(config['stiffness_GPa_voigt_xx_yy_zz_yz_xz_xy'],float)*1e9
    rho=float(config['density_kg_m3'])
    if C.shape!=(6,6) or not np.isfinite(C).all() or not np.allclose(C,C.T) or np.linalg.eigvalsh(C).min()<=0 or rho<=0:
        raise ValueError('Positive density and symmetric positive-definite 6x6 stiffness required.')
    # In-plane normal/shear block in global x,y; engineering shear strain.
    ix=[0,1,5]; other=[2,3,4]
    if np.max(np.abs(C[np.ix_([0,1,5],[3,4])]))>1e-8*np.max(C):
        raise ValueError('Out-of-plane shear coupling needs full vector 3D/SAFE physics.')
    return C,rho


@dataclass
class Plate:
    thickness_m: float
    stiffness: np.ndarray  # plane-strain Cxx,Cyy,Cxy,Gxy block
    density: float
    elements: int = 60
    layers: list | None = None  # [(thickness fraction,C3x3,rho),...], symmetric

    def matrices(self):
        if self.thickness_m<=0 or self.elements<8 or self.density<=0:
            raise ValueError('Positive thickness/density and >=8 thickness elements required.')
        layers=self.layers or [(1.,self.stiffness,self.density)]
        if not np.isclose(sum(layer[0] for layer in layers),1): raise ValueError('Layer fractions must sum to 1.')
        for frac,C,rho in layers:
            if frac<=0 or rho<=0 or np.linalg.eigvalsh(C).min()<=0: raise ValueError('Invalid layer material.')
        for a,b in zip(layers,reversed(layers)):
            if not np.isclose(a[0],b[0]) or not np.allclose(a[1],b[1]) or not np.isclose(a[2],b[2]):
                raise ValueError('A0/S0 classification requires a symmetric plate.')
        self.y=np.linspace(-self.thickness_m/2,self.thickness_m/2,self.elements+1)
        n=2*len(self.y)
        K0=np.zeros((n,n),complex); K1=K0.copy(); K2=K0.copy(); M=K0.copy()
        edges=np.cumsum([layer[0] for layer in layers])
        for e in range(self.elements):
            dy=self.y[e+1]-self.y[e]
            layer=layers[min(np.searchsorted(edges,(e+.5)/self.elements),len(layers)-1)]
            _,C,rho=layer
            indices=np.array([2*e,2*e+1,2*e+2,2*e+3])
            dN=np.array([-1,1])/dy
            for xi in [-1/np.sqrt(3),1/np.sqrt(3)]:
                N=np.array([1-xi,1+xi])/2
                B0=np.zeros((3,4),complex); B1=B0.copy(); Nm=np.zeros((2,4))
                B0[1,1::2]=dN; B0[2,0::2]=dN
                B1[0,0::2]=1j*N; B1[2,1::2]=1j*N
                Nm[0,0::2]=N; Nm[1,1::2]=N
                idx=np.ix_(indices,indices); weight=dy/2
                K0[idx]+=weight*B0.conj().T@C@B0
                K1[idx]+=weight*(B0.conj().T@C@B1+B1.conj().T@C@B0)
                K2[idx]+=weight*B1.conj().T@C@B1
                M[idx]+=weight*rho*Nm.T@Nm
        self.K0,self.K1,self.K2,self.M=K0,K1,K2,M
        return self

    def mode(self,k,kind):
        if not hasattr(self,'M'): self.matrices()
        if kind not in ('A0','S0'): raise ValueError('Mode must be A0 or S0.')
        vals,vecs=eigh(self.K0+k*self.K1+k*k*self.K2,self.M,subset_by_index=[0,min(11,len(self.M)-1)])
        for val,vec in zip(vals,vecs.T):
            if val < -1e-5*max(abs(vals[-1]),1): continue
            reflected=vec.reshape(-1,2)[::-1].copy(); reflected[:,1]*=-1
            parity=(np.vdot(vec,self.M@reflected.ravel())/np.vdot(vec,self.M@vec)).real
            if (kind=='S0' and parity>.9) or (kind=='A0' and parity<-.9):
                # Fix arbitrary eigenvector phase for reproducible excitation.
                v=vec*np.exp(-1j*np.angle(vec[np.argmax(np.abs(vec))]))
                v/=np.max(np.abs(v))
                return float(np.sqrt(max(0,val))),v,parity
        raise RuntimeError('Cannot identify requested fundamental mode; check symmetry and resolution.')

    def first_higher_cutoff_hz(self):
        if not hasattr(self,'M'):self.matrices()
        vals=eigh(self.K0,self.M,eigvals_only=True,subset_by_index=[0,5])
        # Two rigid in-plane translations are the A0/S0 zero-frequency limits.
        return float(np.sqrt(max(0,vals[2]))/(2*np.pi))

    def at_frequency(self,f_hz,kind):
        if f_hz<=0: raise ValueError('Frequency must be positive.')
        target=2*np.pi*f_hz
        # Extremely small k makes the two rigid-translation eigenvalues
        # numerically degenerate and their parity indeterminate.
        low=1.
        high=max(10.,target/np.sqrt(np.linalg.eigvalsh(self.stiffness).min()/self.density)*2)
        while self.mode(high,kind)[0]<target:
            high*=2
        k=brentq(lambda k:self.mode(k,kind)[0]-target,low,high,xtol=1e-7)
        omega,profile,parity=self.mode(k,kind)
        dk=k*1e-3
        cg=(self.mode(k+dk,kind)[0]-self.mode(k-dk,kind)[0])/(2*dk)
        return dict(frequency_hz=f_hz,mode=kind,k_rad_m=k,wavelength_m=2*np.pi/k,
                    phase_velocity_m_s=omega/k,group_velocity_m_s=cg,parity=float(parity),
                    y_m=self.y.copy(),u=profile[0::2],v=profile[1::2])


def project_mode(y_m, u, v, mode):
    """Complex displacement overlap; contains both components, not peak height."""
    y=np.asarray(y_m); u=np.asarray(u); v=np.asarray(v)
    if len(y)<8 or np.any(np.diff(y)<=0): raise ValueError('At least eight increasing thickness samples required.')
    mr=np.interp(y,mode['y_m'],mode['u'].real)+1j*np.interp(y,mode['y_m'],mode['u'].imag)
    mv=np.interp(y,mode['y_m'],mode['v'].real)+1j*np.interp(y,mode['y_m'],mode['v'].imag)
    norm=np.trapezoid(abs(mr)**2+abs(mv)**2,y)
    return np.trapezoid(np.conj(mr)*u+np.conj(mv)*v,y)/norm


def modal_direction_basis(y_m,plate,frequency):
    """Reference modal displacement/traction basis, exp(+iwt), +/-x."""
    y=np.asarray(y_m);C=plate.stiffness
    Z=np.sqrt(plate.density*C[0,0]);scale=2*np.pi*frequency*Z
    columns=[];labels=[]
    for kind in ['S0','A0']:
        d=plate.at_frequency(frequency,kind)
        for direction in ['forward','backward']:
            q=-d['k_rad_m'] if direction=='forward' else d['k_rad_m']
            u=np.conj(d['u']) if direction=='forward' else d['u']
            v=np.conj(d['v']) if direction=='forward' else d['v']
            du=np.gradient(u,d['y_m'],edge_order=2);dv=np.gradient(v,d['y_m'],edge_order=2)
            strains=np.array([1j*q*u,dv,du+1j*q*v]);stress=C@strains
            values=[]
            for profile in [u,v,stress[0]/scale,stress[2]/scale]:
                values.append(np.interp(y,d['y_m'],profile.real)+1j*np.interp(y,d['y_m'],profile.imag))
            columns.append(np.concatenate(values));labels.append(kind+'_'+direction)
    return np.array(columns).T,labels,scale


def separate_modal_directions(y_m,displacement_u,displacement_v,stress_xx,stress_xy,basis,labels,scale):
    """Joint least-squares projection; report nonmodal residual instead of hiding it."""
    y=np.asarray(y_m);n=len(y)
    values=np.array([displacement_u,displacement_v,np.asarray(stress_xx)/scale,np.asarray(stress_xy)/scale])
    finite=np.isfinite(values).all(axis=0)
    if finite.sum()<8:raise ValueError('Insufficient solid thickness points for modal decomposition.')
    selection=np.tile(finite,4)
    weights=np.gradient(y)
    weights=np.sqrt(np.tile(weights,4))[selection]
    matrix=basis[selection]*weights[:,None];rhs=values.ravel()[selection]*weights
    if np.linalg.cond(matrix)>1e8:raise ValueError('Modal direction basis is ill-conditioned at this frequency.')
    amplitudes,_,_,_=np.linalg.lstsq(matrix,rhs,rcond=None)
    residual=np.linalg.norm(matrix@amplitudes-rhs)/np.linalg.norm(rhs)
    return dict(zip(labels,amplitudes)),float(residual)


def reciprocal_modal_directions(y_m,displacement_u,displacement_v,stress_xx,stress_xy,basis,labels,scale):
    """Bilinear Betti/traction projection at a complete uniform solid port.

    Use all four reference states and their numerical Gram matrix instead of
    assuming exact discrete orthogonality. No complex conjugation in Betti's form.
    """
    y=np.asarray(y_m);n=len(y)
    observed=np.array([displacement_u,displacement_v,np.asarray(stress_xx)/scale,np.asarray(stress_xy)/scale])
    if not np.isfinite(observed).all():raise ValueError('Reciprocity requires a complete pore-free solid section.')
    states=basis.T.reshape(len(labels),4,n)
    def betti(a,b):
        return np.trapezoid(a[0]*b[2]+a[1]*b[3]-a[2]*b[0]-a[3]*b[1],y)
    gram=np.array([[betti(test,state) for state in states] for test in states])
    if np.linalg.cond(gram)>1e8:raise ValueError('Reciprocal modal Gram matrix is ill-conditioned.')
    coefficients=np.linalg.solve(gram,np.array([betti(test,observed) for test in states]))
    return dict(zip(labels,coefficients))
