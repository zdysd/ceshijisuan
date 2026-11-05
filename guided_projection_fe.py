"""Betti modal coefficients using COMSOL element quadrature on a solid port."""
import numpy as np
from cfrp_pores.dispersion import modal_direction_basis
from cfrp_pores.comsol import integrate_section


def prepare_projection(model,plate,frequency):
    import jpype
    y=np.linspace(0,plate.thickness_m,801)
    basis,labels,scale=modal_direction_basis(y-plate.thickness_m/2,plate,frequency)
    states=basis.T.reshape(len(labels),4,len(y));c=model.java.component('comp1');expressions=[]
    def W(a,b):return np.trapezoid(a[0]*b[2]+a[1]*b[3]-a[2]*b[0]-a[3]*b[1],y)
    gram=np.array([[W(test,state) for state in states] for test in states])
    if np.linalg.cond(gram)>1e8:raise RuntimeError('FE reciprocal projection is ill-conditioned.')
    for i,state in enumerate(states):
        names=[]
        for j,profile in enumerate(state):
            prefix=f'gp{i}_{j}';names.append(f'(comp1.{prefix}r(y)+i*comp1.{prefix}i(y))')
            for suffix,values in [('r',profile.real),('i',profile.imag)]:
                tag=prefix+suffix
                if tag in list(c.func().tags()):c.func().remove(tag)
                function=c.func().create(tag,'Interpolation');function.set('funcname',tag)
                function.set('source','table');function.set('argunit','m');function.set('fununit','1')
                function.set('table',jpype.JArray(jpype.JString,2)([[str(x),str(v)] for x,v in zip(y,values)]))
        u,v,sx,sy=names
        expressions.append(f'{u}*solid.sx/({scale}[Pa/m])+{v}*solid.sxy/({scale}[Pa/m])-{sx}*u-{sy}*v')
    # Stored solution datasets use the compiled definitions from their solution;
    # synchronize newly added postprocessing functions without resolving fields.
    model.java.sol('sol1').updateSolution()
    return dict(gram=gram,labels=labels,expressions=expressions,h=plate.thickness_m)


def project_port(model,x,projection):
    rhs=integrate_section(model,x,projection['h'],projection['expressions'],quadrature_order=8)
    values=np.linalg.solve(projection['gram'],rhs)
    return dict(zip(projection['labels'],values))
