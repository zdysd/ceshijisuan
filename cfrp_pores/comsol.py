"""COMSOL 6.2 builder. No zero-stiffness or zero-density solid pore domains."""
import json
from pathlib import Path
import time

import numpy as np
from .dispersion import validate_material,project_mode


_client=None


def start_client():
    global _client
    import mph
    # Standalone avoids leaving a background mphserver process behind.
    if _client is None: _client=mph.Client(cores=2,version='6.2')
    return _client


def _box_selection(component,tag,dim,bounds,geometry='geom1'):
    from jpype import JInt
    sel=component.selection().create(tag,'Box')
    sel.set('entitydim',JInt(dim)); sel.set('condition','inside')
    for axis,(lo,hi) in zip('xyz',bounds):
        sel.set(axis+'min',str(lo)); sel.set(axis+'max',str(hi))
    return tag


def _polygon(g,tag,polygon,offset_m=(0.,0.),scale=1e-6):
    outer=g.create(tag,'Polygon')
    pts=np.asarray(polygon.exterior.coords)[:-1]*scale+offset_m
    outer.set('source','vectors'); outer.set('x',pts[:,0].tolist()); outer.set('y',pts[:,1].tolist())
    outer.set('type','solid')
    outer.set('selresult',True)
    if not polygon.interiors: return tag
    holes=[]
    for k,ring in enumerate(polygon.interiors):
        t=tag+'i'+str(k); f=g.create(t,'Polygon')
        pts=np.asarray(ring.coords)[:-1]*scale+offset_m
        f.set('source','vectors'); f.set('x',pts[:,0].tolist()); f.set('y',pts[:,1].tolist()); f.set('type','solid')
        holes.append(t)
    t=tag+'d'; f=g.create(t,'Difference')
    f.selection('input').set([tag]); f.selection('input2').set(holes)
    return t


def create_2d(client,name,config,polygons=(),size_m=(.002,.001),pore_offset_m=(0.,0.),mesh_m=4e-5,
              guided=False,sponge_start_m=None,air=False,lateral_periodic=False,mesh_min_m=None,repair_tolerance_m=None,
              source_x_m=None,displacement_order=None,curvature_factor=.35,periodic_mesh_match=False):
    started=time.perf_counter()
    C,rho=validate_material(config)
    m=client.create(name); j=m.java
    c=j.component().create('comp1',True); g=c.geom().create('geom1',2)
    g.lengthUnit('m')
    if repair_tolerance_m is not None:
        g.absRepairTol(repair_tolerance_m)
        g.feature('fin').set('repairtoltype','absolute')
        g.feature('fin').set('absrepairtol',repair_tolerance_m)
    rect=g.create('rect1','Rectangle'); rect.set('size',list(size_m));
    solid_tag='rect1'
    if source_x_m is not None:
        if not 0<source_x_m<size_m[0]:raise ValueError('Internal source must lie inside the specimen.')
        rect.set('size',[source_x_m,size_m[1]])
        right=g.create('rect2','Rectangle');right.set('pos',[source_x_m,0])
        right.set('size',[size_m[0]-source_x_m,size_m[1]])
        union=g.create('source_partition','Union');union.selection('input').set(['rect1','rect2'])
        union.set('intbnd',True);solid_tag='source_partition'
    pore_tags=[_polygon(g,'pore'+str(i),p,pore_offset_m) for i,p in enumerate(polygons)]
    if pore_tags:
        diff=g.create('subtract','Difference')
        diff.selection('input').set([solid_tag]); diff.selection('input2').set(pore_tags)
        if air:
            diff.set('keepsubtract',True)
    g.run()
    if g.getNDomains()<1:raise RuntimeError('Boolean geometry has no material domains; do not solve an empty model.')
    # Entity bounding boxes can be rounded at single precision. A fixed 1 nm
    # tolerance silently misses outer edges of longer specimens (e.g. 37.5 mm).
    padding=max(1e-9,2e-7*max(size_m))
    for tag,axis,value in [('left',0,0.),('right',0,size_m[0]),('bottom',1,0.),('top',1,size_m[1])]:
        bounds=[(-padding,size_m[0]+padding),(-padding,size_m[1]+padding)]
        bounds[axis]=(value-padding,value+padding)
        _box_selection(c,tag,1,bounds)
    if source_x_m is not None:
        _box_selection(c,'excitation',1,[(source_x_m-padding,source_x_m+padding),(-padding,size_m[1]+padding)])
    checked_selections={}
    if source_x_m is not None or displacement_order is not None:
        measure=g.measureFinal()
        expected_sides=[('left',size_m[1]),('right',size_m[1]),('top',size_m[0]),('bottom',size_m[0])]
        if source_x_m is not None:expected_sides.append(('excitation',size_m[1]))
        for tag,expected in expected_sides:
            entities=c.selection(tag).entities(1)
            if not len(entities):raise RuntimeError(f'Empty {tag} boundary selection.')
            measure.selection().geom('geom1',1);measure.selection().set(entities)
            actual=float(measure.getLength())
            if abs(actual/expected-1)>1e-6:
                raise RuntimeError(f'Incomplete {tag} boundary selection.')
            checked_selections[tag]=dict(entities=[int(n) for n in entities],length_m=actual,expected_length_m=expected)
    p=c.physics().create('solid','SolidMechanics','geom1')
    if displacement_order is not None:p.prop('ShapeProperty').set('order_displacement',str(displacement_order))
    le=p.feature('lemm1')
    le.set('SolidModel','Anisotropic'); le.set('AnisotropicOption','AnisotropicVo')
    le.set('DVo_mat','userdef')
    import jpype
    # COMSOL exposes its symmetric matrix property as a flattened array.
    le.set('DVo',jpype.JArray(jpype.JString)([f'{x:.16g}[Pa]' for x in C.ravel(order='F')]))
    le.set('rho_mat','userdef'); le.set('rho',str(rho))
    le.set('geometricNonlinearity','linear')
    for tag in (['right'] if guided else (['left','right'] if lateral_periodic else ['left','right','top','bottom'])):
        feature=p.create('absorb_'+tag,'LowReflectingBoundary',1)
        feature.selection().named(tag)
        if not guided:
            # COMSOL's default P/S speed reduction is not a directional
            # anisotropic impedance. Match principal-axis outgoing waves.
            normal=0 if tag in ['left','right'] else 1
            cl=np.sqrt(C[normal,normal]/rho)
            cs=np.sqrt(C[5,5]/rho)
            impedances=[rho*cs,rho*cs,rho*np.sqrt(C[4,4]/rho)]
            impedances[normal]=rho*cl
            # The automatically created sys1 is a BOUNDARY coordinate system;
            # its second axis is normal in 2D. d_i here is in global axes.
            feature.set('coordinateSystem','GlobalSystem')
            feature.set('DampingType','userdef')
            feature.set('d_i',[f'{z}[Pa*s/m]' for z in impedances])
    if lateral_periodic:
        ids=list(c.selection('top').entities(1))+list(c.selection('bottom').entities(1))
        per=p.create('periodic','PeriodicCondition',1)
        per.selection().set(jpype.JArray(jpype.JInt)(ids))
        if source_x_m is not None:
            # Internal source partitions split both sides into multiple edges.
            per.create('destination','DestinationDomains',1).selection().named('top')
    if sponge_start_m is not None:
        d=le.create('sponge','Damping',2)
        d.set('DampingType','RayleighDamping')
        d.set('alpha_dM',f'if(x>{sponge_start_m}[m],2e7[1/s]*((x-{sponge_start_m}[m])/{size_m[0]-sponge_start_m}[m])^3,0[1/s])')
        d.set('beta_dK','0[s]')
    mesh=c.mesh().create('mesh1')
    mesh.create('size1','Size')
    mesh.feature('size1').set('custom',True)
    mesh.feature('size1').set('hmax',str(mesh_m))
    mesh.feature('size1').set('hmin',str(mesh_m/30 if mesh_min_m is None else mesh_min_m))
    mesh.feature('size1').set('hgrad',1.35)
    mesh.feature('size1').set('hcurve',curvature_factor)
    if periodic_mesh_match:
        if not lateral_periodic:raise ValueError('Matching side meshes requires lateral periodicity.')
        edge=mesh.create('periodic_source_edges','Edge');edge.selection().named('bottom')
        measure=g.measureFinal()
        def ordered_edges(name):
            values=[]
            for entity in c.selection(name).entities(1):
                measure.selection().geom('geom1',1);measure.selection().set(jpype.JArray(jpype.JInt)([entity]))
                values.append((float(measure.getBoundingBox()[0]),int(entity)))
            return [entity for _,entity in sorted(values)]
        source_edges=ordered_edges('bottom');destination_edges=ordered_edges('top')
        if len(source_edges)!=len(destination_edges):raise RuntimeError('Periodic edge partitions differ.')
        for index,(source,destination) in enumerate(zip(source_edges,destination_edges)):
            copy=mesh.create(f'periodic_copy_{index}','CopyEdge')
            copy.selection('source').set(jpype.JArray(jpype.JInt)([source]))
            copy.selection('destination').set(jpype.JArray(jpype.JInt)([destination]))
    mesh.create('tri','FreeTri')
    mesh.run()
    if periodic_mesh_match:
        vertices=np.asarray(mesh.getVertex())
        bottom=np.sort(vertices[0,abs(vertices[1])<padding]);top=np.sort(vertices[0,abs(vertices[1]-size_m[1])<padding])
        if len(bottom)!=len(top) or not np.allclose(bottom,top,rtol=0,atol=1e-10):
            raise RuntimeError('The actual top/bottom mesh vertices do not match.')
        checked_selections['periodic_mesh']=dict(vertices_per_side=len(bottom),max_x_difference_m=float(np.max(abs(bottom-top))))
    m._build_seconds=time.perf_counter()-started
    m._boundary_selection_checks=checked_selections
    return m


def add_load(m,expressions):
    p=m.java.component('comp1').physics('solid')
    f=p.create('source','BoundaryLoad',1); f.selection().named('left')
    f.set('FperArea',expressions)


def add_frequency_study(m,frequency_hz):
    j=m.java; j.param().set('f0',f'{frequency_hz}[Hz]')
    j.study().create('std1'); j.study('std1').create('freq','Frequency')
    j.study('std1').feature('freq').set('plist','f0')
    j.study('std1').createAutoSequences('all')
    # The default fully coupled frequency solver is explicitly wired to direct.
    sol=j.sol('sol1'); step=sol.feature('s1')
    if 'd1' not in list(step.feature().tags()): step.create('d1','Direct')
    step.feature('d1').set('linsolver','pardiso')
    if 'fc1' in list(step.feature().tags()): step.feature('fc1').set('linsolver','d1')


def add_modal_displacement(m,mode,thickness_m,amplitude_m=1e-8):
    import jpype
    j=m.java;c=j.component('comp1')
    for name,values in [('ur',mode['u'].real),('ui',mode['u'].imag),('vr',mode['v'].real),('vi',mode['v'].imag)]:
        f=c.func().create(name,'Interpolation')
        f.set('source','table');f.set('funcname',name)
        table=[[str(y+thickness_m/2),str(z)] for y,z in zip(mode['y_m'],values)]
        f.set('table',jpype.JArray(jpype.JString,2)(table))
        f.set('argunit','m');f.set('fununit','1')
    p=c.physics('solid');f=p.create('source','Displacement1',1)
    f.selection().named('left');f.set('Direction',['prescribed','prescribed','free'])
    f.set('U0',[f'{amplitude_m}[m]*(ur(y)+i*ui(y))',f'{amplitude_m}[m]*(vr(y)+i*vi(y))','0'])


def add_transient_study(m,frequency_hz,end_time_s,dt_s):
    j=m.java
    j.param().set('f0',f'{frequency_hz}[Hz]'); j.param().set('dt',f'{dt_s}[s]'); j.param().set('tend',f'{end_time_s}[s]')
    pulse='if(t<5/f0,1e6[Pa]*sin(2*pi*f0*t)*(1-cos(2*pi*f0*t/5))/2,0[Pa])'
    add_load(m,[pulse,'0','0'])
    j.study().create('std1'); j.study('std1').create('time','Transient')
    j.study('std1').feature('time').set('tlist','range(0,dt,tend)')
    j.study('std1').createAutoSequences('all')
    step=j.sol('sol1').feature('t1')
    step.set('timemethod','genalpha'); step.set('tstepsgenalpha','manual')
    step.set('timestepgenalpha','dt'); step.set('rhoinf',1.)
    step.set('rtol',1e-5)
    if 'd1' not in list(step.feature().tags()): step.create('d1','Direct')
    step.feature('d1').set('linsolver','pardiso')
    if 'fc1' in list(step.feature().tags()): step.feature('fc1').set('linsolver','d1')


def interpolate(m,coords,expressions=('u','v')):
    """Evaluate complex FE solution at physical coordinates, retaining missing points."""
    import jpype
    j=m.java
    if 'sample' in list(j.result().numerical().tags()): j.result().numerical().remove('sample')
    n=j.result().numerical().create('sample','Interp')
    datasets=[str(t) for t in list(j.result().dataset().tags()) if str(j.result().dataset(t).getType())=='Solution' and str(j.result().dataset(t).getString('solution'))=='sol1']
    if len(datasets)!=1:raise RuntimeError('Interpolation needs a unique explicitly selected solution dataset.')
    n.set('data',datasets[0])
    n.set('expr',list(expressions)); n.set('coord',jpype.JArray(jpype.JDouble,2)(np.asarray(coords,float).tolist()))
    real=np.asarray(n.getData(),float)
    try: imaginary=np.asarray(n.getImagData(),float)
    except Exception:
        # Only real time-dependent solutions may legitimately lack imag data.
        if any(m.java.study(s).feature(t).getType()=='Frequency' for s in list(m.java.study().tags()) for t in list(m.java.study(s).feature().tags())):
            raise
        imaginary=np.zeros_like(real)
    return real+1j*imaginary


def receiver_amplitudes(m,x_m,width_m,samples=401,mode=None,bulk_impedance=None,plate=None):
    if mode is not None and plate is not None:
        from .dispersion import modal_direction_basis,separate_modal_directions
        y=np.linspace(width_m*1e-5,width_m*(1-1e-5),samples)
        frequency=m.java.param().evaluate('f0')
        basis,labels,scale=modal_direction_basis(y-width_m/2,plate,frequency)
        values=[];missing=[];residuals=[]
        for x in x_m:
            d=interpolate(m,[np.full(samples,x),y],['u','v','solid.sx','solid.sxy'])
            coefficient,residual=separate_modal_directions(y-width_m/2,*[v[-1] for v in d],basis,labels,scale)
            values.append(coefficient[mode['mode']+'_forward']);residuals.append(residual)
            missing.append(1-float(np.isfinite(d[:,-1,:]).all(axis=0).mean()))
        m._modal_residuals=residuals
        return np.array(values),missing
    if mode is None and bulk_impedance is not None:
        from .attenuation import forward_bulk_displacement
        freq=m.java.param().evaluate('f0')
        values=[];missing=[]
        for x in x_m:
            integrated=integrate_section(m,x,width_m,['u','solid.sx','1'])
            if integrated[2].real<=0:raise RuntimeError('Receiver section has no solid material.')
            mean_u=integrated[0]/integrated[2].real
            mean_stress=integrated[1]/width_m
            values.append(forward_bulk_displacement(mean_u,mean_stress,freq,bulk_impedance))
            missing.append(1-integrated[2].real/width_m)
        return np.array(values),missing
    y=np.linspace(width_m*1e-5,width_m*(1-1e-5),samples)
    values=[];missing=[]
    for x in x_m:
        data=interpolate(m,[np.full(samples,x),y],['u','v','solid.sx'] if bulk_impedance is not None else ['u','v'])
        # getData has axes expression, inner solution, point.
        u,v=data[0,-1],data[1,-1]
        finite=np.isfinite(u)&np.isfinite(v)
        missing.append(1-float(finite.mean()))
        if not finite.any(): raise RuntimeError('No solid points on receiver section.')
        u=np.where(finite,u,0.); v=np.where(finite,v,0.)
        # Intrinsic displacement of the solid phase: remove receiver occupancy
        # changes from missing void pixels, and report the missing fraction.
        if mode is None:
            mean_u=np.trapezoid(u,y)/np.trapezoid(finite.astype(float),y)
            if bulk_impedance is not None:
                from .attenuation import forward_bulk_displacement
                stress=np.where(finite,data[2,-1],0.)
                mean_stress=np.trapezoid(stress,y)/np.ptp(y)
                mean_u=forward_bulk_displacement(mean_u,mean_stress,m.java.param().evaluate('f0'),bulk_impedance)
            values.append(mean_u)
        else:
            centered=y-width_m/2
            pu=np.interp(centered,mode['y_m'],mode['u'].real)+1j*np.interp(centered,mode['y_m'],mode['u'].imag)
            pv=np.interp(centered,mode['y_m'],mode['v'].real)+1j*np.interp(centered,mode['y_m'],mode['v'].imag)
            norm=np.trapezoid(finite*(abs(pu)**2+abs(pv)**2),y)
            values.append(np.trapezoid(np.conj(pu)*u+np.conj(pv)*v,y)/norm)
    return np.asarray(values),missing


def integrate_section(m,x_m,width_m,expressions,quadrature_order=4):
    """FE quadrature on actual cut-line segments, including exact void gaps."""
    import jpype
    j=m.java
    if 'section' in list(j.result().dataset().tags()):j.result().dataset().remove('section')
    d=j.result().dataset().create('section','CutLine2D');d.set('data','dset1')
    d.set('genpoints',jpype.JArray(jpype.JDouble,2)([[x_m,0.],[x_m,width_m]]))
    if 'section_integral' in list(j.result().numerical().tags()):j.result().numerical().remove('section_integral')
    n=j.result().numerical().create('section_integral','IntLine');n.set('data','section');n.set('expr',expressions)
    n.set('intorder',str(quadrature_order))
    values=np.asarray(n.getReal()).ravel()+1j*np.asarray(n.getImag()).ravel()
    if len(values)!=len(expressions):raise RuntimeError('Section integration requires one frequency solution.')
    return values


def save_and_solve(m,path,save=True):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    m.solve()
    if len(m.java.sol('sol1').getU())==0:raise RuntimeError('Empty solution (zero DOFs), not a successful physical simulation.')
    elapsed=time.perf_counter()-started
    j=m.java
    if 'wave' not in list(j.result().tags()):
        dim=j.component('comp1').geom('geom1').getSDim()
        pg=j.result().create('wave','PlotGroup2D' if dim==2 else 'PlotGroup3D')
        pg.label('Displacement magnitude - illustrative material')
        surface=pg.create('surface','Surface')
        surface.set('expr','sqrt(abs(u)^2+abs(v)^2'+(')' if dim==2 else '+abs(w)^2)'))
    if save:m.save(str(path.resolve()))
    mesh=m.java.component('comp1').mesh('mesh1')
    return dict(solve_seconds=elapsed,build_seconds=getattr(m,'_build_seconds',None),
                dofs=len(m.java.sol('sol1').getU()),mesh_elements=int(mesh.getNumElem()),model=str(path.resolve()),
                comsol_version='6.2',material_status='illustrative_completion_not_measured_CFRP')


def compact_transient(m,path,times_s):
    """Keep the final FE state after complete receiver history has been exported."""
    import jpype
    sol=m.java.sol('sol1')
    last=jpype.JInt(len(times_s))
    state=sol.getU(last)
    try:derivative=sol.getUDot(last)
    except jpype.JException:derivative=None  # derivative storage can be disabled
    sol.clearSolutionData()
    sol.setU(state)
    if derivative is not None:sol.setUDot(derivative)
    sol.setPNames(['t']);sol.setPVals(jpype.JArray(jpype.JDouble)([float(times_s[-1])]))
    sol.createSolution()
    actual=np.asarray(m.evaluate('t')).ravel()
    if actual.size!=1 or not np.isclose(actual[0],times_s[-1]):
        raise RuntimeError('Compacted solution did not retain the correct final time.')
    m.save(str(Path(path).resolve()))
    return dict(saved_field_frames=1,receiver_history_frames=len(times_s),final_time_s=float(times_s[-1]),
                derivative_retained=derivative is not None)


def retain_transient_frames(m,path,times_s,count=96):
    """Retain selected actual FE states and derivatives, including both endpoints.

    COMSOL 6.2 solver getU/setU indices are ONE-based. Full receiver histories
    must be exported before calling this storage-only operation.
    """
    import jpype
    times=np.asarray(times_s,float)
    if count<2 or len(times)<2:raise ValueError('At least two frames required.')
    selected=np.unique(np.linspace(0,len(times)-1,min(count,len(times)),dtype=int))
    sol=m.java.sol('sol1')
    states=[sol.getU(jpype.JInt(int(i)+1)) for i in selected]
    try:derivatives=[sol.getUDot(jpype.JInt(int(i)+1)) for i in selected]
    except jpype.JException:derivatives=None
    sol.clearSolutionData()
    for i,state in enumerate(states):
        sol.setU(jpype.JInt(i+1),state)
        if derivatives is not None:sol.setUDot(jpype.JInt(i+1),derivatives[i])
    sol.setPNames(['t']);sol.setPVals(jpype.JArray(jpype.JDouble)(times[selected].tolist()));sol.createSolution()
    np.testing.assert_allclose(np.asarray(m.evaluate('t')).ravel(),times[selected],rtol=1e-12,atol=1e-18)
    m.save(str(Path(path).resolve()))
    return dict(saved_field_frames=len(selected),receiver_history_frames=len(times),
                final_time_s=float(times[-1]),derivative_retained=derivatives is not None,
                retained_time_indices=selected.tolist())


def create_3d(client,name,config,stl_path=None,size_m=(.0008,.0004,.0004),mesh_m=4e-5):
    import jpype
    C,rho=validate_material(config)
    model=client.create(name);j=model.java;c=j.component().create('comp1',True)
    g=c.geom().create('geom1',3);g.lengthUnit('m')
    b=g.create('block','Block');b.set('size',list(size_m))
    if stl_path:
        print('3D importing closed Gaussian surfaces',flush=True)
        f=g.create('voids','Import');f.set('type','mesh')
        # CAD route retained for inspection; mesh-based route avoids its slivers.
        f.set('simplifymesh',False);f.set('formsolid',True)
        f.set('meshfilename',str(Path(stl_path).resolve()))
        g.run('voids')
        # An STL containing disconnected pores produces voids(1), voids(2), ...
        # Selecting the feature tag itself is not a selection of those solids.
        objects=[str(t) for t in list(g.objectNames()) if str(t)!='block']
        model._imported_pore_objects=len(objects)
        if not objects:raise RuntimeError('STL import did not produce any pore solids.')
        d=g.create('subtract','Difference');d.selection('input').set(['block']);d.selection('input2').set(objects)
    print('3D building Boolean geometry',flush=True)
    g.run()
    if g.getNDomains()<1:raise RuntimeError('3D subtraction produced an empty specimen; check STL units and surface orientation.')
    measurement=g.measureFinal();measurement.selection().geom('geom1',3);measurement.selection().all()
    model._actual_porosity=1-float(measurement.getVolume())/np.prod(size_m)
    print('3D geometry complete; porosity',model._actual_porosity,flush=True)
    for tag,axis,value in [('left',0,0.),('right',0,size_m[0]),('ymin',1,0.),('ymax',1,size_m[1]),('zmin',2,0.),('zmax',2,size_m[2])]:
        bounds=[(-1e-9,s+1e-9) for s in size_m];bounds[axis]=(value-1e-9,value+1e-9)
        _box_selection(c,tag,2,bounds)
    p=c.physics().create('solid','SolidMechanics','geom1');le=p.feature('lemm1')
    le.set('SolidModel','Anisotropic');le.set('AnisotropicOption','AnisotropicVo')
    le.set('DVo_mat','userdef');le.set('DVo',jpype.JArray(jpype.JString)([str(x)+'[Pa]' for x in C.ravel(order='F')]))
    le.set('rho_mat','userdef');le.set('rho',str(rho));le.set('geometricNonlinearity','linear')
    for tag in ['right','ymin','ymax','zmin','zmax']:
        f=p.create('abs_'+tag,'LowReflectingBoundary',2);f.selection().named(tag)
        normal=0 if tag=='right' else (1 if tag.startswith('y') else 2)
        shear_map={(0,1):5,(0,2):4,(1,2):3}
        impedance=[]
        for axis in range(3):
            stiffness=C[normal,normal] if axis==normal else C[shear_map[tuple(sorted((axis,normal)))],shear_map[tuple(sorted((axis,normal)))]]
            impedance.append(f'{np.sqrt(rho*stiffness)}[Pa*s/m]')
        f.set('coordinateSystem','GlobalSystem');f.set('DampingType','userdef');f.set('d_i',impedance)
    f=p.create('source','BoundaryLoad',2);f.selection().named('left');f.set('FperArea',['1e6[Pa]','0','0'])
    mesh=c.mesh().create('mesh1');s=mesh.create('size1','Size');s.set('custom',True)
    s.set('hmax',str(mesh_m));s.set('hmin',str(mesh_m/100));s.set('hgrad',1.35)
    if stl_path:model.save(str((Path(stl_path).resolve().parent/'comsol_geometry.mph')))
    mesh.create('tet','FreeTet');print('3D meshing',flush=True);mesh.run();print('3D mesh elements',mesh.getNumElem(),flush=True)
    return model


def add_air_coupling(m,air_domains):
    """Explicit fluid domains, never an elastic solid with G=0."""
    import jpype
    c=m.java.component('comp1');solid=c.physics('solid')
    domains=list(range(1,c.geom('geom1').getNDomains()+1))
    solid.selection().set(jpype.JArray(jpype.JInt)([d for d in domains if d not in air_domains]))
    ac=c.physics().create('acpr','PressureAcoustics','geom1')
    ac.selection().set(jpype.JArray(jpype.JInt)(air_domains))
    candidates=[ac.feature(t) for t in list(ac.feature().tags()) if 'rho_mat' in list(ac.feature(t).properties()) and 'c_mat' in list(ac.feature(t).properties())]
    if len(candidates)!=1:
        raise RuntimeError('Cannot uniquely identify pressure-acoustics fluid model: '+str([(str(t),str(ac.feature(t).getType())) for t in list(ac.feature().tags())]))
    fluid=candidates[0];fluid.set('rho_mat','userdef');fluid.set('rho','1.24[kg/m^3]')
    fluid.set('c_mat','userdef');fluid.set('c','344[m/s]')
    coupling=c.multiphysics().create('asb','AcousticStructureBoundary',1);coupling.selection().all()
    size=c.mesh('mesh1').create('air_size','Size')
    size.selection().geom('geom1',2);size.selection().set(jpype.JArray(jpype.JInt)(air_domains))
    size.set('custom',True);size.set('hmax','8e-6');size.set('hmin','1e-6')
    # A local Size must precede FreeTri in the sequence.
    c.mesh('mesh1').feature().move('air_size',1)
    c.mesh('mesh1').run()
