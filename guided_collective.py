"""Thin-plate modal validation and diameter/frequency collective-response study.

Uses example CFRP stiffness, Gaussian-field-derived irregular contours, explicit
closed cavities, uniform modal ports, and an incident-normalized 2-mode matrix.
No experimental detectability claim is inferred from pore diameter/wavelength.
"""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import argparse,json,hashlib,sys,time
from functools import lru_cache
import numpy as np
from shapely import affinity
from cfrp_pores.comsol import create_2d,add_frequency_study,add_modal_displacement,save_and_solve,interpolate
from cfrp_pores.dispersion import Plate,validate_material,modal_direction_basis,separate_modal_directions
from cfrp_pores.pores import load_polygons,reposition,polygon_statistics,write_csv

sys.stdout.reconfigure(encoding='utf-8')
MODES=['S0','A0'];LABELS=['S0_forward','S0_backward','A0_forward','A0_backward']
OUT=ROOT/'output/guided_collective'


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def settings():
    cfg=read(ROOT/'config.json');C,rho=validate_material(cfg)
    plate=Plate(.001,C[np.ix_([0,1,5],[0,1,5])],rho,120)
    # This reference plate is fixed for the study; reuse identical eigenproblems.
    plate.at_frequency=lru_cache(maxsize=32)(plate.at_frequency)
    wavelength=plate.at_frequency(50000.,'S0')['wavelength_m']
    start=2*wavelength;end=start+.016
    sponge=end+wavelength;length=sponge+4*wavelength
    # Keep the already-built geometric coordinates bit-identical across BLAS
    # implementations; material, mesh, and contours remain checked separately.
    layout_path=OUT/'study_settings.json'
    if layout_path.exists():
        layout=read(layout_path)
        start=layout['start'];end=start+.016
        sponge=layout['sponge'];length=layout['length']
    return dict(config=cfg,C=C,rho=rho,plate=plate,h=.001,start=start,
                sponge=sponge,length=length,absorb_length=length-sponge,
                pre_positions=[start-.012,start-.010,start-.008],
                post_positions=[end+.008,end+.010,end+.012],frequencies=[50000.,100000.,200000.])


def grouped_pores(diameter_um,length_mm,seed,porosity=.015):
    """Use original irregular contours, normalize diameter, then place without overlap.

    This is a controlled diameter-group derivative of the Gaussian field model;
    its final position correlation is NOT the original random-field correlation.
    All requested pores are retained. Skin/end buffers keep cavities internal.
    """
    library,_=load_polygons(ROOT/'output/sample/geometry.json')
    library=sorted(library,key=lambda p:p.area,reverse=True)[:64]
    rng=np.random.default_rng(seed)
    area=porosity*length_mm*1000*1000
    n=max(3,round(area/(np.pi*diameter_um**2/4*(1+.4**2/12))))
    factors=rng.uniform(.8,1.2,n);diameters=diameter_um*factors
    diameters*=np.sqrt(area/(np.pi*np.sum(diameters**2)/4))
    polygons=[]
    for i,d in enumerate(diameters):
        p=library[int(rng.integers(len(library)))]
        p=affinity.translate(p,-p.centroid.x,-p.centroid.y)
        factor=d/np.sqrt(4*p.area/np.pi)
        polygons.append(affinity.scale(p,xfact=factor,yfact=factor,origin=(0,0)))
    placed=reposition(polygons,(length_mm*1000-100,900),seed+10000)
    placed=[affinity.translate(p,50,50) for p in placed]
    assert np.isclose(sum(p.area for p in placed),area,rtol=1e-11)
    rows=polygon_statistics(placed)
    return placed,dict(nominal_diameter_um=diameter_um,porosity=porosity,length_mm=length_mm,seed=seed,
                       pore_count=n,actual_porosity=sum(p.area for p in placed)/(length_mm*1e6),
                       diameter_min_um=min(r['equivalent_diameter_um'] for r in rows),
                       diameter_mean_um=float(np.mean([r['equivalent_diameter_um'] for r in rows])),
                       diameter_max_um=max(r['equivalent_diameter_um'] for r in rows),
                       skin_um=50,shape_source='Gaussian-field contours; normalized size groups and independent uniform placement',
                       original_field_position_correlation_preserved=False)


def export_geometry(folder,pores,metadata):
    folder.mkdir(parents=True,exist_ok=True)
    payload=dict(unit='um',size_um=[metadata['length_mm']*1000,1000],polygons=[
        dict(exterior=np.asarray(p.exterior.coords).tolist(),interiors=[np.asarray(r.coords).tolist() for r in p.interiors]) for p in pores])
    (folder/'geometry.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
    (folder/'metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    write_csv(folder/'pores.csv',polygon_statistics(pores))


def dephase_average(coefficients,positions,k):
    return np.mean(np.asarray(coefficients)*np.exp(1j*k*np.asarray(positions)))


def modal_port_data(model,s,f,plate=None,samples=401):
    plate=plate or s['plate'];h=s['h']
    y=np.linspace(h*1e-5,h*(1-1e-5),samples)
    basis,labels,scale=modal_direction_basis(y-h/2,plate,f)
    positions=s['pre_positions']+s['post_positions'];coefficients=[];residuals=[]
    for x in positions:
        data=interpolate(model,[np.full(len(y),x),y],['u','v','solid.sx','solid.sxy'])
        coefficient,residual=separate_modal_directions(y-h/2,*[v[-1] for v in data],basis,labels,scale)
        coefficients.append([coefficient[label] for label in LABELS]);residuals.append(residual)
    coefficients=np.asarray(coefficients);waves={kind:plate.at_frequency(f,kind) for kind in MODES}
    result={}
    for port,index in [('pre',slice(0,3)),('post',slice(3,6))]:
        points=np.array(positions)[index]
        result[port+'_forward']=np.array([dephase_average(coefficients[index,2*i],points,waves[m]['k_rad_m']) for i,m in enumerate(MODES)])
        result[port+'_backward']=np.array([dephase_average(coefficients[index,2*i+1],points,-waves[m]['k_rad_m']) for i,m in enumerate(MODES)])
    # Modal power flux for unit forward displacement amplitude, per out-of-plane width.
    powers=[]
    for i,m in enumerate(MODES):
        col=basis[:,2*i].reshape(4,len(y));u,v=col[:2];sx,sy=col[2:]*scale
        flux=np.trapezoid(-.5*np.real(sx*np.conj(2j*np.pi*f*u)+sy*np.conj(2j*np.pi*f*v)),y)
        if flux<=0:raise RuntimeError('Forward modal power flux has the wrong sign.')
        powers.append(float(flux))
    result.update(coefficients=coefficients,residuals=residuals,powers=powers,
                  wavelength_m=[waves[m]['wavelength_m'] for m in MODES])
    return result


def set_source(model,mode,h):
    import jpype
    c=model.java.component('comp1')
    mode=dict(mode,u=np.conj(mode['u']),v=np.conj(mode['v']))
    if 'source' not in list(c.physics('solid').feature().tags()):
        add_modal_displacement(model,mode,h)
    else:
        for name,values in [('ur',mode['u'].real),('ui',mode['u'].imag),('vr',mode['v'].real),('vi',mode['v'].imag)]:
            c.func(name).set('table',jpype.JArray(jpype.JString,2)([[str(y+h/2),str(z)] for y,z in zip(mode['y_m'],values)]))


def matched_rectangle_area(client,s,repair_tolerance_m=None):
    """Subtract a kernel-matched rectangle, avoiding roundoff in tiny void areas."""
    key=(s['length'],s['h'],repair_tolerance_m)
    cache=getattr(client,'_guide_area_cache',{})
    if key in cache:return cache[key]
    model=client.create('area_reference')
    try:
        c=model.java.component().create('comp1',True);g=c.geom().create('geom1',2);g.lengthUnit('m')
        if repair_tolerance_m is not None:
            g.absRepairTol(repair_tolerance_m)
            g.feature('fin').set('repairtoltype','absolute');g.feature('fin').set('absrepairtol',repair_tolerance_m)
        rect=g.create('rect1','Rectangle');rect.set('size',[s['length'],s['h']]);g.run()
        measure=g.measureFinal();measure.selection().geom('geom1',2);measure.selection().all()
        area=float(measure.getArea())
    finally:client.remove(model)
    cache[key]=area;client._guide_area_cache=cache
    return area


def direct_pore_area(client,s,pore,repair_tolerance_m=None):
    """Measure the tiny object directly; large-minus-small area cancellation is unsafe."""
    import jpype
    from cfrp_pores.comsol import _polygon
    model=client.create('pore_area_reference')
    try:
        c=model.java.component().create('comp1',True);g=c.geom().create('geom1',2);g.lengthUnit('m')
        if repair_tolerance_m is not None:g.absRepairTol(repair_tolerance_m)
        tag=_polygon(g,'pore0',pore,(s['start'],0.));g.run(tag)
        measure=g.measure();measure.selection().init(jpype.JInt(2))
        measure.selection().set(tag,jpype.JArray(jpype.JInt)([1]))
        return float(measure.getArea())
    finally:client.remove(model)


def run_case(client,folder,s,pores=(),frequencies=None,mesh=.00015,keep_model=True,force=False,repair_tolerance_m=None,micro_divisor=30):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    frequencies=frequencies or s['frequencies']
    geometry_hash=hashlib.sha256(b''.join(p.wkb for p in pores)).hexdigest()
    signature=dict(geometry_hash=geometry_hash,frequencies=frequencies,mesh=mesh,
                   material=s['config'],layout=[s['start'],s['sponge'],s['length']],
                   solver_schema='uniform_ports_2source_matrix_v1')
    if repair_tolerance_m is not None:signature['absolute_repair_tolerance_m']=repair_tolerance_m
    if micro_divisor!=30:signature['micro_mesh_divisor']=micro_divisor
    default_post=[s['start']+.016+.008,s['start']+.016+.010,s['start']+.016+.012]
    if s['post_positions']!=default_post:signature['port_positions_m']=s['pre_positions']+s['post_positions']
    fingerprint=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
    cache=folder/'results.json';rawpath=folder/'signals.npz'
    if cache.exists() and not force:
        record=read(cache)
        if record.get('fingerprint')!=fingerprint:
            cached_signature=dict(record.get('signature',{}))
            tightened=cached_signature.pop('absolute_repair_tolerance_m',None)
            compatible=(repair_tolerance_m is None and tightened==1e-9 and
                        cached_signature==signature and record.get('finite_voids')==len(pores))
            if not compatible:raise ValueError('Cached case differs; use another folder: '+str(folder))
        return record
    minimum=min((np.sqrt(4*p.area/np.pi)*1e-6/micro_divisor for p in pores),default=mesh/30)
    minimum=min(mesh/30,minimum)
    print('BUILD',folder.as_posix(),'pores',len(pores),flush=True)
    model=create_2d(client,'guided_collective',s['config'],pores,size_m=(s['length'],s['h']),
                    pore_offset_m=(s['start'],0.),mesh_m=mesh,mesh_min_m=minimum,
                    guided=True,sponge_start_m=s['sponge'],repair_tolerance_m=repair_tolerance_m)
    geometry=model.java.component('comp1').geom('geom1')
    measure=geometry.measureFinal();measure.selection().geom('geom1',2);measure.selection().all()
    measured_void_area=matched_rectangle_area(client,s,repair_tolerance_m)-float(measure.getArea())
    requested_void_area=sum(p.area for p in pores)*1e-12
    measured_void_area_method='Matched rectangle minus final solid area'
    if len(pores)==1:
        measured_void_area=direct_pore_area(client,s,pores[0],repair_tolerance_m)
        measured_void_area_method='Direct tiny-pore object area; avoids subtractive cancellation'
    finite_voids=int(geometry.getNFiniteVoids())
    if finite_voids!=len(pores):
        client.remove(model)
        if repair_tolerance_m is None:
            print('RETRY with absolute CAD tolerance 1 nm',folder.as_posix(),flush=True)
            return run_case(client,folder,s,pores,frequencies,mesh,keep_model,force,
                            repair_tolerance_m=1e-9,micro_divisor=micro_divisor)
        raise RuntimeError('Final CAD pore count differs from input; tighten geometric repair tolerance.')
    if requested_void_area and abs(measured_void_area/requested_void_area-1)>.005:
        client.remove(model)
        raise RuntimeError('COMSOL repair changed pore area by >0.5%; use absolute tolerance in a fresh case folder.')
    model.java.component('comp1').physics('solid').feature('lemm1').feature('sponge').set(
        'alpha_dM',f'if(x>{s["sponge"]}[m],3*2*pi*f0*((x-{s["sponge"]}[m])/{s["absorb_length"]}[m])^3,0[1/s])')
    arrays={};rows=[]
    try:
        add_frequency_study(model,frequencies[0])
        for f in frequencies:
            model.java.param().set('f0',f'{f}[Hz]')
            launches=[]
            for kind in MODES:
                wave=s['plate'].at_frequency(f,kind);set_source(model,wave,s['h'])
                started=time.perf_counter();model.solve();elapsed=time.perf_counter()-started
                data=modal_port_data(model,s,f)
                launches.append(data)
                arrays[f'{int(f)}_{kind}_coefficients']=data['coefficients']
                arrays[f'{int(f)}_{kind}_residuals']=np.array(data['residuals'])
                print('SOLVED',folder.name,int(f),kind,'seconds',round(elapsed,2),'residual',max(data['residuals']),flush=True)
            incoming=np.column_stack([d['pre_forward'] for d in launches])
            outgoing=np.column_stack([d['post_forward'] for d in launches])
            back=np.column_stack([d['pre_backward'] for d in launches])
            right_in=np.column_stack([d['post_backward'] for d in launches])
            if np.linalg.cond(incoming)>100:raise RuntimeError('Independent modal excitations are ill-conditioned.')
            T=outgoing@np.linalg.inv(incoming);R=back@np.linalg.inv(incoming)
            right=right_in@np.linalg.inv(incoming)
            power=np.array(launches[0]['powers']);weight=power[:,None]/power[None,:]
            energy=(weight*(abs(T)**2+abs(R)**2)).sum(axis=0)
            arrays[f'{int(f)}_T']=T;arrays[f'{int(f)}_R']=R;arrays[f'{int(f)}_right_in']=right
            row=dict(frequency_hz=f,modes=MODES,transmission_real=T.real.tolist(),transmission_imag=T.imag.tolist(),
                     reflection_real=R.real.tolist(),reflection_imag=R.imag.tolist(),
                     right_incoming_ratio_max=float(np.max(abs(right))),
                     uniform_port_residual_max=max(max(d['residuals']) for d in launches),
                     energy_balance_by_incident_mode=energy.tolist(),modal_power_flux=power.tolist(),
                     wavelength_m=launches[0]['wavelength_m'])
            rows.append(row)
            np.savez_compressed(rawpath,positions_m=s['pre_positions']+s['post_positions'],**arrays)
        if keep_model:model.save(str((folder/'model.mph').resolve()))
        record=dict(fingerprint=fingerprint,signature=signature,results=rows,
                    dofs=len(model.java.sol('sol1').getU()),mesh_elements=int(model.java.component('comp1').mesh('mesh1').getNumElem()),
                    saved_model_state='Last frequency, A0 excitation; all modal probe arrays retained.',
                    requested_void_area_m2=requested_void_area,measured_void_area_m2=measured_void_area,
                    measured_void_area_method=measured_void_area_method,finite_voids=finite_voids,
                    material_status=s['config']['material_status'])
        cache.write_text(json.dumps(record,indent=2),encoding='utf-8')
        return record
    finally:client.remove(model)


def complex_T(record,f):
    row=next(r for r in record['results'] if r['frequency_hz']==f)
    return np.asarray(row['transmission_real'])+1j*np.asarray(row['transmission_imag'])


def contrast(record,reference,f):
    return np.diag(complex_T(record,f))/np.diag(complex_T(reference,f))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['baseline','sweep','single'])
    parser.add_argument('--samples',type=int,default=6);parser.add_argument('--cores',type=int,default=4);args=parser.parse_args()
    OUT.mkdir(exist_ok=True);s=settings();client=mph.Client(cores=args.cores,version='6.2')
    public={k:v for k,v in s.items() if k not in ['config','C','rho','plate']}
    public.update(diameter_groups_um=[20,40,80],porosity=.015,thresholds=[.0001,.001,.01],
                  primary_threshold=.001,modal_residual_tolerance=.01,
                  numerical_amplitude_error_target=.0002,
                  shape_model='Controlled groups of Gaussian-field-derived irregular contours; position field is regenerated independently.')
    setting_path=OUT/'study_settings.json'
    if setting_path.exists():
        if read(setting_path)!=public:raise ValueError('Study settings changed; use a fresh output directory.')
    else:setting_path.write_text(json.dumps(public,indent=2),encoding='utf-8')
    reference=run_case(client,OUT/'reference',s)
    if args.stage=='baseline':return
    if args.stage=='sweep':
        for seed in range(101,101+args.samples):
            for diameter in [20,40,80]:
                folder=OUT/f'd{diameter}/L8/seed_{seed}'
                pores,meta=grouped_pores(diameter,8,seed);export_geometry(folder/'geometry',pores,meta)
                run_case(client,folder,s,pores)
    else:
        for diameter in [20,40,80]:
            pores,meta=grouped_pores(diameter,8,101)
            ordered=sorted(pores,key=lambda p:p.area)
            for label,pore in [('minimum',ordered[0]),('median',ordered[len(ordered)//2]),('maximum',ordered[-1])]:
                folder=OUT/f'd{diameter}/single_{label}'
                export_geometry(folder/'geometry',[pore],dict(meta,single_case=label,actual_porosity=pore.area/8e6,pore_count=1))
                record=run_case(client,folder,s,[pore])
                values=[dict(frequency_hz=f,relative_transmission_change=np.abs(1-abs(contrast(record,reference,f))).tolist()) for f in s['frequencies']]
                (folder/'contrast.json').write_text(json.dumps(values,indent=2),encoding='utf-8')


if __name__=='__main__':main()
