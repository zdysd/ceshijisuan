"""Use the installed COMSOL Python runtime (ASCII path) on this Windows host.

JPype's native library cannot currently load reliably from the Chinese-named
workspace venv. Import it from the existing COMSOL runtime BEFORE adding the
workspace scientific dependencies. No server environment is modified.
"""
import mph  # imports JPype from the COMSOL runtime first
import site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))

import argparse
import json
import traceback
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from cfrp_pores.comsol import start_client,create_2d,create_3d,add_load,add_frequency_study,add_transient_study,save_and_solve,receiver_amplitudes,interpolate,add_modal_displacement,add_air_coupling,compact_transient
from cfrp_pores.pores import PoreSpec,export_sample,load_polygons,write_csv,simplify_geometry
from cfrp_pores.dispersion import Plate,validate_material
from cfrp_pores.attenuation import fit_attenuation,fit_quality


def bulk(args):
    cfg=json.loads(Path(args.config).read_text())
    root=Path(args.output); root.mkdir(parents=True,exist_ok=True)
    client=start_client()
    results=[]
    roi_size=tuple(cfg['pore_roi_um'])
    supplied=None
    if args.geometry:
        supplied,roi_size=load_polygons(args.geometry)
    length,width=np.array(roi_size)*1e-6
    receivers=.0003+np.linspace(.15,.95,5)*length
    reference=None
    for kind in ([args.case] if args.case else ['no_pore','single_pore','random_pores']):
        if kind=='single_pore':
            from shapely.geometry import Point
            polygons=[Point(roi_size[0]/2,roi_size[1]/2).buffer(20,quad_segs=24)]
        elif kind=='random_pores':
            if supplied is None:
                spec=PoreSpec(size_um=roi_size,pixel_um=cfg['pixel_um'],correlation_um=tuple(cfg['correlation_lengths_um']),porosity=args.porosity,seed=args.seed,angle_deg=cfg['orientation_deg'])
                export_sample(root/'random_geometry',spec)
                polygons,_=load_polygons(root/'random_geometry/geometry.json')
            else:polygons=supplied
        else: polygons=[]
        model=None
        try:
            print('BUILD',kind,flush=True)
            model=create_2d(client,'bulk_'+kind,cfg,polygons,size_m=(length+.001,width),pore_offset_m=(.0003,0),mesh_m=args.mesh)
            if args.transient:
                dt=args.dt_factor/(60*args.frequency)
                add_transient_study(model,args.frequency,2e-6,dt)
            else:
                add_load(model,['1e6[Pa]','0','0']); add_frequency_study(model,args.frequency)
            record=save_and_solve(model,root/(kind+'.mph'),save=not args.transient)
            record.update(case=kind,frequency_hz=args.frequency,mesh_max_m=args.mesh,status='solved')
            if not args.transient:
                C,rho=validate_material(cfg)
                a,missing=receiver_amplitudes(model,receivers,width,bulk_impedance=np.sqrt(rho*C[0,0]))
                if reference is None: reference=a
                record.update(fit_attenuation(receivers,a,reference))
                record.update(fit_quality(record))
                record['missing_receiver_fraction']=missing
                write_csv(root/(kind+'_receivers.csv'),[dict(x_m=float(x),real=float(z.real),imag=float(z.imag)) for x,z in zip(receivers,a)])
            else:
                # Actual output at centerline, for stability/arrival/step checks.
                data=interpolate(model,[receivers,np.full(len(receivers),width/2)])
                times=np.asarray(model.evaluate('t')).ravel()
                np.savez_compressed(root/(kind+'_transient.npz'),times_s=times,displacement=data.real,receivers_m=receivers)
                if args.keep_full_fields:model.save(str((root/(kind+'.mph')).resolve()))
                else:
                    record.update(compact_transient(model,root/(kind+'.mph'),times))
                    retained=interpolate(model,[receivers,np.full(len(receivers),width/2)])
                    if not np.allclose(retained[:,0,:].real,data[:,-1,:].real,rtol=1e-8,atol=1e-15):
                        raise RuntimeError('Final displacement changed during storage compaction.')
            results.append(record)
            print(json.dumps(record),flush=True)
        except Exception as e:
            results.append(dict(case=kind,status='failed',error=str(e))); traceback.print_exc()
            (root/'results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
            raise
        finally:
            if model is not None: client.remove(model)
        (root/'results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')


def guided(args):
    cfg=json.loads(Path(args.config).read_text());C,rho=validate_material(cfg)
    thickness=cfg['plate_thickness_m']
    plate=Plate(thickness,C[np.ix_([0,1,5],[0,1,5])],rho,60)
    root=Path(args.output);root.mkdir(parents=True,exist_ok=True)
    frequency=args.frequency
    if frequency>=plate.first_higher_cutoff_hz():
        raise ValueError('Selected frequency exceeds the first higher Lamb-mode cutoff; this two-mode demonstration needs a lower frequency.')
    client=start_client();records=[]
    # Correlation lengths are deliberately coarser than the bulk demonstration.
    # This is recorded; these two demo cases must not be compared as identical pores.
    if args.geometry:
        polygons,size=load_polygons(args.geometry)
        if not np.isclose(size[1]*1e-6,thickness):raise ValueError('Guided pore geometry height must equal the full plate thickness.')
    else:
        size=(40000,thickness*1e6)
        spec=PoreSpec(size,10,(100,50),args.porosity,args.seed)
        export_sample(root/'random_geometry',spec)
        polygons,_=load_polygons(root/'random_geometry/geometry.json')
    roi_length=size[0]*1e-6
    length=roi_length+.04;sponge_start=roi_length+.015
    receivers=.008+np.linspace(.1,.9,5)*roi_length
    polygons,simplification=simplify_geometry(polygons,1.)
    if abs(simplification['relative_void_area_change'])>.02:
        raise RuntimeError('Guided contour simplification changed void area by more than 2%.')
    (root/'simplification.json').write_text(json.dumps(simplification,indent=2))
    for kind in ['S0','A0']:
        mode=plate.at_frequency(frequency,kind)
        # SAFE uses exp(+ikx-iwt); COMSOL uses exp(+iwt).
        mode['u']=np.conj(mode['u']);mode['v']=np.conj(mode['v'])
        reference=None
        for label,pores in [('no_pore',[]),('random_pores',polygons)]:
            model=None
            try:
                print('BUILD',kind,label,flush=True)
                model=create_2d(client,f'{kind}_{label}',cfg,pores,size_m=(length,thickness),pore_offset_m=(.008,0.),
                        mesh_m=args.mesh if args.mesh!=4e-5 else .00015,guided=True,sponge_start_m=sponge_start)
                add_modal_displacement(model,mode,thickness)
                add_frequency_study(model,frequency)
                record=save_and_solve(model,root/f'{kind}_{label}.mph')
                a,missing=receiver_amplitudes(model,receivers,thickness,401,mode,plate=plate)
                if reference is None: reference=a
                record.update(case=label,mode=kind,frequency_hz=frequency,status='solved',
                              phase_velocity_m_s=mode['phase_velocity_m_s'],
                              missing_receiver_fraction=missing,**fit_attenuation(receivers,a,reference))
                record['nonmodal_residual_max']=max(model._modal_residuals)
                record['modal_decomposition_accepted']=record['nonmodal_residual_max']<=.1
                record.update(fit_quality(record))
                write_csv(root/f'{kind}_{label}_receivers.csv',[dict(x_m=float(x),real=float(z.real),imag=float(z.imag)) for x,z in zip(receivers,a)])
                records.append(record);print(json.dumps(record),flush=True)
                (root/'results.json').write_text(json.dumps(records,indent=2))
            finally:
                if model is not None:client.remove(model)


def three_dimensional_cad(args):
    cfg=json.loads(Path(args.config).read_text())
    root=Path(args.output);root.mkdir(parents=True,exist_ok=True)
    # Small method-verification block; width is deliberately a validation variable.
    spec=PoreSpec((600,300,300),10,(80,40,120),args.porosity,args.seed)
    export_sample(root/'random_geometry',spec)
    client=start_client();records=[];reference=None
    receivers=np.array([.00012,.00021,.00030,.00039,.00048])
    ys,zs=np.meshgrid(np.linspace(1e-7,.0003-1e-7,15),np.linspace(1e-7,.0003-1e-7,15))
    for label,path in [('no_pore',None),('random_pores',root/'random_geometry/pores.stl')]:
        model=None
        try:
            print('BUILD 3D',label,flush=True)
            model=create_3d(client,'three_'+label,cfg,path,size_m=(.0006,.0003,.0003),mesh_m=args.mesh)
            add_frequency_study(model,args.frequency)
            record=save_and_solve(model,root/(label+'.mph'))
            record['actual_geometry_porosity']=model._actual_porosity
            values=[]
            for x in receivers:
                data=interpolate(model,[np.full(ys.size,x),ys.ravel(),zs.ravel()],['u'])
                # A zero contribution in voids, normalized over total section.
                values.append(np.nan_to_num(data[0,-1],nan=0.).mean())
            values=np.array(values)
            if reference is None:reference=values
            record.update(case=label,status='solved',frequency_hz=args.frequency,**fit_attenuation(receivers,values,reference))
            records.append(record);print(json.dumps(record),flush=True)
            (root/'results.json').write_text(json.dumps(records,indent=2))
            write_csv(root/(label+'_receivers.csv'),[dict(x_m=x,real=z.real,imag=z.imag) for x,z in zip(receivers,values)])
        finally:
            if model is not None:client.remove(model)


def three_dimensional(args):
    from cfrp_pores.tetra import volume_mesh
    from run_mesh_model import solve_mesh_model
    cfg=json.loads(Path(args.config).read_text());root=Path(args.output);root.mkdir(parents=True,exist_ok=True)
    records=[]
    for case,phi in [('no_pore',0.),('random_pores',args.porosity)]:
        folder=root/(case+'_geometry')
        spec=PoreSpec((600,300,300),10,(80,40,120),phi,args.seed)
        export_sample(folder,spec);volume_mesh(folder,spec.size_um,args.mesh*1e6)
        record=solve_mesh_model(folder,root/(case+'.mph'),cfg,args.frequency)
        record['case']=case;records.append(record)
        (root/'results.json').write_text(json.dumps(records,indent=2));print(json.dumps(record),flush=True)


def air_comparison(args):
    from shapely.geometry import Point
    cfg=json.loads(Path(args.config).read_text());client=start_client()
    root=Path(args.output);root.mkdir(parents=True,exist_ok=True);records=[];signals=[]
    points=np.array([.00045,.00065,.00085,.00105,.00125])
    pore=Point(500,500).buffer(20,quad_segs=24)
    for label,air in [('cavity',False),('air',True)]:
        m=None
        try:
            print('BUILD',label,flush=True)
            m=create_2d(client,'air_check_'+label,cfg,[pore],pore_offset_m=(.0003,0),mesh_m=4e-5,air=air)
            if air:
                ids=list(m.java.component('comp1').selection('geom1_pore0_dom').entities(2))
                if len(ids)!=1:raise RuntimeError('Single-pore air-domain selection is not unique.')
                add_air_coupling(m,ids)
            add_load(m,['1e6[Pa]','0','0']);add_frequency_study(m,args.frequency)
            record=save_and_solve(m,root/(label+'.mph'))
            z,missing=receiver_amplitudes(m,points,.001)
            record.update(case=label,status='solved',frequency_hz=args.frequency,missing_receiver_fraction=missing)
            records.append(record);signals.append(z)
            (root/'results.json').write_text(json.dumps(records,indent=2))
        finally:
            if m is not None:client.remove(m)
    relative=np.abs(np.abs(signals[1])/np.abs(signals[0])-1)
    summary=dict(max_relative_amplitude_difference=float(relative.max()),
                 amplitude_difference_below_1_percent=bool(relative.max()<=.01),
                 note='Single 20 um radius pore at this frequency only; not a guarantee for every pore size or a cavity resonance.')
    (root/'comparison.json').write_text(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='output/comsol_bulk')
    parser.add_argument('--config',default='config.json')
    parser.add_argument('--frequency',type=float)
    parser.add_argument('--geometry',help='Generated 2D geometry.json; overrides pore generator for bulk/guided comparisons.')
    parser.add_argument('--porosity',type=float,default=.015)
    parser.add_argument('--seed',type=int,default=101)
    parser.add_argument('--case',choices=['no_pore','single_pore','random_pores'])
    parser.add_argument('--dt-factor',type=float,default=1.)
    parser.add_argument('--keep-full-fields',action='store_true',help='Store all time frames; can produce multi-GB files.')
    parser.add_argument('--mesh',type=float,default=4e-5)
    parser.add_argument('--transient',action='store_true')
    parser.add_argument('--guided',action='store_true')
    parser.add_argument('--three-d',action='store_true')
    parser.add_argument('--cad-three-d',action='store_true',help='Diagnostic CAD route; mesh-based route is preferred for STL pores.')
    parser.add_argument('--air',action='store_true')
    args=parser.parse_args()
    args.frequency=args.frequency or (200000. if args.guided else 5000000.)
    folder=Path(args.output);folder.mkdir(parents=True,exist_ok=True)
    status_file=folder/'status.json'
    status_file.write_text(json.dumps(dict(status='running',arguments=vars(args)),indent=2))
    try:
        if args.air:air_comparison(args)
        elif args.cad_three_d:three_dimensional_cad(args)
        elif args.three_d:three_dimensional(args)
        elif args.guided: guided(args)
        else: bulk(args)
    except Exception as error:
        status_file.write_text(json.dumps(dict(status='failed',arguments=vars(args),error=str(error),traceback=traceback.format_exc()),indent=2),encoding='utf-8')
        raise
    else:
        status_file.write_text(json.dumps(dict(status='solved_not_automatically_validated',arguments=vars(args)),indent=2))
