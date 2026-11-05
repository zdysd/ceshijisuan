"""Independent, bounded 200 kHz / 48 mm detectability verification.

P0 and geometry preparation are offline. Native solves use one serial worker;
its runtime resource policy may permit overlap without changing the old queue.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,csv,hashlib,json,re,sys
from datetime import datetime,timezone
from pathlib import Path
import site
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import numpy as np
from scipy.stats import beta
from scipy.spatial.distance import pdist
from shapely import affinity
from shapely.geometry import Polygon,box
from cfrp_pores.pores import load_polygons as _load_polygons
from cfrp_pores.portability import project_path

def load_polygons(path):return _load_polygons(project_path(path))
from cfrp_pores.transient import source_band

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'output/detectability_200kHz_L48'
OLD=ROOT/'output/guided_collective'
MODES=['S0','A0'];TAU=.001;F0=200000.

def read(p):return json.loads(project_path(p).read_text(encoding='utf-8-sig'))
def digest(p):return hashlib.sha256(project_path(p).read_bytes()).hexdigest()
def signature(v):return hashlib.sha256(json.dumps(v,sort_keys=True,ensure_ascii=True).encode()).hexdigest()
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()),encoding='utf-8');q.replace(p)
def stamp():return datetime.now(timezone.utc).isoformat()
def freeze(p,v):
    if Path(p).exists():
        if read(p)!=v:raise RuntimeError(f'Frozen configuration differs: {p}')
    else:write(p,v)
def geometry_hash(polys):return hashlib.sha256(b''.join(p.wkb for p in polys)).hexdigest()
def geometry_payload(polys,length=48.):
    return dict(unit='um',size_um=[length*1000,1000],polygons=[dict(exterior=list(p.exterior.coords),interiors=[list(r.coords) for r in p.interiors]) for p in polys])
def geometry_check(polys,group=True):
    inside=box(0,0,48000,1000)
    overlaps=[(i,j) for i,a in enumerate(polys) for j,b in enumerate(polys[:i]) if a.intersects(b)]
    phi=sum(p.area for p in polys)/48e6
    good=all(p.is_valid and inside.contains(p) for p in polys) and not overlaps
    if group:good=good and abs(phi-.015)<=read(ROOT/'config.json')['geometry_porosity_absolute_tolerance']
    return dict(passed=bool(good),pore_count=len(polys),area_fraction=phi,overlaps=overlaps,all_strictly_inside=all(inside.contains(p) for p in polys))
def cp_lower(k,n):return 0. if k==0 else float(beta.ppf(.05,k,n-k+1))
def classify(loss,u,*,single=False):
    if u is None:return 'unverified_uncertainty'
    if not np.isfinite(loss) or not np.isfinite(u) or u<0:raise ValueError('Finite loss and nonnegative uncertainty required')
    value=abs(loss) if single else loss
    if value-u>TAU:return 'above'
    if value+u<TAU:return 'below'
    return 'indeterminate'

def dev_path(seed):return OLD/(f'd80/L48/seed_{seed}' if seed<=106 else f'extended_samples/d80/L48/seed_{seed}')
def source_record(folder):
    folder=project_path(folder);record=read(folder/'results.json');sig=record['signature']
    if signature(sig)!=record['fingerprint']:raise RuntimeError(f'Bad fingerprint: {folder}')
    with np.load(folder/'signals.npz') as z:
        positions=z['positions_m'].copy();row=next(v for v in record['results'] if v['frequency_hz']==F0)
        k=2*np.pi/np.asarray(row['wavelength_m']);inc=[];out=[];back=[];right=[];res=[]
        for mode in MODES:
            c=z[f'200000_{mode}_coefficients'];res.extend(z[f'200000_{mode}_residuals'].tolist())
            def avg(x,sl,sign):return np.mean(x[sl]*np.exp(sign*1j*positions[sl,None]*k),axis=0)
            inc.append(avg(c[:,[0,2]],slice(0,3),1));out.append(avg(c[:,[0,2]],slice(3,6),1))
            back.append(avg(c[:,[1,3]],slice(0,3),-1));right.append(avg(c[:,[1,3]],slice(3,6),-1))
        incoming=np.column_stack(inc);cond=float(np.linalg.cond(incoming))
        matrices={key:np.linalg.solve(incoming.T,np.column_stack(a).T).T for key,a in [('T',out),('R',back),('right',right)]}
        for key,raw in [('T','200000_T'),('R','200000_R'),('right','200000_right_in')]:
            np.testing.assert_allclose(matrices[key],z[raw],rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(matrices['T'],np.array(row['transmission_real'])+1j*np.array(row['transmission_imag']),rtol=1e-10,atol=1e-12)
    powers=np.asarray(row['modal_power_flux']);energy=np.sum(powers[:,None]/powers[None,:]*(abs(matrices['T'])**2+abs(matrices['R'])**2),axis=0)
    np.testing.assert_allclose(energy,row['energy_balance_by_incident_mode'],rtol=1e-10,atol=1e-12)
    declared=sig.get('port_positions_m',read(OLD/'study_settings.json')['pre_positions']+read(OLD/'study_settings.json')['post_positions'])
    np.testing.assert_allclose(positions,declared,rtol=0,atol=1e-12)
    return dict(folder=folder,record=record,signature=sig,positions=positions,**matrices,powers=powers,
                residual=max(res),right_max=float(abs(matrices['right']).max()),condition=cond,energy=energy)

def compatible(a,b):
    return a['signature']['material']==b['signature']['material'] and a['signature']['layout']==b['signature']['layout'] and \
        a['signature']['solver_schema']==b['signature']['solver_schema'] and np.allclose(a['positions'],b['positions'],rtol=0,atol=1e-12)

def metrics(a,ref,case_id,category,seed=None,u=None,coverage='missing'):
    if not compatible(a,ref):raise ValueError('Incompatible geometry layout/material/ports for normalization')
    t=np.diag(a['T'])/np.diag(ref['T']);rows=[]
    for i,mode in enumerate(MODES):
        loss=float(1-abs(t[i]));score=abs(loss);uncertainty=None if u is None else float(u[i]);j=1-i
        rows.append(dict(case_id=case_id,category=category,seed=seed,mode=mode,frequency_hz=F0,
            T_real=float(t[i].real),T_imag=float(t[i].imag),loss=loss,score=score,phase_deg=float(np.angle(t[i])*180/np.pi),
            reflection_S0_real=float(a['R'][0,i].real),reflection_S0_imag=float(a['R'][0,i].imag),
            reflection_A0_real=float(a['R'][1,i].real),reflection_A0_imag=float(a['R'][1,i].imag),
            reflected_change_max=float(abs(a['R'][:,i]-ref['R'][:,i]).max()),
            conversion_real=float(a['T'][j,i].real),conversion_imag=float(a['T'][j,i].imag),
            conversion_change=float(abs(a['T'][j,i]-ref['T'][j,i])/abs(ref['T'][i,i])),
            incident_condition=a['condition'],fit_residual=a['residual'],right_incoming=a['right_max'],
            power_balance=float(a['energy'][i]),port_fit_passed=a['residual']<=.01,
            reference_right_passed=ref['right_max']<=.001,power_balance_passed=abs(a['energy'][i]-1)<=.002,
            u_num=uncertainty,uncertainty_coverage=coverage,classification=classify(loss,uncertainty,single=category in ['matched_single','legacy_single']),
            raw_threshold_exceeded=score>TAU,conservative_above=((score if category in ['matched_single','legacy_single'] else loss)-uncertainty>TAU) if uncertainty is not None else None,
            source=str(a['folder'].resolve()),reference_source=str(ref['folder'].resolve()),
            source_fingerprint=a['record']['fingerprint'],reference_fingerprint=ref['record']['fingerprint']))
    return rows

def manifest_case(folder,status,reason,geometry_source=None):
    folder=Path(folder);paths=[folder/n for n in ['results.json','signals.npz'] if (folder/n).exists()]
    for p in [folder/'geometry/geometry.json',folder/'geometry/metadata.json']:
        if p.exists():paths.append(p)
    record=read(folder/'results.json') if (folder/'results.json').exists() else {}
    return dict(path=str(folder.resolve()),reuse=status,reason=reason,
        fingerprint=record.get('fingerprint'),signature=record.get('signature'),
        saved_model_state=record.get('saved_model_state'),geometry_source=geometry_source,
        files=[dict(path=str(p.resolve()),sha256=digest(p),bytes=p.stat().st_size) for p in paths],
        model_mph=dict(exists=(folder/'model.mph').exists(),use='Last-frequency A0 only; not a source of missing S0 fields'))

def audit():
    OUT.mkdir(parents=True,exist_ok=True)
    refs={'L48':source_record(OLD/'reference_L48'),'legacy':source_record(OLD/'reference'),
          'fine48':source_record(OLD/'checks/reference_L48_fine')}
    manifest=[manifest_case(OLD/'reference_L48','可直接复用数据','200 kHz matched 48 mm six ports'),
        manifest_case(OLD/'checks/reference_L48_fine','可直接复用数据','Matched fine-grid 48 mm reference'),
        manifest_case(OLD/'reference','可直接复用数据','Legacy single-pore normalization only; not interchangeable with L48 ports')]
    rows=[];dev=[];mesh=[]
    audit_map={v['case'].replace('\\','/'):v for v in read(OLD/'geometry_audit.json')}
    for seed in range(101,113):
        folder=dev_path(seed);a=source_record(folder)
        polygons,_=load_polygons(folder/'geometry/geometry.json');assert geometry_hash(polygons)==a['signature']['geometry_hash']
        geom=geometry_check(polygons);assert geom['passed']
        native_count=a['record'].get('finite_voids');historical=audit_map.get(folder.relative_to(OLD).as_posix())
        if native_count is None and historical and historical['passed']:native_count=historical['cad_cavities']
        if native_count!=len(polygons):raise RuntimeError(f'Missing or inconsistent CAD count for {folder}')
        native_phi=a['record']['measured_void_area_m2']/(.048*.001)
        if abs(native_phi-.015)>1e-4:raise RuntimeError('Native area fraction failed')
        local=metrics(a,refs['L48'],f'development_seed_{seed}','development_group',seed)
        rows+=local;dev.append(dict(seed=seed,T=[complex(v['T_real'],v['T_imag']) for v in local],min_loss=min(v['loss'] for v in local)))
        entry=manifest_case(folder,'可直接复用数据','200 kHz raw dual-launch coefficients reconstructed; numerical validation still incomplete')
        entry.update(geometry=geom,native_pore_count=native_count,native_area_fraction=native_phi,native_count_source='results.json' if 'finite_voids' in a['record'] else 'geometry_audit.json')
        manifest.append(entry)
    fine=source_record(OLD/'checks/d80_L48_fine');coarse=source_record(dev_path(101))
    assert fine['signature']['geometry_hash']==coarse['signature']['geometry_hash']
    a=np.diag(coarse['T'])/np.diag(refs['L48']['T']);b=np.diag(fine['T'])/np.diag(refs['fine48']['T'])
    for i,m in enumerate(MODES):
        delta=float(abs(abs(a[i])-abs(b[i])));relative=delta/max(abs(1-abs(a[i])),abs(1-abs(b[i])),1e-12);phase=float(abs(np.angle(a[i]/b[i]))*180/np.pi)
        mesh.append(dict(seed=101,mode=m,magnitude_delta=delta,relative_loss_change=relative,phase_error_deg=phase,
                         passed=bool(delta<=.0002 and relative<=.05 and phase<=1),scope='This seed and this mesh pair only; not total u_num'))
        next(v for v in rows if v['seed']==101 and v['mode']==m)['measured_mesh_difference']=delta
    manifest.append(manifest_case(OLD/'checks/d80_L48_fine','可直接复用数据','Same geometry and ports as seed101; paired mesh sensitivity only',str(dev_path(101)/'geometry/geometry.json')))
    singles=[]
    for folder in sorted((OLD/'d80').iterdir()):
        if not folder.is_dir() or not (folder.name.startswith('single_') or folder.name.startswith('global_maximum_')):continue
        a=source_record(folder);local=metrics(a,refs['legacy'],folder.name,'legacy_single')
        rows+=local;singles+=local
        manifest.append(manifest_case(folder,'可复用几何或设置但须补算','Legacy ports and longitudinal pore location differ from required L48 controls; saved MPH only A0'))
    strongest=max(singles,key=lambda x:(x['score'],x['case_id'],x['mode']))
    order=sorted(dev,key=lambda x:(x['min_loss'],x['seed']))
    reps=sorted(set([101,order[0]['seed'],order[(len(order)-1)//2]['seed'],order[-1]['seed']]))
    stats=[]
    for i,mode in enumerate(MODES):
        values=np.array([v['T'][i] for v in dev]);loss=1-abs(values);k=int(np.sum(abs(loss)>TAU))
        stats.append(dict(mode=mode,n=12,raw_score_successes=k,raw_loss_successes=int(np.sum(loss>TAU)),
            development_only=True,one_sided95_lower_descriptive=cp_lower(k,12),
            coherent_amplitude_loss=float(1-abs(values.mean())),mean_individual_loss=float(loss.mean()),
            strongest_legacy_single_score=max(v['score'] for v in singles if v['mode']==mode),
            strongest_legacy_single=max((v for v in singles if v['mode']==mode),key=lambda v:v['score'])['case_id'],
            sensitivities=[dict(tau=t,raw_score_successes=int(np.sum(abs(loss)>t))) for t in [.0001,.001,.01]],
            independent_validation_successes=None))
    support=[]
    for name in ['验证结果_中文.md','PROTOCOL.md','study_settings.json','extended_attenuation.json','single_controls_summary.json',
                 'numerical_checks.json','geometry_audit.json','fe_port_check.json','projection_crosscheck.json','basis_convergence.json']:
        p=OLD/name;item=dict(path=str(p.resolve()),sha256=digest(p),bytes=p.stat().st_size)
        if name in ['fe_port_check.json','projection_crosscheck.json']:
            data=read(p);item.update(reuse='method/background only; not P1 L48 dual-incident error budget',
                cases=sorted(set(v['case'] for v in data)),observable='Existing A0/direct reciprocity checks; not a full two-launch L48 matched control')
        elif name=='basis_convergence.json':item.update(reuse='200 kHz SAFE basis convergence; not complete 17-point new-band validation')
        elif name in ['extended_attenuation.json','numerical_checks.json']:item.update(reuse='Background and explicitly matched center-frequency controls; alpha stability is not the detectability gate')
        support.append(item)
    tv=ROOT/'output/transient_verified';guide=read(tv/'validation.json')
    manifest.append(dict(path=str(tv),reuse='可复用几何或设置但须补算',reason='8 mm transient validation reuses implementation/settings only; different geometry/ports/layout, not evidence for 48 mm',
        guided_checks=[dict(case=c['case'],complete=c['complete'],passed=c['passed']) for c in guide['checks'] if c['kind']=='guided'],
        validation_sha256=digest(tv/'validation.json')))
    for name in ['absorber_L48','matched_single_controls_L48','frequency_band_L48','transient_L48','independent_30']:
        manifest.append(dict(path=str(OUT/name),reuse='缺失',reason='Required by the frozen new validation protocol; no accepted compatible evidence found'))
    summary=dict(updated_at=stamp(),stage='P0',passed=True,development_statistics=stats,mesh_checks=mesh,
        representative_seeds=reps,representative_order=[v['seed'] for v in order],strongest_legacy_single=strongest['case_id'],
        raw_matrices_reconstructed_from_two_launches=True,legacy_single_not_L48_control=True,
        source_protocol=str(ROOT/'200kHz_S0_A0_48mm_验证执行提示词.md'))
    write(OUT/'p0_audit.json',summary);write(OUT/'reuse_manifest.json',dict(cases=manifest,supporting_evidence=support))
    write(OUT/'per_case_metrics.json',rows);write_csv(rows)
    return summary

def write_csv(rows):
    fields=sorted(set().union(*(r.keys() for r in rows)))
    with (OUT/'per_case_metrics.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def prepare():
    import guided_collective as gc
    from guided_length_study import length_settings
    audit_result=read(OUT/'p0_audit.json')
    protocol_path=OUT/'protocol.json'
    if protocol_path.exists():return read(protocol_path)
    # Inspect filenames/geometry metadata only; never inspect a candidate's T.
    seen=set();evidence=[]
    for p in (ROOT/'output').rglob('results.json'):
        if p.is_relative_to(OUT):continue
        match=re.search(r'seed_(\d+)',str(p))
        if match:seen.add(int(match[1]));evidence.append(str(p))
        else:
            metadata=p.parent/'geometry/metadata.json'
            if metadata.exists():
                seed=read(metadata).get('seed')
                if seed is not None:seen.add(int(seed))
    seeds=[];candidate=113
    while len(seeds)<30:
        if candidate not in seen:seeds.append(candidate)
        candidate+=1
    seed_manifest=dict(seeds=seeds,n=30,excluded_preexisting_response_seeds=sorted(s for s in seen if s>=113),
                       screening='Response filenames and geometry metadata only; no candidate T read',screened_paths_sha256=signature(sorted(evidence)))
    write(OUT/'independent_seeds.json',seed_manifest)
    s=length_settings(gc.settings(),48)
    layout={k:s[k] for k in ['h','start','sponge','length','absorb_length','pre_positions','post_positions']}
    pool=[];geometries=[]
    for seed in seeds:
        polys,meta=gc.grouped_pores(80.,48.,seed,.015);check=geometry_check(polys)
        if not check['passed']:raise RuntimeError('Frozen geometry invalid')
        path=OUT/f'frozen_geometries/seed_{seed}/geometry.json';write(path,geometry_payload(polys));write(path.with_name('metadata.json'),meta)
        geometries.append(dict(seed=seed,path=str(path),geometry_hash=geometry_hash(polys),file_sha256=digest(path),audit=check))
        pool.extend((p,dict(seed=seed,index=i,source=str(path),split='independent_geometry_only')) for i,p in enumerate(polys))
    for seed in range(101,113):
        for length in [8,24,48]:
            path=OLD/(f'd80/L{length}/seed_{seed}/geometry/geometry.json' if seed<=106 else f'extended_samples/d80/L{length}/seed_{seed}/geometry/geometry.json')
            if not path.exists():continue
            polys,_=load_polygons(path)
            pool.extend((p,dict(seed=seed,length_mm=length,index=i,source=str(path),split='development')) for i,p in enumerate(polys))
    largest=max(pool,key=lambda v:v[0].area)
    def span(p):return float(pdist(np.asarray(p.convex_hull.exterior.coords)[:-1]).max())
    widest=max(pool,key=lambda v:span(v[0]))
    old_path=OLD/'d80'/audit_result['strongest_legacy_single']/'geometry/geometry.json';old,_=load_polygons(old_path)
    chosen=[];seen_shapes=set()
    for criterion,(p,origin) in [('max_area',largest),('max_span',widest),('strongest_development_response',(old[0],dict(source=str(old_path),split='development_single')))]:
        centered=affinity.translate(p,-p.centroid.x,-p.centroid.y)
        # Round only for translation-invariant de-duplication; saved geometry stays unrounded.
        key=signature(dict(exterior=np.round(centered.exterior.coords,8).tolist(),interiors=[np.round(x.coords,8).tolist() for x in centered.interiors]))
        if key in seen_shapes:continue
        seen_shapes.add(key);chosen.append((criterion,centered,origin))
    controls=[]
    for shape,p,origin in chosen:
        for rotation in [0,90]:
            turned=affinity.rotate(p,rotation,origin=(0,0))
            for depth in ['center','upper','lower']:
                dy=500-turned.centroid.y if depth=='center' else 950-turned.bounds[3] if depth=='upper' else 50-turned.bounds[1]
                placed=affinity.translate(turned,24000-turned.centroid.x,dy)
                case_id=f'single_{shape}_{depth}_r{rotation}';check=geometry_check([placed],False)
                valid=check['passed'] and abs(placed.area/p.area-1)<1e-12
                path=OUT/f'single_controls/{case_id}/geometry.json'
                if valid:write(path,geometry_payload([placed]))
                controls.append(dict(case_id=case_id,criterion=shape,origin=origin,rotation_deg=rotation,depth=depth,
                    usable=bool(valid),reason=None if valid else 'Not fully embedded; no clipping permitted',area_um2=p.area,
                    equivalent_diameter_um=float(np.sqrt(4*p.area/np.pi)),span_um=span(p),geometry_hash=geometry_hash([placed]),geometry_path=str(path)))
    write(OUT/'single_controls.json',controls);write(OUT/'frozen_geometry_manifest.json',geometries)
    models=[dict(case_id='reference',kind='reference',geometry_path=None,geometry_hash=geometry_hash([]))]
    models += [dict(case_id=f'group_seed_{seed}',kind='group',seed=seed,geometry_path=str(dev_path(seed)/'geometry/geometry.json'),
                    geometry_hash=read(dev_path(seed)/'results.json')['signature']['geometry_hash']) for seed in audit_result['representative_seeds']]
    models += [dict(case_id=v['case_id'],kind='single',geometry_path=v['geometry_path'],geometry_hash=v['geometry_hash']) for v in controls if v['usable']]
    write(OUT/'p1_models.json',models)
    protocol=dict(schema='detectability_200kHz_L48_v1',frequency_hz=F0,modes=MODES,thickness_m=.001,length_mm=48.,porosity=.015,nominal_diameter_um=80.,
        tau=TAU,sensitivity_thresholds=[.0001,.01],score='abs(1-abs(T))',loss='1-abs(T)',
        normalization='Two independent launches; outgoing @ inverse(actual incoming); same-mode diagonal divided by matched reference diagonal',
        development_seeds=list(range(101,113)),independent_seeds=seeds,n_validation=30,
        exact_binomial=dict(method='Clopper-Pearson one-sided 95% lower',target=.9,lower_30of30=cp_lower(30,30),lower_29of30=cp_lower(29,30)),
        representative_selection='Include 101, then minimum, lower-middle sorted rank, maximum of min(signed S0 loss, signed A0 loss); ascending seed breaks ties',
        representative_seeds=audit_result['representative_seeds'],contour_scope='Fixed upstream polygon library; inserting vertices is not contour convergence',
        numerical=dict(amplitude=.0002,phase_deg=1.,relative_group_loss=.05,residual=.01,right_reference=.001,power_balance=.002,
            max_extra_refinements=2,u_num_rule='sum of measured mesh, absorber and port-sampling differences; explicitly label representative-only coverage'),
        transient=dict(frequencies_hz=np.linspace(*source_band(F0),17).tolist(),cycles=5,window='Hann',initial_record_s=280e-6,hard_limit_s=560e-6,
            output='every internal step',initial_dt_s=1/(60*source_band(F0)[1]),half_step_required=True,tail=.0001,window_delta=.0002,arrival_relative=.05),
        layout=layout,material=read(ROOT/'config.json'),geometry_manifest_sha256=digest(OUT/'frozen_geometry_manifest.json'),
        single_controls_sha256=digest(OUT/'single_controls.json'),seed_manifest_sha256=digest(OUT/'independent_seeds.json'),
        source_geometry_sha256=digest(ROOT/'output/sample/geometry.json'),instruction_sha256=digest(ROOT/'200kHz_S0_A0_48mm_验证执行提示词.md'),
        experimental_validation='pending: no measured noise, sensor voltage model, material calibration or experimental specimens',
        solver_cores=2,old_P_queue_policy='Leave all existing jobs, limits, checkpoints and dependencies unchanged; serialize new FE work after protected jobs exit')
    freeze(protocol_path,protocol);write(OUT/'protocol.sha256.json',dict(sha256=digest(protocol_path)))
    return protocol

def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['audit','prepare','report','run'],default='audit',nargs='?');args=p.parse_args()
    if args.stage=='audit':audit();prepare()
    elif args.stage=='prepare':prepare()
    elif args.stage=='run':
        from detectability_worker import main as worker
        return worker()
    from detectability_report import render
    render()

if __name__=='__main__':main()
