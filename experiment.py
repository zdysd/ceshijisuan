"""Reproducible sample generation, calibration, diagnostics and analysis CLI."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
plt.rcParams['axes.unicode_minus']=False
import numpy as np
from shapely.geometry import Polygon

from cfrp_pores.pores import PoreSpec,random_field,ranked_mask,mask_statistics,export_sample,load_polygons,reposition,polygon_statistics,write_csv
from cfrp_pores.dispersion import Plate,validate_material
from cfrp_pores.attenuation import fit_attenuation,uncertainty


def preview(folder):
    folder=Path(folder)
    data=np.load(folder/'sample.npz'); meta=json.loads((folder/'metadata.json').read_text())
    size=meta['spec']['size_um']
    field,mask=data['field'],data['mask']
    if mask.ndim==3: field=field[:,:,field.shape[2]//2]; mask=mask[:,:,mask.shape[2]//2]
    fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    axs[0].imshow(field.T,origin='lower',extent=(0,size[0],0,size[1]),cmap='RdBu_r')
    axs[0].set_title('高斯随机场')
    axs[1].imshow(mask.T,origin='lower',extent=(0,size[0],0,size[1]),cmap='gray_r')
    axs[1].set_title(f'排序二值孔隙：P={100*meta["pixel_porosity"]:.3f}%')
    if (folder/'geometry.json').exists():
        polygons,_=load_polygons(folder/'geometry.json')
        for p in polygons:
            xy=np.asarray(p.exterior.coords); axs[2].fill(xy[:,0],xy[:,1],color='black')
            for ring in p.interiors:
                xy=np.asarray(ring.coords); axs[2].fill(xy[:,0],xy[:,1],color='white')
        axs[2].set_xlim(0,size[0]);axs[2].set_ylim(0,size[1]);axs[2].set_aspect('equal')
        axs[2].set_title(f'矢量孔隙轮廓：P={100*meta["vector_porosity"]:.3f}%')
    else:
        axs[2].text(.1,.5,'三维相关孔隙\n左侧显示中间截面\n完整孔形见 pores.stl');axs[2].axis('off')
    for ax in axs[:(3 if (folder/'geometry.json').exists() else 2)]:
        ax.set_xlabel('横向位置 x（μm）');ax.set_ylabel('纵向位置 y（μm）')
    fig.savefig(folder/'preview.png',dpi=160);plt.close(fig)


def calibration(spec,target_x_um,target_y_um,seeds=(101,102,103),iterations=3):
    """Simultaneously match axis-aligned component extents (explicit definition).

    Boundary-cut pores are excluded from size calibration, not from geometry.
    Do not claim these extents equal the paper's W_eq without its definition.
    Fixed padding preserves the same noise realization during parameter search.
    """
    records=[]; center=np.array(spec.correlation_um,float)
    padding=12*max(center)
    best=None
    for iteration in range(iterations):
        span=.6/(iteration+1)
        candidates=[center*np.array([a,b]) for a in [1-span,1,1+span] for b in [1-span,1,1+span]]
        for lengths in candidates:
            if min(lengths)<4*spec.pixel_um or max(lengths)>padding/3: continue
            sx,sy=[],[]
            for seed in seeds:
                trial=replace(spec,seed=seed,correlation_um=tuple(lengths))
                field=random_field(trial,padding_um=padding)
                rows,_=mask_statistics(ranked_mask(field,spec.porosity),trial)
                full=[r for r in rows if not r['boundary_truncated']]
                sx.extend(r['extent_0_um'] for r in full);sy.extend(r['extent_1_um'] for r in full)
            if not sx: continue
            mx,my=np.mean(sx),np.mean(sy)
            error_x=mx/target_x_um-1;error_y=my/target_y_um-1
            score=error_x**2+error_y**2
            rec=dict(iteration=iteration,a_um=float(lengths[0]),b_um=float(lengths[1]),
                     mean_x_extent_um=float(mx),mean_y_extent_um=float(my),relative_error_x=float(error_x),
                     relative_error_y=float(error_y),score=float(score),uncut_pores=len(sx))
            records.append(rec)
            if best is None or score<best['score']: best=rec
        if best is None: raise RuntimeError('No complete pores available for calibration.')
        center=np.array([best['a_um'],best['b_um']])
    return dict(best=best,metric='axis-aligned full component extent; not assumed identical to paper W_eq',
                seeds=list(seeds),target_x_um=target_x_um,target_y_um=target_y_um,
                matched_within_5_percent=max(abs(best['relative_error_x']),abs(best['relative_error_y']))<=.05),records


def generate(args):
    cfg=json.loads(Path(args.config).read_text())
    size=tuple(args.size or cfg['pore_roi_um'])
    corr=tuple(args.correlation or cfg['correlation_lengths_um'])
    spec=PoreSpec(size,args.pixel or cfg['pixel_um'],corr,args.porosity,args.seed,args.angle)
    summary=export_sample(args.output,spec);preview(args.output)
    if len(size)==2 and args.distributions and spec.porosity>0:
        polygons,_=load_polygons(Path(args.output)/'geometry.json')
        for name,cluster in [('dispersed',False),('clustered',True)]:
            moved=reposition(polygons,size,args.seed+1000,cluster)
            payload=[dict(exterior=np.asarray(p.exterior.coords).tolist(),interiors=[np.asarray(r.coords).tolist() for r in p.interiors]) for p in moved]
            folder=Path(args.output)/name;folder.mkdir(exist_ok=True)
            (folder/'geometry.json').write_text(json.dumps(dict(unit='um',size_um=size,polygons=payload),indent=2))
            write_csv(folder/'vector_pores.csv',polygon_statistics(moved))
            centers=np.array([[p.centroid.x,p.centroid.y] for p in moved])
            counts=np.histogramdd(centers,bins=[np.linspace(0,s,5) for s in size])[0]
            (folder/'metadata.json').write_text(json.dumps(dict(seed=args.seed+1000,distribution=name,
                vector_porosity=sum(p.area for p in moved)/np.prod(size),same_contours=True,
                quadrat_variance_to_mean=float(counts.var(ddof=1)/counts.mean())),indent=2))
    print(json.dumps(summary,indent=2))


def dispersion(args):
    cfg=json.loads(Path(args.config).read_text());C,rho=validate_material(cfg)
    plate=Plate(cfg['plate_thickness_m'],C[np.ix_([0,1,5],[0,1,5])],rho,args.elements)
    root=Path(args.output);root.mkdir(parents=True,exist_ok=True)
    rows=[]
    for f in np.linspace(5e4,5e5,24):
        for name in ['A0','S0']:
            d=plate.at_frequency(float(f),name)
            rows.append({k:v for k,v in d.items() if k not in ['y_m','u','v']})
    write_csv(root/'dispersion.csv',rows)
    for f in cfg['guided_frequencies_hz']:
        for name in ['A0','S0']:
            d=plate.at_frequency(f,name)
            write_csv(root/f'{name}_{f}_profile.csv',[dict(y_m=y,u_real=u.real,u_imag=u.imag,v_real=v.real,v_imag=v.imag) for y,u,v in zip(d['y_m'],d['u'],d['v'])])
    plot_dispersion(rows,root)
    (root/'metadata.json').write_text(json.dumps(dict(material_status=cfg['material_status'],thickness_m=plate.thickness_m,
           discretization='linear SAFE through thickness',elements=args.elements,plane_strain=True,
           note='Fundamental branches only; not a proof that higher modes are cut off.'),indent=2))


def plot_dispersion(rows,folder):
    """Render saved dispersion values without recalculating wave solutions."""
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for name,label in [('A0','A0（反对称模态）'),('S0','S0（对称模态）')]:
        selected=[r for r in rows if r['mode']==name]
        for ax,key,ylabel,title in zip(axs,['phase_velocity_m_s','group_velocity_m_s'],
                                      ['相速度（m/s）','群速度（m/s）'],
                                      ['基阶模态相速度频散','基阶模态群速度频散']):
            ax.plot([float(r['frequency_hz'])/1000 for r in selected],[float(r[key]) for r in selected],label=label)
            ax.set_xlabel('频率（kHz）');ax.set_ylabel(ylabel);ax.set_title(title);ax.legend();ax.grid(alpha=.3)
    fig.savefig(Path(folder)/'dispersion.png',dpi=180);plt.close(fig)


def analyze(args):
    import csv
    groups={}
    with Path(args.input).open(encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            key=(row['sample'],row['mode'],float(row['frequency_hz']))
            groups.setdefault(key,[]).append(row)
    rows=[]
    for (sample,mode,f),group in groups.items():
        group=sorted(group,key=lambda r:float(r['x_m']))
        x=[float(r['x_m']) for r in group]
        z=[complex(float(r['real']),float(r['imag'])) for r in group]
        has_ref='reference_real' in group[0]
        ref=[complex(float(r['reference_real']),float(r['reference_imag'])) for r in group] if has_ref else None
        rows.append(dict(sample=sample,mode=mode,frequency_hz=f,**fit_attenuation(x,z,ref)))
    root=Path(args.output);root.mkdir(parents=True,exist_ok=True)
    write_csv(root/'attenuation.csv',rows)
    summary=[]
    for mode,f in sorted(set((r['mode'],r['frequency_hz']) for r in rows)):
        values=[r['alpha_np_m'] for r in rows if r['mode']==mode and r['frequency_hz']==f]
        if len(values)>1: summary.append(dict(mode=mode,frequency_hz=f,**uncertainty(values)))
    write_csv(root/'uncertainty.csv',summary)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default='config.json')
    subs=p.add_subparsers(dest='command',required=True)
    g=subs.add_parser('generate');g.add_argument('--size',nargs='+',type=float);g.add_argument('--correlation',nargs='+',type=float)
    g.add_argument('--pixel',type=float);g.add_argument('--porosity',type=float,default=.015)
    g.add_argument('--seed',type=int,default=101);g.add_argument('--angle',type=float,default=0)
    g.add_argument('--distributions',action='store_true');g.add_argument('--output',default='output/sample')
    d=subs.add_parser('dispersion');d.add_argument('--elements',type=int,default=60);d.add_argument('--output',default='output/dispersion')
    a=subs.add_parser('analyze');a.add_argument('input');a.add_argument('--output',default='output/analysis')
    c=subs.add_parser('calibrate');c.add_argument('--target-x',type=float,required=True);c.add_argument('--target-y',type=float,required=True)
    c.add_argument('--output',default='output/calibration')
    args=p.parse_args()
    if args.command=='generate':generate(args)
    elif args.command=='dispersion':dispersion(args)
    elif args.command=='analyze':analyze(args)
    else:
        cfg=json.loads(Path(args.config).read_text());spec=PoreSpec(tuple(cfg['pore_roi_um']),cfg['pixel_um'],tuple(cfg['correlation_lengths_um']))
        summary,records=calibration(spec,args.target_x,args.target_y)
        folder=Path(args.output);folder.mkdir(parents=True,exist_ok=True)
        write_csv(folder/'search.csv',records);(folder/'calibration.json').write_text(json.dumps(summary,indent=2))
        print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
