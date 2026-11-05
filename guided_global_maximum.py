"""Test the actual largest contour across each completed diameter-group family."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
from shapely import affinity
from cfrp_pores.pores import load_polygons
from guided_collective import OUT,settings,read,run_case,export_geometry,grouped_pores


def main():
    s=settings();client=mph.Client(cores=2,version='6.2')
    for d in [20,40,80]:
        candidates=[]
        for seed in range(101,113 if d==80 else 107):
            for length in ([8,24,48] if d==80 else [8]):
                base=OUT if seed<=106 else OUT/'extended_samples'
                source=base/f'd{d}/L{length}/seed_{seed}'
                if (source/'geometry/metadata.json').exists():
                    meta=read(source/'geometry/metadata.json');pores,_=load_polygons(source/'geometry/geometry.json')
                else:pores,meta=grouped_pores(d,length,seed)
                candidates.append((meta['diameter_max_um'],meta,max(pores,key=lambda p:p.area),source))
        _,meta,pore,source=max(candidates,key=lambda row:row[0])
        for label,y in [('center',500-pore.centroid.y),('upper',950-pore.bounds[3]),('lower',50-pore.bounds[1])]:
            shifted=affinity.translate(pore,xoff=4000-pore.centroid.x,yoff=y)
            folder=OUT/f'd{d}/global_maximum_{label}'
            export_geometry(folder/'geometry',[shifted],dict(meta,length_mm=8,pore_count=1,
                actual_porosity=pore.area/8e6,source_case=str(source.relative_to(OUT)),
                single_case=label,note='Largest actual equivalent-diameter contour among all generated group cases.'))
            run_case(client,folder,s,[shifted],repair_tolerance_m=1e-9)


if __name__=='__main__':main()
