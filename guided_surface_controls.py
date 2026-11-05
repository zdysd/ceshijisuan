"""Check the largest pore near both plate surfaces, avoiding neutral-axis bias."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
from shapely import affinity
from guided_collective import settings,OUT,grouped_pores,run_case,export_geometry


def main():
    s=settings();client=mph.Client(cores=2,version='6.2')
    for d in [20,40,80]:
        pores,meta=grouped_pores(d,8,101);pore=max(pores,key=lambda p:p.area)
        for label,offset in [('upper',950-pore.bounds[3]),('lower',50-pore.bounds[1])]:
            shifted=affinity.translate(pore,yoff=offset)
            folder=OUT/f'd{d}/single_surface_{label}'
            export_geometry(folder/'geometry',[shifted],dict(meta,single_case=label,pore_count=1,
                actual_porosity=pore.area/8e6,note='Largest contour retained; translated to 50 um surface ligament.'))
            run_case(client,folder,s,[shifted],repair_tolerance_m=1e-9)


if __name__=='__main__':main()
