"""Cross-check least-squares and reciprocal projections on saved A0 fields."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import json,sys
import numpy as np
from guided_collective import settings,OUT,dephase_average
from cfrp_pores.comsol import interpolate
from cfrp_pores.dispersion import modal_direction_basis,separate_modal_directions,reciprocal_modal_directions
sys.stdout.reconfigure(encoding='utf-8')


def main():
    s=settings();f=200000.;y=np.linspace(1e-8,.001-1e-8,801)
    basis,labels,scale=modal_direction_basis(y-.0005,s['plate'],f)
    k=s['plate'].at_frequency(f,'A0')['k_rad_m']
    client=mph.Client(cores=4,version='6.2');rows=[];references={}
    for name in ['reference','d20/L8/seed_101','d40/L8/seed_101','d80/single_maximum']:
        path=OUT/name/'model.mph'
        if not path.exists():continue
        model=client.load(str(path));projection={'least_squares':[],'reciprocity':[]}
        try:
            g=model.java.component('comp1').geom('geom1');measure=g.measureFinal();measure.selection().geom('geom1',2);measure.selection().all()
            void_area=s['length']*s['h']-float(measure.getArea())
            for x in s['pre_positions']+s['post_positions']:
                data=interpolate(model,[np.full(len(y),x),y],['u','v','solid.sx','solid.sxy'])
                field=[v[-1] for v in data]
                ls,residual=separate_modal_directions(y-.0005,*field,basis,labels,scale)
                reciprocal=reciprocal_modal_directions(y-.0005,*field,basis,labels,scale)
                projection['least_squares'].append(ls['A0_forward']);projection['reciprocity'].append(reciprocal['A0_forward'])
        finally:client.remove(model)
        normalized={}
        for method,values in projection.items():
            values=np.array(values);normalized[method]=dephase_average(values[3:],s['post_positions'],k)/dephase_average(values[:3],s['pre_positions'],k)
            if name=='reference':references[method]=normalized[method]
            normalized[method]/=references[method]
        requested=0.
        metadata_path=OUT/name/'geometry/metadata.json'
        if metadata_path.exists():
            meta=json.loads(metadata_path.read_text());requested=meta['actual_porosity']*meta['length_mm']*1e-6
        row=dict(case=name,frequency_hz=f,mode='A0',measured_void_area_m2=void_area,requested_void_area_m2=requested,
                 relative_void_area_error=float(abs(void_area/requested-1)) if requested else 0.,
                 least_squares_reduction=float(1-abs(normalized['least_squares'])),
                 reciprocity_reduction=float(1-abs(normalized['reciprocity'])),
                 absolute_magnitude_difference=float(abs(abs(normalized['least_squares'])-abs(normalized['reciprocity']))))
        rows.append(row);print(json.dumps(row),flush=True)
    (OUT/'projection_crosscheck.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')


if __name__=='__main__':main()
