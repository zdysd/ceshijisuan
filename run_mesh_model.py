"""COMSOL on imported volume meshes; bypass CAD fitting of pore surfaces."""
import mph,site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))
import json,sys,jpype
import numpy as np
from cfrp_pores.comsol import _box_selection,add_frequency_study,save_and_solve,interpolate,start_client
from cfrp_pores.dispersion import validate_material
sys.stdout.reconfigure(encoding='utf-8')


def solve_mesh_model(root,output,config,frequency=5e6,size=(.0006,.0003,.0003)):
    root=Path(root);output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    C,rho=validate_material(config);client=start_client();model=client.create('CFRP_3D_mesh_based');j=model.java
    c=j.component().create('comp1',True);g=c.geom().create('geom1',3);g.lengthUnit('m')
    mesh=c.mesh().create('mesh1','geom1');mesh.geometricModel(True)
    imp=mesh.create('import','Import');imp.set('filename',str((root/'solid_volume.bdf').resolve()));mesh.run();c.geometricModel('mesh1')
    print('mesh model',c.geometricModel(),mesh.getNumElem(),flush=True)
    p=c.physics().create('solid','SolidMechanics','geom1');le=p.feature('lemm1')
    le.set('SolidModel','Anisotropic');le.set('AnisotropicOption','AnisotropicVo');le.set('DVo_mat','userdef')
    le.set('DVo',jpype.JArray(jpype.JString)([str(v)+'[Pa]' for v in C.ravel(order='F')]))
    le.set('rho_mat','userdef');le.set('rho',str(rho));le.set('geometricNonlinearity','linear')
    for tag,axis,value in [('left',0,0.),('right',0,size[0]),('ymin',1,0.),('ymax',1,size[1]),('zmin',2,0.),('zmax',2,size[2])]:
        bounds=[(-1e-9,s+1e-9) for s in size];bounds[axis]=(value-1e-9,value+1e-9)
        _box_selection(c,tag,2,bounds)
        ids=list(c.selection(tag).entities(2))
        if not ids:raise RuntimeError('Imported mesh boundary is not selectable: '+tag)
        if tag=='left':
            load=p.create('source','BoundaryLoad',2);load.selection().named(tag);load.set('FperArea',['1e6[Pa]','0','0'])
        else:
            f=p.create('abs_'+tag,'LowReflectingBoundary',2);f.selection().named(tag)
            mapping={(0,1):5,(0,2):4,(1,2):3};z=[]
            for d in range(3):
                v=C[axis,axis] if d==axis else C[mapping[tuple(sorted((d,axis)))],mapping[tuple(sorted((d,axis)))]]
                z.append(str(np.sqrt(rho*v))+'[Pa*s/m]')
            f.set('DampingType','userdef');f.set('coordinateSystem','GlobalSystem');f.set('d_i',z)
    add_frequency_study(model,frequency)
    try:
        info=save_and_solve(model,output)
        values=interpolate(model,[[.2*size[0],.5*size[0],.8*size[0]],[size[1]/2]*3,[size[2]/2]*3],['u','v','w'])
        if not np.isfinite(values).all():raise RuntimeError('Imported-mesh solution has invalid receiver displacements.')
        info.update(actual_porosity=json.loads((root/'volume_mesh.json').read_text())['porosity'],mesh_based_geometry=True,
                    status='solved_local_method_validation',geometry_note='Same Gaussian surfaces; TetGen volume mesh, no CAD fitting.')
        np.savez_compressed(output.with_suffix('.receivers.npz'),displacements=values)
        return info
    finally:client.remove(model)


if __name__=='__main__':
    cfg=json.loads(Path('config.json').read_text())
    info=solve_mesh_model('output/comsol_3d_refined/random_geometry','output/comsol_3d_mesh_based/random_pores.mph',cfg)
    Path('output/comsol_3d_mesh_based/results.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
