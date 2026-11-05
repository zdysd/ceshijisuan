"""TetGen PLC volume mesh: preserve Gaussian pore surfaces, bypass CAD fitting."""
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components
import meshpy.tet as tet


def volume_mesh(folder,size_um=(600,300,300),edge_um=50):
    folder=Path(folder);path=folder/'pores.stl'
    triangles=(np.array([[float(v) for v in line.split()[1:]] for line in path.read_text().splitlines() if line.startswith('vertex ')]).reshape(-1,3,3)*1e6) if path.exists() else np.empty((0,3,3))
    points,indices=np.unique(np.round(triangles.reshape(-1,3),8),axis=0,return_inverse=True)
    faces=indices.reshape(-1,3)
    graph=coo_array((np.ones(faces.size*2),(
        np.concatenate([faces[:,0],faces[:,1],faces[:,2],faces[:,1],faces[:,2],faces[:,0]]),
        np.concatenate([faces[:,1],faces[:,2],faces[:,0],faces[:,0],faces[:,1],faces[:,2]]))),shape=(len(points),len(points))).tocsr()
    n,labels=connected_components(graph,directed=False)
    holes=[]
    for label in range(n):
        face=faces[labels[faces[:,0]]==label]
        vertices=points[face]
        volume=np.einsum('ij,ij->i',vertices[:,0],np.cross(vertices[:,1],vertices[:,2])).sum()/6
        if volume>0:
            a,b,c=vertices[0];normal=np.cross(b-a,c-a);normal/=np.linalg.norm(normal)
            holes.append((a+b+c)/3-normal*1e-3)
    lx,ly,lz=size_um
    cube=np.array([[0,0,0],[lx,0,0],[lx,ly,0],[0,ly,0],[0,0,lz],[lx,0,lz],[lx,ly,lz],[0,ly,lz]],float)
    outer=[[0,3,2],[0,2,1],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[3,7,6],[3,6,2],[0,4,7],[0,7,3],[1,2,6],[1,6,5]]
    all_points=np.vstack([cube,points]);facets=outer+(faces[:,::-1]+8).tolist()
    markers=[1,1,2,2,3,3,4,4,5,5,6,6]+[7+int(labels[f[0]]) for f in faces]
    info=tet.MeshInfo();info.set_points(all_points);info.set_facets(facets,markers=markers);info.set_holes(holes)
    print('TetGen PLC',len(points),'surface vertices,',len(holes),'closed cavities',flush=True)
    mesh=tet.build(info,options=tet.Options('pq2.0f'),max_volume=edge_um**3/6,verbose=False)
    nodes=np.asarray(mesh.points);elements=np.asarray(mesh.elements,dtype=int)
    corners=nodes[elements]
    volumes=abs(np.einsum('ij,ij->i',corners[:,1]-corners[:,0],np.cross(corners[:,2]-corners[:,0],corners[:,3]-corners[:,0])))/6
    if np.any(volumes<=0):raise RuntimeError('TetGen produced zero-volume cells.')
    phi=1-volumes.sum()/np.prod(size_um)
    target=json.loads((folder/'metadata.json').read_text())['stl_porosity']
    if abs(phi-target)>1e-5:raise RuntimeError('Tetrahedral volume does not preserve the cavity volume.')
    # Explicit NASTRAN boundary triangles preserve outer face markers.
    boundary_faces=np.asarray(mesh.faces);boundary_markers=np.asarray(mesh.face_markers)
    bdf=folder/'solid_volume.bdf'
    with bdf.open('w',encoding='ascii') as stream:
        stream.write('BEGIN BULK\n')
        for i,node in enumerate(nodes*1e-6,1):stream.write(f'GRID,{i},,'+','.join(f'{v:.12g}' for v in node)+'\n')
        for i,tetra in enumerate(elements+1,1):stream.write(f'CTETRA,{i},1,'+','.join(map(str,tetra))+'\n')
        offset=len(elements)+1
        for face,marker in zip(boundary_faces,boundary_markers):
            if marker:
                stream.write(f'CTRIA3,{offset},{int(marker)+1},'+','.join(map(str,face+1))+'\n');offset+=1
        stream.write('ENDDATA\n')
    np.savez_compressed(folder/'tetrahedra.npz',nodes_m=nodes*1e-6,elements=elements)
    (folder/'volume_mesh.json').write_text(json.dumps(dict(nodes=len(nodes),tetrahedra=len(elements),porosity=phi,
        cavities=len(holes),min_cell_volume_um3=float(volumes.min()),method='TetGen PLC; no pore deletion; no CAD surface fitting'),indent=2))
    return bdf
