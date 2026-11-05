from dataclasses import asdict, dataclass
import csv
import json
from pathlib import Path

import numpy as np
from scipy import ndimage, fft
from scipy.spatial import cKDTree
from shapely import affinity
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
from skimage import measure


@dataclass(frozen=True)
class PoreSpec:
    size_um: tuple = (1000., 1000.)  # x,y[,z], in physical units
    pixel_um: float = 2.
    correlation_um: tuple = (20., 10.)
    porosity: float = .015  # fraction, NOT percent
    seed: int = 101
    angle_deg: float = 0.

    def validate(self):
        if len(self.size_um) not in (2, 3) or len(self.correlation_um) != len(self.size_um):
            raise ValueError('Use two or three physical dimensions and matching correlations.')
        if self.pixel_um <= 0 or min(self.size_um) <= 0 or min(self.correlation_um) <= 0:
            raise ValueError('Lengths must be positive.')
        if not 0 <= self.porosity <= .5:
            raise ValueError('Porosity must be a fraction between 0 and 0.5.')
        if self.pixel_um > min(self.correlation_um) / 4:
            raise ValueError('Random field needs at least four pixels per shortest correlation length.')
        if len(self.size_um) == 3 and self.angle_deg != 0:
            raise ValueError('3D uses principal-axis correlations; arbitrary 3D rotations are not implemented.')
        for size in self.size_um:
            if not np.isclose(size / self.pixel_um, round(size / self.pixel_um)):
                raise ValueError('Physical dimensions must be integer multiples of pixel_um.')


def random_field(spec, padding_um=None):
    """Real Gaussian field with R(r)=exp(-sum((r_i/a_i)^2)).

    Array axes are physical x,y[,z]. Spectral filtering of real white noise
    preserves Hermitian symmetry. Padding >3 longest correlations prevents
    periodic wrap-around in the central specimen; it is NOT a periodic cell.
    """
    spec.validate()
    shape = tuple(round(s / spec.pixel_um) for s in spec.size_um)
    padding_um=max(spec.correlation_um)*3 if padding_um is None else padding_um
    if padding_um<3*max(spec.correlation_um): raise ValueError('Padding must be at least three longest correlations.')
    margin = int(np.ceil(padding_um / spec.pixel_um))
    padded = tuple(fft.next_fast_len(n + 2 * margin) for n in shape)
    if np.prod(padded) > 20_000_000:
        raise ValueError('Padded field exceeds 20 million cells; reduce the local domain or increase pixel size with a resolution check.')
    rng = np.random.default_rng(spec.seed)
    noise = rng.standard_normal(padded)
    ks = np.meshgrid(*[2*np.pi*fft.fftfreq(n, d=spec.pixel_um) for n in padded], indexing='ij', sparse=True)
    if len(shape) == 2:
        theta = np.deg2rad(spec.angle_deg)
        ku = np.cos(theta)*ks[0] + np.sin(theta)*ks[1]
        kv = -np.sin(theta)*ks[0] + np.cos(theta)*ks[1]
        exponent = (spec.correlation_um[0]*ku)**2 + (spec.correlation_um[1]*kv)**2
    else:
        exponent = sum((a*k)**2 for a, k in zip(spec.correlation_um, ks))
    # Squared transfer is exp(-a^2 k^2/4), the target covariance spectrum.
    field = fft.ifftn(fft.fftn(noise) * np.exp(-exponent/8)).real
    crop = tuple(slice((n-s)//2, (n-s)//2+s) for n, s in zip(padded, shape))
    field = field[crop].copy()
    return (field-field.mean()) / field.std()


def ranked_mask(field, porosity):
    if not np.isfinite(field).all() or not 0 <= porosity <= .5:
        raise ValueError('Finite field and porosity in [0,.5] required.')
    count = int(round(field.size * porosity))
    mask = np.zeros(field.size, dtype=bool)
    if count:
        # stable tie-breaking, reproducible even if a field has repeated values
        mask[np.argsort(field.ravel(), kind='stable')[:count]] = True
    return mask.reshape(field.shape)


def mask_statistics(mask, spec):
    labels, n = ndimage.label(mask)  # face connectivity, shared edges only
    objects = ndimage.find_objects(labels)
    rows, centers = [], []
    for label, sl in enumerate(objects, 1):
        coords = np.argwhere(labels[sl] == label)
        coords += np.array([s.start for s in sl])
        points = (coords + .5) * spec.pixel_um
        center = points.mean(axis=0)
        centers.append(center)
        extent = (np.ptp(coords, axis=0)+1)*spec.pixel_um
        # Pixel cells have finite moment of inertia, even for a one-pixel pore.
        covariance = np.atleast_2d(np.cov(points, rowvar=False, bias=True)) if len(points) > 1 else np.zeros((mask.ndim, mask.ndim))
        covariance += np.eye(mask.ndim)*spec.pixel_um**2/12
        eig, vectors = np.linalg.eigh(covariance)
        axes = np.sqrt((16 if mask.ndim == 2 else 20)*eig)
        size = len(coords)*spec.pixel_um**mask.ndim
        equivalent = 2*np.sqrt(size/np.pi) if mask.ndim == 2 else (6*size/np.pi)**(1/3)
        touches = any(np.any(coords[:, d] == 0) or np.any(coords[:, d] == mask.shape[d]-1) for d in range(mask.ndim))
        row = dict(pore_id=label, pixel_count=len(coords), equivalent_diameter_um=float(equivalent),
                   aspect_ratio=float(axes[-1]/axes[0]), boundary_truncated=bool(touches))
        row['area_um2' if mask.ndim == 2 else 'volume_um3'] = float(size)
        for d in range(mask.ndim):
            row[f'center_{d}_um'] = float(center[d])
            row[f'extent_{d}_um'] = float(extent[d])
        if mask.ndim == 2:
            row['orientation_deg'] = float(np.rad2deg(np.arctan2(vectors[1,-1],vectors[0,-1])) % 180)
        rows.append(row)
    centers = np.array(centers)
    nearest = [] if n < 2 else cKDTree(centers).query(centers, k=2)[0][:,1].tolist()
    # This is a descriptor, not a universal significance test for clustering.
    bins = [np.linspace(0, s, 5) for s in spec.size_um]
    quadrat = np.histogramdd(centers, bins=bins)[0].ravel() if n else np.zeros(4**mask.ndim)
    diameters = [r['equivalent_diameter_um'] for r in rows]
    summary = dict(spec=asdict(spec), pixel_porosity=float(mask.mean()), pore_count=n,
                   boundary_pore_count=sum(r['boundary_truncated'] for r in rows),
                   diameter_mean_um=float(np.mean(diameters)) if n else None,
                   diameter_p95_um=float(np.percentile(diameters,95)) if n else None,
                   diameter_max_um=max(diameters) if n else None,
                   nearest_neighbor_mean_um=float(np.mean(nearest)) if nearest else None,
                   quadrat_variance_to_mean=float(quadrat.var(ddof=1)/quadrat.mean()) if n else None)
    return rows, summary


def axial_correlations(field, mask, pixel_um, max_lag=80):
    """Nonperiodic, overlap-normalized correlations of field and binary pores."""
    rows = []
    for axis in range(field.ndim):
        for lag in range(min(max_lag, field.shape[axis]//3)+1):
            for name, data in [('field', field), ('binary', mask.astype(float))]:
                centered = data-data.mean()
                variance = np.mean(centered**2)
                left, right = [slice(None)]*data.ndim, [slice(None)]*data.ndim
                if lag:
                    left[axis], right[axis] = slice(None,-lag),slice(lag,None)
                value = float(np.mean(centered[tuple(left)]*centered[tuple(right)])/variance) if variance else None
                rows.append(dict(axis=axis,lag_um=lag*pixel_um,kind=name,correlation=value))
    return rows


def contours_from_field(field, threshold, pixel_um, size_um):
    """Iso-contours of the smooth field, preserving islands and boundary cuts."""
    if not np.isfinite(threshold):
        return []
    outside = max(float(field.max()+1), threshold+1)
    contours = measure.find_contours(np.pad(field,1,constant_values=outside), threshold, fully_connected='low')
    rings = []
    for contour in contours:
        vertices = (contour-.5)*pixel_um
        polygon = Polygon(vertices)
        if not polygon.is_valid:
            raise ValueError('Invalid extracted contour; refine field resolution instead of silently repairing it.')
        if polygon.area:
            rings.append(polygon)
    # Nesting depth distinguishes void shells from solid islands in a void.
    rings.sort(key=lambda p:p.area, reverse=True)
    parents = []
    for i, ring in enumerate(rings):
        containing = [j for j in range(i) if rings[j].contains(ring.representative_point())]
        parents.append(min(containing,key=lambda j:rings[j].area) if containing else None)
    depth = []
    for parent in parents:
        depth.append(0 if parent is None else depth[parent]+1)
    holes = []
    for i, ring in enumerate(rings):
        if depth[i] % 2:
            continue
        interiors = [list(rings[j].exterior.coords) for j in range(len(rings)) if parents[j] == i and depth[j]%2]
        p = Polygon(ring.exterior.coords,interiors).intersection(box(0,0,*size_um))
        holes.extend([p] if p.geom_type == 'Polygon' else list(p.geoms))
    return [p for p in holes if p.area > 0]


def geometry_from_field(field, spec, tolerance=1e-4):
    """Adjust iso-level to match actual vector area; do not confuse pixel fraction with vector area."""
    if field.ndim != 2:
        raise ValueError('2D vector contours required.')
    if spec.porosity == 0:
        return [], None
    target_area = np.prod(spec.size_um)*spec.porosity
    lo, hi = float(field.min()-1e-6),float(np.quantile(field,min(.99,spec.porosity+.02)))
    best = []
    for _ in range(30):
        q = (lo+hi)/2
        polygons = contours_from_field(field,q,spec.pixel_um,spec.size_um)
        area = sum(p.area for p in polygons)
        best = polygons
        if abs(area-target_area)/np.prod(spec.size_um) <= tolerance:
            return polygons,q
        if area < target_area:
            lo = q
        else:
            hi = q
    raise RuntimeError('Vector porosity did not converge; no geometry is exported silently.')


def polygon_statistics(polygons):
    rows = []
    for i,p in enumerate(polygons):
        rect = np.array(p.minimum_rotated_rectangle.exterior.coords)
        sides = np.linalg.norm(np.diff(rect,axis=0),axis=1)
        rows.append(dict(pore_id=i+1,area_um2=p.area,equivalent_diameter_um=2*np.sqrt(p.area/np.pi),
                         length_um=float(sides.max()),width_um=float(sides.min()),
                         centroid_x_um=p.centroid.x,centroid_y_um=p.centroid.y,
                         vertex_count=len(p.exterior.coords),solid_islands=len(p.interiors)))
    return rows


def simplify_geometry(polygons,tolerance_um):
    simplified=[p.simplify(tolerance_um,preserve_topology=True) for p in polygons]
    if any(not p.is_valid or p.is_empty for p in simplified):
        raise ValueError('Simplification changed geometry validity.')
    if any(len(a.interiors)!=len(b.interiors) for a,b in zip(polygons,simplified)):
        raise ValueError('Simplification changed solid-island topology.')
    if any(a.intersects(b) for i,a in enumerate(simplified) for b in simplified[:i]):
        raise ValueError('Simplification caused pores to overlap.')
    original=sum(p.area for p in polygons)
    change=(sum(p.area for p in simplified)-original)/original if original else 0.
    return simplified,dict(tolerance_um=tolerance_um,relative_void_area_change=change,
                          original_vertices=sum(len(p.exterior.coords) for p in polygons),
                          simplified_vertices=sum(len(p.exterior.coords) for p in simplified))


def reposition(polygons, size_um, seed, clustered=False, max_attempts=100000):
    """Translate identical contours without changing area, orientation or size."""
    rng = np.random.default_rng(seed)
    domain = box(0,0,*size_um)
    placed = []
    clusters = rng.uniform(.2,.8,(3,2))*np.array(size_um)
    for p in sorted(polygons,key=lambda p:p.area,reverse=True):
        centered = affinity.translate(p,-p.centroid.x,-p.centroid.y)
        for _ in range(max_attempts):
            xy = (clusters[rng.integers(3)]+rng.normal(size=2)*np.array(size_um)*.09) if clustered else rng.uniform(size=2)*size_um
            candidate = affinity.translate(centered,*xy)
            if domain.contains(candidate) and all(candidate.distance(other)>1e-6 for other in placed):
                placed.append(candidate)
                break
        else:
            raise RuntimeError('Cannot place the same contours without overlap; do not drop or resize pores.')
    return placed


def write_csv(path, rows):
    path = Path(path)
    if not rows:
        path.write_text('',encoding='utf-8')
        return
    with path.open('w',newline='',encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export_sample(path, spec, field=None):
    path=Path(path); path.mkdir(parents=True,exist_ok=True)
    field=random_field(spec) if field is None else field
    mask=ranked_mask(field,spec.porosity)
    rows,summary=mask_statistics(mask,spec)
    np.savez_compressed(path/'sample.npz',field=field,mask=mask)
    write_csv(path/'pixel_pores.csv',rows)
    write_csv(path/'correlations.csv',axial_correlations(field,mask,spec.pixel_um))
    if mask.ndim == 2:
        polygons,q=geometry_from_field(field,spec)
        payload=[dict(exterior=np.asarray(p.exterior.coords).tolist(),interiors=[np.asarray(r.coords).tolist() for r in p.interiors]) for p in polygons]
        (path/'geometry.json').write_text(json.dumps(dict(unit='um',size_um=spec.size_um,polygons=payload),indent=2),encoding='utf-8')
        write_csv(path/'vector_pores.csv',polygon_statistics(polygons))
        summary.update(vector_porosity=sum(p.area for p in polygons)/np.prod(spec.size_um),vector_threshold=q,vector_pore_count=len(polygons))
    else:
        volume,q=export_gaussian_stl(path/'pores.stl',field,spec)
        summary['stl_porosity']=volume/(np.prod(spec.size_um)*1e-18)
        summary['stl_porosity_matches_target']=bool(abs(summary['stl_porosity']-spec.porosity)<=1e-4)
        summary['stl_threshold']=q
        summary['stl_note']='Smooth Gaussian-field iso-surface, volume-calibrated separately from ranked voxels; boundary closures lie inside the finite ROI. Validate field resolution and actual COMSOL volume.'
    (path/'metadata.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return summary


def export_stl(path, mask, pixel_um):
    """Closed voxel-derived 3D surfaces in metres, not a 2D extrusion."""
    if not mask.any():
        return 0.
    vertices,faces,_,_=measure.marching_cubes(np.pad(mask.astype(float),1),.5,spacing=(pixel_um*1e-6,)*3,allow_degenerate=False)
    vertices=vertices.astype(np.float64)
    vertices-=pixel_um*1e-6*.5
    # Binary iso-surfaces lie on a half-voxel lattice. Snap float32 marching
    # output back to that lattice to avoid artificial picometre CAD edges.
    step=pixel_um*1e-6/2
    vertices=np.round(vertices/step)*step
    unique,inverse=np.unique(vertices,axis=0,return_inverse=True)
    faces=inverse[faces];vertices=unique
    lengths=np.linalg.norm(np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]]),axis=1)
    faces=faces[lengths>step**2*1e-10]
    edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
    _,counts=np.unique(edges,axis=0,return_counts=True)
    if np.any(counts!=2):raise ValueError('3D pore surface is nonmanifold; refine the field rather than patching away pores.')
    signed=np.einsum('ij,ij->i',vertices[faces[:,0]],np.cross(vertices[faces[:,1]],vertices[faces[:,2]])).sum()/6
    if signed<0:
        faces=faces[:,[0,2,1]]  # outward shells, inward surfaces around solid islands
    volume=abs(signed)
    with Path(path).open('w',encoding='ascii') as stream:
        stream.write('solid random_pores\n')
        for face in faces:
            a,b,c=vertices[face]
            normal=np.cross(b-a,c-a); length=np.linalg.norm(normal)
            normal=normal/length if length else normal
            stream.write('facet normal '+ ' '.join(map(str,normal))+'\nouter loop\n')
            for v in [a,b,c]: stream.write('vertex '+' '.join(map(str,v))+'\n')
            stream.write('endloop\nendfacet\n')
        stream.write('endsolid random_pores\n')
    return float(volume)


def export_gaussian_stl(path,field,spec,tolerance=1e-4):
    """3D analogue of continuous-field contour calibration, no voxel staircase."""
    if spec.porosity==0:return 0.,None
    # Cell-centred cubic refinement resolves smooth minima before extracting
    # surfaces; piecewise-linear cones at coarse extrema make spurious tiny
    # CAD edges. The same realization is refined, no new random field is drawn.
    refinement=2
    field=ndimage.zoom(field,refinement,order=3,mode='reflect',grid_mode=True)
    geometry_pixel_um=spec.pixel_um/refinement
    outside=float(field.max()+5*np.ptp(field))
    padded=np.pad(field,1,constant_values=outside)
    target=spec.porosity*np.prod(spec.size_um)*1e-18
    lo=float(field.min()+1e-8);hi=float(np.quantile(field,min(.49,spec.porosity+.1)))
    for _ in range(35):
        q=(lo+hi)/2
        vertices,faces,_,_=measure.marching_cubes(padded,q,allow_degenerate=False)
        vertices=(vertices.astype(float)-.5)*geometry_pixel_um*1e-6
        signed=np.einsum('ij,ij->i',vertices[faces[:,0]],np.cross(vertices[faces[:,1]],vertices[faces[:,2]])).sum()/6
        volume=abs(signed)
        if abs(volume-target)/(np.prod(spec.size_um)*1e-18)<=tolerance:break
        if volume<target:lo=q
        else:hi=q
    else:raise RuntimeError('Gaussian 3D iso-volume did not meet porosity tolerance.')
    if signed<0:faces=faces[:,[0,2,1]]
    edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
    if np.any(np.unique(edges,axis=0,return_counts=True)[1]!=2):raise ValueError('Gaussian pore surface is not closed and manifold.')
    if vertices.min() < -1e-12 or np.any(vertices.max(axis=0)>np.array(spec.size_um)*1e-6+1e-12):
        raise ValueError('Closed Gaussian pores exceed their physical ROI.')
    with Path(path).open('w',encoding='ascii') as stream:
        stream.write('solid gaussian_pores\n')
        for face in faces:
            a,b,c=vertices[face];normal=np.cross(b-a,c-a);normal/=np.linalg.norm(normal)
            stream.write('facet normal '+' '.join(map(str,normal))+'\nouter loop\n')
            for v in [a,b,c]:stream.write('vertex '+' '.join(map(str,v))+'\n')
            stream.write('endloop\nendfacet\n')
        stream.write('endsolid gaussian_pores\n')
    return float(volume),q


def load_polygons(path):
    from .portability import project_path
    data=json.loads(project_path(path).read_text(encoding='utf-8'))
    if data.get('unit')!='um':raise ValueError('Geometry must explicitly use micrometres.')
    size=tuple(data['size_um'])
    if len(size)!=2 or min(size)<=0:raise ValueError('Positive 2D geometry dimensions required.')
    polygons=[Polygon(p['exterior'],p['interiors']) for p in data['polygons']]
    if any(not p.is_valid or p.is_empty for p in polygons):raise ValueError('Invalid pore geometry.')
    domain=box(0,0,*size).buffer(1e-7)
    if any(not domain.covers(p) for p in polygons):raise ValueError('Pore coordinates exceed declared region.')
    if any(a.intersection(b).area>1e-10 for i,a in enumerate(polygons) for b in polygons[:i]):raise ValueError('Input pores overlap.')
    return polygons,size
