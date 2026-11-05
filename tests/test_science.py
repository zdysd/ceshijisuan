import json
import unittest
import tempfile
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon,Point
from cfrp_pores.pores import PoreSpec,random_field,ranked_mask,geometry_from_field,contours_from_field,reposition,mask_statistics,axial_correlations,export_stl,export_gaussian_stl
from cfrp_pores.attenuation import fit_attenuation,converged_alpha,uncertainty,gated_spectra,coherent_mean,forward_bulk_displacement
from cfrp_pores.dispersion import Plate,project_mode


class PoreTests(unittest.TestCase):
    def test_rank_scaling_and_exact_fraction(self):
        f=np.random.default_rng(4).normal(size=(200,100))
        original=ranked_mask(f,.015)
        np.testing.assert_array_equal(original,ranked_mask(f*37,.015))
        self.assertEqual(original.sum(),300)
        self.assertFalse(ranked_mask(f,0).any())

    def test_gaussian_correlation_length(self):
        spec=PoreSpec((600,600),2,(30,12),.015,22)
        f=random_field(spec)
        rows=axial_correlations(f,ranked_mask(f,.015),2,25)
        for axis,a in enumerate(spec.correlation_um):
            actual=[r['correlation'] for r in rows if r['kind']=='field' and r['axis']==axis and r['lag_um']==a][0]
            self.assertAlmostEqual(actual,np.exp(-1),delta=.12)
        np.testing.assert_array_equal(f,random_field(spec))

    def test_geometry_area_and_islands(self):
        spec=PoreSpec((200,200),2,(20,10),.03)
        polygons,q=geometry_from_field(random_field(spec),spec)
        self.assertAlmostEqual(sum(p.area for p in polygons)/40000,.03,delta=1e-4)
        x,y=np.meshgrid(np.arange(101),np.arange(101),indexing='ij')
        field=(np.sqrt((x-50)**2+(y-50)**2)-25)**2
        rings=contours_from_field(field,25,1,(101,101))
        self.assertEqual(len(rings),1);self.assertEqual(len(rings[0].interiors),1)

    def test_reposition_keeps_shapes_and_no_overlap(self):
        pores=[Point(0,0).buffer(r,quad_segs=12) for r in [5,7,8,9,10]]
        for cluster in [False,True]:
            moved=reposition(pores,(200,200),9,cluster)
            np.testing.assert_allclose(sorted(p.area for p in pores),sorted(p.area for p in moved),rtol=1e-12)
            for i,a in enumerate(moved):
                for b in moved[:i]:self.assertFalse(a.intersects(b))

    def test_3d_is_not_an_extrusion(self):
        spec=PoreSpec((120,100,80),2,(16,12,8),.03,9)
        field=random_field(spec); mask=ranked_mask(field,.03)
        self.assertFalse(np.array_equal(mask[:,:,10],mask[:,:,20]))
        self.assertAlmostEqual(mask.mean(),.03,delta=1/field.size)
        rows,meta=mask_statistics(mask,spec);self.assertGreater(meta['pore_count'],0)

    def test_stl_is_closed_and_outward(self):
        mask=np.zeros((8,8,8),bool);mask[2:6,2:6,2:6]=True
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'pore.stl';volume=export_stl(path,mask,2)
            v=np.array([[float(x) for x in line.split()[1:]] for line in path.read_text().splitlines() if line.startswith('vertex ')]).reshape(-1,3,3)
            signed=np.einsum('ij,ij->i',v[:,0],np.cross(v[:,1],v[:,2])).sum()/6
            self.assertGreater(signed,0);self.assertAlmostEqual(signed/volume,1,places=8)
            normal=np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0])
            self.assertTrue(np.all(np.einsum('ij,ij->i',normal,v.mean(axis=1)-8e-6)>0))

    def test_gaussian_surface_volume_preserves_porosity(self):
        spec=PoreSpec((100,100,100),5,(30,20,30),.03,81)
        with tempfile.TemporaryDirectory() as folder:
            volume,q=export_gaussian_stl(Path(folder)/'pores.stl',random_field(spec),spec)
            self.assertAlmostEqual(volume/(100**3*1e-18),.03,delta=1e-4)


class AttenuationTests(unittest.TestCase):
    def test_known_exponential_and_reference(self):
        x=np.array([.01,.02,.03,.04]);ref=2*np.exp((-3+400j)*x)
        a=ref*np.exp((-7+10j)*x)
        fitted=fit_attenuation(x,a,ref)
        self.assertAlmostEqual(fitted['alpha_np_m'],7,places=10)
        self.assertAlmostEqual(fitted['alpha_db_m'],7*8.685889638,places=6)
        self.assertAlmostEqual(fitted['phase_slope_rad_m'],10,places=8)

    def test_reject_invalid_receivers(self):
        for x,a in [([0,1],[1,1]),([0,0,1],[1,1,1]),([0,1,2],[1,0,1])]:
            with self.assertRaises(ValueError):fit_attenuation(x,a)

    def test_zero_and_statistics(self):
        self.assertTrue(converged_alpha(0,0,.1)['passed'])
        self.assertFalse(converged_alpha(10,20,.1)['passed'])
        small=converged_alpha(8,2,.0008)
        self.assertTrue(small['passed'])
        self.assertTrue(small['absolute_cumulative_passed'])
        self.assertFalse(small['relative_coefficient_passed'])
        ci=uncertainty([1,2,3,4]);self.assertLess(ci['ci_lower'],2.5);self.assertGreater(ci['ci_upper'],2.5)

    def test_complex_average_is_not_magnitude_average(self):
        np.testing.assert_allclose(coherent_mean([[1,1],[-1,-1]]),[0,0])

    def test_counterpropagating_wave_separation(self):
        x=np.linspace(0,.003,30);f=5e6;Z=4e6;k=13000.
        forward=(2+3j)*np.exp(-1j*k*x);back=(1-2j)*np.exp(1j*k*x)
        stress=2j*np.pi*f*Z*(-forward+back)
        np.testing.assert_allclose(forward_bulk_displacement(forward+back,stress,f,Z),forward,atol=1e-12)

    def test_gates_preserve_time_phase(self):
        t=np.linspace(0,.01,20001);f=2000.;tau=[.003,.003125]
        signals=np.array([np.exp(-((t-q)/.0003)**2)*np.cos(2*np.pi*f*(t-q)) for q in tau])
        z=gated_spectra(t,signals,[(q-.001,q+.001) for q in tau],f)
        self.assertAlmostEqual(np.angle(z[1]/z[0]),-np.pi/2,places=5)


class DispersionTests(unittest.TestCase):
    def test_isotropic_low_frequency_limits_and_projection(self):
        E=70e9;nu=.3;rho=2700.;mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu))
        C=np.array([[lam+2*mu,lam,0],[lam,lam+2*mu,0],[0,0,mu]])
        p=Plate(.001,C,rho,60)
        s=p.at_frequency(10000,'S0');a=p.at_frequency(10000,'A0')
        extensional=np.sqrt(E/(rho*(1-nu**2)))
        self.assertAlmostEqual(s['phase_velocity_m_s']/extensional,1,delta=.005)
        bending=E*.001**3/(12*(1-nu**2))
        thin_k=(rho*.001*(2*np.pi*10000)**2/bending)**.25
        self.assertAlmostEqual(a['k_rad_m']/thin_k,1,delta=.04)
        self.assertGreater(s['parity'],.99);self.assertLess(a['parity'],-.99)
        amplitude=2+3j
        self.assertAlmostEqual(abs(project_mode(s['y_m'],amplitude*s['u'],amplitude*s['v'],s)-amplitude),0,places=8)


if __name__=='__main__':unittest.main()
