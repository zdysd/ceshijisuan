import unittest
import numpy as np
from guided_collective import dephase_average,grouped_pores
from guided_statistics import coherent_interval,attenuation_interval


class GuidedCollectiveTests(unittest.TestCase):
    def test_incident_ports_keep_phase_reference(self):
        x=np.array([.01,.012,.014]);k=1300.;a=2+3j
        self.assertAlmostEqual(abs(dephase_average(a*np.exp(-1j*k*x),x,k)-a),0,places=12)

    def test_group_area_and_internal_cavities(self):
        counts=[]
        for d in [20,40,80]:
            p,m=grouped_pores(d,1,101)
            self.assertAlmostEqual(m['actual_porosity'],.015,places=12)
            self.assertTrue(all(min(q.bounds[:2])>=49.99999 and max(q.bounds[2:])<=950.00001 for q in p))
            self.assertTrue(all(not a.intersects(b) for i,a in enumerate(p) for b in p[:i]))
            counts.append(len(p))
        self.assertGreater(counts[0],counts[1]);self.assertGreater(counts[1],counts[2])

    def test_length_fit_retains_interfaces_and_phase(self):
        x=np.array([.004,.008,.016]);alpha=2.5;phase=40.
        a=np.array([(.97+z*1e-4)*np.exp(-(alpha+1j*phase)*x) for z in [-1,0,1,2]])
        result=attenuation_interval(x,a,100)
        self.assertAlmostEqual(result['fit']['alpha_np_m'],alpha,places=10)
        self.assertTrue(result['positive_attenuation_identified'])
        self.assertTrue(result['exponential_fit_passed'])

    def test_coherent_loss_can_include_phase_decorrelation(self):
        a=np.exp(1j*np.array([-.2,-.1,.1,.2]))
        r=coherent_interval(a,100)
        self.assertGreater(r['amplitude_reduction'],0)
        self.assertAlmostEqual(r['mean_individual_amplitude_reduction'],0,places=12)
        self.assertAlmostEqual(r['individual_magnitude_standard_deviation'],0,places=12)


if __name__=='__main__':unittest.main()
