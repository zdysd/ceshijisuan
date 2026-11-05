"""Offline checks of the scientific rules and dependency/resource gates."""
import tempfile,unittest
from pathlib import Path
from unittest.mock import patch,MagicMock
import numpy as np
import detectability_200kHz_L48 as d
import detectability_worker as w
import detectability_followon as follow

class DetectabilityTests(unittest.TestCase):
    def test_incident_mixing_is_removed(self):
        incident=np.array([[1.,.23j],[.17,.8]],complex)
        T=np.array([[.998+.01j,.03],[-.02j,.997-.02j]])
        R=np.array([[.02j,.01],[.004,.03j]])
        launches=[dict(pre_forward=incident[:,i],post_forward=(T@incident)[:,i],pre_backward=(R@incident)[:,i],post_backward=np.zeros(2),powers=[1.,2.],residuals=[0.]) for i in range(2)]
        a=w.assemble(launches,np.zeros(6));np.testing.assert_allclose(a['T'],T,atol=1e-15);np.testing.assert_allclose(a['R'],R,atol=1e-15)
        self.assertGreater(abs((T@incident)[0,0]/incident[0,0]-T[0,0]),1e-4)

    def test_frequency_band_keeps_transient_absorber_and_refreshes_modes(self):
        model=MagicMock();plate=MagicMock();settings=dict(plate=plate,h=.001)
        for frequency,mode in [(131250.,'S0'),(268750.,'A0')]:
            wave=object();plate.at_frequency.return_value=wave
            with patch.object(w.gc,'set_source') as source:
                w.frequency_source(model,settings,frequency,mode)
                model.java.param().set.assert_called_with('f0','200000.0[Hz]')
                model.java.study('std1').feature('freq').set.assert_called_with('plist',str(frequency))
                plate.at_frequency.assert_called_with(frequency,mode)
                source.assert_called_once_with(model,wave,.001)

    def test_signed_loss_and_uncertainty_not_discarded(self):
        self.assertEqual(d.classify(-.002,.0001),'below')
        self.assertEqual(d.classify(.001,.0001),'indeterminate')
        self.assertEqual(d.classify(.002,None),'unverified_uncertainty')
        self.assertEqual(d.classify(.0013,.0002),'above')

    def test_exact_binomial_fixed_cohort(self):
        self.assertAlmostEqual(d.cp_lower(30,30),.05**(1/30),places=12)
        self.assertGreater(d.cp_lower(30,30),.9)
        self.assertLess(d.cp_lower(29,30),.9)
        self.assertLess(d.cp_lower(11,12),.67)

    def test_legacy_single_ports_rejected_for_48mm_reference(self):
        a=d.source_record(d.OLD/'d80/single_maximum');b=d.source_record(d.OLD/'reference_L48')
        self.assertFalse(d.compatible(a,b))
        with self.assertRaises(ValueError):d.metrics(a,b,'bad','single')

    def test_frozen_geometry_and_protocol(self):
        protocol=d.read(d.OUT/'protocol.json')
        self.assertEqual(d.digest(d.OUT/'protocol.json'),d.read(d.OUT/'protocol.sha256.json')['sha256'])
        self.assertEqual(len(set(protocol['independent_seeds'])),30)
        self.assertFalse(set(protocol['independent_seeds']) & set(range(101,113)))
        self.assertEqual(d.digest(d.OUT/'frozen_geometry_manifest.json'),protocol['geometry_manifest_sha256'])
        for v in d.read(d.OUT/'frozen_geometry_manifest.json'):
            self.assertEqual(d.digest(v['path']),v['file_sha256']);self.assertTrue(v['audit']['passed'])
        for v in d.read(d.OUT/'single_controls.json'):
            if not v['usable']:continue
            p,_=d.load_polygons(v['geometry_path']);self.assertAlmostEqual(p[0].centroid.x,24000,places=7)
            self.assertAlmostEqual(p[0].area,v['area_um2'],places=6)
            if v['depth']=='upper':self.assertAlmostEqual(1000-p[0].bounds[3],50,places=7)
            if v['depth']=='lower':self.assertAlmostEqual(p[0].bounds[1],50,places=7)

    def test_protected_resource_gate_precedes_new_session(self):
        parent=d.OUT/'_tests';parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as name:
            folder=Path(name).resolve();self.assertTrue(folder.is_relative_to(d.ROOT.resolve()))
            with patch.object(d,'OUT',folder),patch.object(w,'verify_inputs'),patch.object(w,'wait_resources',side_effect=RuntimeError('protected worker still running')),patch.object(rt_module(),'start_session') as native,patch('detectability_report.render'):
                with self.assertRaisesRegex(RuntimeError,'protected worker'):w.main()
                native.assert_not_called()

    def test_authorized_overlap_keeps_memory_and_shutdown_guards(self):
        parent=d.OUT/'_tests';parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as name:
            folder=Path(name).resolve()
            d.write(folder/'resource_policy.json',dict(allow_concurrent_with_protected_p=True,
                minimum_free_gib_before_session=8,minimum_free_gib_before_solve=4))
            with patch.object(d,'OUT',folder),patch.object(w,'free_physical_gib',return_value=10):
                self.assertTrue(w.resource_decision([dict(pid=123)],initial=True)[0])
                self.assertFalse(w.resource_decision([dict(pid=None)],initial=True)[0])
            with patch.object(d,'OUT',folder),patch.object(w,'free_physical_gib',return_value=6):
                self.assertFalse(w.resource_decision([dict(pid=123)],initial=True)[0])
                self.assertTrue(w.resource_decision([dict(pid=123)],initial=False)[0])

    def test_idle_schedule_limits_overlap_to_p1(self):
        parent=d.OUT/'_tests';parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as name:
            folder=Path(name).resolve()
            d.write(folder/'resource_policy.json',dict(allow_concurrent_with_protected_p=True,
                overlap_stages=['P1_level_0','P1_level_1'],require_user_idle_minutes=20,
                minimum_free_gib_before_session=12,minimum_free_gib_before_solve=8,
                minimum_free_gib_before_solve_with_protected_p=12))
            with patch.object(d,'OUT',folder),patch.object(w,'free_physical_gib',return_value=14):
                with patch.object(w,'user_idle_minutes',return_value=5):
                    self.assertFalse(w.resource_decision([dict(pid=123)],initial=True,stage='P1_level_1')[0])
                with patch.object(w,'user_idle_minutes',return_value=25):
                    self.assertTrue(w.resource_decision([dict(pid=123)],initial=True,stage='P1_level_1')[0])
                    self.assertFalse(w.resource_decision([dict(pid=123)],initial=True,stage='P1_level_2')[0])
                    self.assertFalse(w.resource_decision([dict(pid=123)],initial=True,stage='P2_band')[0])
                    self.assertTrue(w.resource_decision([],initial=True,stage='P2_band')[0])
            with patch.object(d,'OUT',folder),patch.object(w,'free_physical_gib',return_value=10),patch.object(w,'user_idle_minutes',return_value=25):
                self.assertFalse(w.resource_decision([dict(pid=123)],initial=False,stage='P1_level_1')[0])
                self.assertTrue(w.resource_decision([],initial=False,stage='P1_level_1')[0])

    def test_failed_P1_prevents_P2_P3(self):
        parent=d.OUT/'_tests';parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as name:
            folder=Path(name).resolve();self.assertTrue(folder.is_relative_to(d.ROOT.resolve()))
            with patch.object(d,'OUT',folder),patch.object(w,'verify_inputs'),patch.object(w,'wait_resources'),patch.object(rt_module(),'start_session',return_value=object()),patch.object(w,'run_p1',side_effect=RuntimeError('P1 failed')),patch.object(follow,'run') as later,patch('detectability_report.render'):
                with self.assertRaisesRegex(RuntimeError,'P1 failed'):w.main()
                later.assert_not_called()

    def test_guided_time_cap_precedes_native_work(self):
        with patch.object(w,'settings',return_value={}),patch.object(w,'build') as build:
            with self.assertRaisesRegex(RuntimeError,'560 us'):follow.transient_case(None,{},1,False,561e-6)
            build.assert_not_called()

    def test_shutdown_or_incomplete_source_does_not_release_queue(self):
        parent=d.OUT/'_tests';parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as name:
            folder=Path(name).resolve();self.assertTrue(folder.is_relative_to(d.ROOT.resolve()))
            with patch.object(d,'ROOT',folder):
                source=folder/'output/transient_verified';source.mkdir(parents=True)
                d.write(source/'continuation_state.json',dict(stage='waiting_for_bulk_transient_controls'))
                self.assertTrue(w.protected_jobs())
                d.write(source/'continuation_state.json',dict(stage='complete'))
                d.write(source/'PAUSED_FOR_SHUTDOWN.json',dict(stage='ready_for_shutdown'))
                self.assertTrue(w.protected_jobs())
                (source/'PAUSED_FOR_SHUTDOWN.json').unlink()
                self.assertEqual(w.protected_jobs(),[])

def rt_module():return w.rt

if __name__=='__main__':unittest.main()
