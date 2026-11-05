import unittest
import numpy as np
from cfrp_pores.transient import hann_pulse,source_band,fourier_history,scattering,modal_bank,analyze_ports,tail_fraction,field_frame_indices,minimum_inplane_speed,continuation_block_ends,bulk_adjustment
from cfrp_pores.dispersion import Plate
from cfrp_pores.attenuation import fit_attenuation


class TransientTests(unittest.TestCase):
    def test_bulk_refinement_prioritizes_incomplete_windows(self):
        checks=[dict(kind='bulk',case=c,complete=True,passed=False,differences={},frequency_band_consistency=dict(passed=False)) for c in ['reference','group']]
        gate=dict(checks=checks,thresholds=dict(tail_energy=1e-4,absolute_amplitude=2e-4,phase_degrees=1.,relative_nonzero_loss=.05))
        rows=[dict(kind='bulk',case='reference',variant='base',tail_energy_fraction=1.25e-4,window_delta_max=9.92e-4)]
        self.assertEqual(bulk_adjustment(gate,rows),dict(action='extend_window',cases=[['base','reference']]))
        rows[0].update(tail_energy_fraction=1e-8,window_delta_max=1e-5)
        self.assertEqual(bulk_adjustment(gate,rows)['action'],'refine_time')
        checks[0]['complete']=False
        with self.assertRaises(ValueError):bulk_adjustment(gate,rows)

    def test_bulk_refinement_does_not_hide_unknown_failures(self):
        checks=[dict(kind='bulk',case=c,complete=True,passed=False,differences={},frequency_band_consistency=dict(passed=True)) for c in ['reference','group']]
        gate=dict(checks=checks,thresholds=dict(tail_energy=1e-4,absolute_amplitude=2e-4,phase_degrees=1.,relative_nonzero_loss=.05))
        with self.assertRaises(RuntimeError):bulk_adjustment(gate,[])
        for c in checks:c['passed']=True
        self.assertEqual(bulk_adjustment(gate,[])['action'],'accepted')

    def test_extension_keeps_former_final_short_block(self):
        self.assertEqual(continuation_block_ends([30,60,90,120,150,180,210,240,270,280],[30,60,90,120,140]),
                         [30,60,90,120,140,150,180,210,240,270,280])
        self.assertEqual(continuation_block_ends([30,60,90,120,140],[30,60]),[30,60,90,120,140])
        with self.assertRaises(ValueError):continuation_block_ends([30,60],[30,70])

    def test_shortest_anisotropic_wave_uses_acoustic_tensor(self):
        rho=1815.;C=np.array([[10.478,3.029,0],[3.029,10.478,0],[0,0,7.170]])*1e9
        self.assertAlmostEqual(minimum_inplane_speed(C,rho),np.sqrt((C[0,0]-C[0,1])/(2*rho)),places=7)
        lam,mu=3e9,2e9
        isotropic=np.array([[lam+2*mu,lam,0],[lam,lam+2*mu,0],[0,0,mu]])
        self.assertAlmostEqual(minimum_inplane_speed(isotropic,rho),np.sqrt(mu/rho),places=7)

    def test_animation_cadence_does_not_alias_carrier(self):
        for fmax,duration in [(6718750.,13e-6),(268750.,140e-6)]:
            times=np.arange(0,duration,1/(20*fmax));indices=field_frame_indices(times,fmax)
            self.assertEqual(indices[0],0);self.assertEqual(indices[-1],len(times)-1)
            self.assertLessEqual(max(np.diff(times[indices]))*fmax,1/8+1e-12)

    def test_pulse_band_covers_both_quadratures(self):
        low,high=source_band(200000.)
        self.assertLess(low,200000);self.assertGreater(high,200000)
        self.assertLess(high,400000)
        np.testing.assert_array_equal(hann_pulse([-1,0,25e-6,1],200000),[0,0,0,0])

    def test_absolute_phase_survives_gate(self):
        f=100.;t=np.arange(0,1,.00005);delay=.023
        pulse=hann_pulse(t-delay,f)
        spec=fourier_history(t,pulse[:,None],[f])[0,0]
        original=fourier_history(t,hann_pulse(t,f)[:,None],[f])[0,0]
        np.testing.assert_allclose(spec/original,np.exp(-2j*np.pi*f*delay),atol=1e-10)

    def test_mixed_incidence_matrix_recovers_conversion(self):
        incoming=np.array([[2,.2j],[.1,1.5]])
        expected=np.array([[.8+.1j,.02j],[.04,.7-.2j]])
        reflection=np.array([[.1,.02],[.03,.15]])
        T,R,right,_=scattering(incoming,expected@incoming,reflection@incoming,np.zeros((2,2)))
        np.testing.assert_allclose(T,expected);np.testing.assert_allclose(R,reflection);self.assertEqual(abs(right).max(),0)
        with self.assertRaises(ValueError):scattering(np.ones((2,2)),expected,reflection,reflection)

    def test_plane_wave_transmission_keeps_direction_and_reference(self):
        rho=1000.;speed=2000.;C=np.eye(6)*rho*speed**2
        f=100.;t=np.linspace(0,.1,4001);x=np.array([.1,.2,.3,.4,.5,.6]);y=np.linspace(0,.001,17)
        k=2*np.pi*f/speed;Z=rho*speed;values=np.zeros((len(t),6,4,len(y)))
        for p,xx in enumerate(x):
            phase=np.exp(-1j*k*xx);amplitude=phase*(1 if p<3 else .8)
            values[:,p,0,:]=np.real(amplitude*np.exp(2j*np.pi*f*t))[:,None]
            values[:,p,2,:]=np.real(-2j*np.pi*f*Z*amplitude*np.exp(2j*np.pi*f*t))[:,None]
        h=dict(times=t,positions=x,y=y,values=values)
        result=analyze_ports([h],[f],'bulk',(C,rho))
        np.testing.assert_allclose(result['T'][0],[[.8]],atol=1e-12)
        self.assertLess(abs(result['R']).max(),1e-12)

    def test_guided_synthetic_fields_recover_two_modes(self):
        rho=1800.;C3=np.array([[10e9,3e9,0],[3e9,10e9,0],[0,0,7e9]])
        plate=Plate(.001,C3,rho,20);f=200000.;y=np.linspace(1e-8,.001-1e-8,33)
        bank=modal_bank(plate,[f],y);basis,labels,scale,modes,power=bank[0]
        t=np.linspace(0,20/f,1601);positions=np.array([.003,.004,.005,.015,.016,.017])
        incoming=np.array([[1,.07j],[.03,.8]]);expected=np.array([[.9,.02],[.04,.85]])
        histories=[]
        for launch in range(2):
            values=np.zeros((len(t),6,4,len(y)))
            for p,x in enumerate(positions):
                amps=incoming[:,launch] if p<3 else (expected@incoming)[:,launch]
                field=sum(basis[:,2*m]*amps[m]*np.exp(-1j*modes[m]['k_rad_m']*x) for m in range(2)).reshape(4,len(y))
                field[2:]*=scale
                values[:,p]=np.real(np.exp(2j*np.pi*f*t)[:,None,None]*field)
            histories.append(dict(times=t,values=values,y=y,positions=positions))
        result=analyze_ports(histories,[f],'guided',(np.eye(6),rho),bank)
        np.testing.assert_allclose(result['T'][0],expected,atol=1e-10)
        self.assertLess(result['residual'].max(),1e-10)

    def test_tail_identifies_truncated_pulse(self):
        t=np.linspace(0,.1,1001);v=np.zeros((len(t),6,4,17));v[-30:,:,:2]=1
        self.assertGreater(tail_fraction(dict(times=t,values=v)),.9)

    def test_interface_loss_does_not_become_bulk_attenuation(self):
        lengths=np.array([.008,.024,.048])
        result=fit_attenuation(lengths,.8*np.exp(-2.5*lengths+3j*lengths))
        self.assertAlmostEqual(result['alpha_np_m'],2.5,places=11)
        self.assertAlmostEqual(result['log_amplitude_intercept'],np.log(.8),places=11)


if __name__=='__main__':unittest.main()
