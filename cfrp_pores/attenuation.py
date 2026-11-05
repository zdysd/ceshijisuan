"""Amplitude attenuation, reference subtraction and sample uncertainty."""
import numpy as np
from scipy import stats


def fit_attenuation(distance_m, amplitudes, reference=None):
    """Fit spectral/modal amplitudes at >=3 distinct receiver positions.

    Supply a reference at the SAME receivers for excess attenuation.
    Amplitudes can be complex. Distances are metres; output is Np/m, dB/m.
    Negative fitted slopes are reported rather than clipped to a desired result.
    """
    x=np.asarray(distance_m,float)
    a=np.asarray(amplitudes)
    if x.ndim != 1 or a.shape != x.shape or len(x)<3 or len(np.unique(x))!=len(x):
        raise ValueError('At least three distinct positions with matching amplitudes required.')
    if not np.isfinite(x).all() or not np.isfinite(a).all() or np.any(np.abs(a)<=0):
        raise ValueError('Positions and amplitudes must be finite; amplitudes must be nonzero.')
    if reference is not None:
        ref=np.asarray(reference)
        if ref.shape!=a.shape or not np.isfinite(ref).all() or np.any(np.abs(ref)<=0):
            raise ValueError('Reference must match positions and have finite nonzero amplitudes.')
        a=a/ref
    logs=np.log(np.abs(a))
    if np.ptp(logs)<1e-14:
        return dict(alpha_np_m=0.,alpha_db_m=0.,alpha_db_mm=0.,slope_stderr_np_m=0.,r_squared=1.,
                    log_amplitude_intercept=float(logs.mean()),
                    log_amplitude_residual_max=0.,phase_slope_rad_m=float(stats.linregress(x,np.unwrap(np.angle(a))).slope),
                    distance_span_m=float(np.ptp(x)),reference_subtracted=reference is not None)
    result=stats.linregress(x,logs)
    prediction=result.intercept+result.slope*x
    phase=stats.linregress(x,np.unwrap(np.angle(a)))
    alpha=-result.slope
    return dict(alpha_np_m=float(alpha),alpha_db_m=float(alpha*20/np.log(10)),
                log_amplitude_intercept=float(result.intercept),
                alpha_db_mm=float(alpha*20/np.log(10)/1000),
                slope_stderr_np_m=float(result.stderr),r_squared=float(result.rvalue**2),
                log_amplitude_residual_max=float(np.max(np.abs(np.log(np.abs(a))-prediction))),
                phase_slope_rad_m=float(phase.slope),distance_span_m=float(np.ptp(x)),
                reference_subtracted=reference is not None)


def uncertainty(values,confidence=.95):
    values=np.asarray(values,float)
    if len(values)<2 or not np.isfinite(values).all():
        raise ValueError('At least two independent finite samples are required for a confidence interval.')
    mean=float(values.mean()); sd=float(values.std(ddof=1))
    half=float(stats.t.ppf((1+confidence)/2,len(values)-1)*sd/np.sqrt(len(values)))
    return dict(n=len(values),mean=mean,standard_deviation=sd,
                confidence=confidence,ci_lower=mean-half,ci_upper=mean+half)


def converged_alpha(first,second,length_m,relative=.05,near_zero_np=.01):
    # Near-zero is defined by cumulative attenuation, avoiding division by zero.
    diff=abs(first-second)
    near_zero=max(abs(first),abs(second))*length_m <= near_zero_np
    relative_change=float(diff/max(abs(first),abs(second))) if max(abs(first),abs(second)) else 0.
    relative_passed=relative_change<=relative
    absolute_passed=diff*length_m<=near_zero_np
    passed=absolute_passed if near_zero else relative_passed
    return dict(passed=bool(passed),near_zero=bool(near_zero),delta_cumulative_np=float(diff*length_m),
                relative_change=relative_change,relative_coefficient_passed=bool(relative_passed),
                absolute_cumulative_passed=bool(absolute_passed))


def gated_spectra(times_s, signals, gates_s, frequency_hz):
    """Fourier integral at an arbitrary frequency with common absolute time.

    Each receiver gate must contain its complete first-arrival pulse, excluding
    echoes. Do not reset each pulse to t=0: that destroys propagation phase.
    """
    times=np.asarray(times_s,float); data=np.asarray(signals)
    if data.shape!=(len(gates_s),len(times)) or np.any(np.diff(times)<=0):
        raise ValueError('Signals must have shape receivers x times and strictly increasing times.')
    result=[]
    for signal,(start,end) in zip(data,gates_s):
        selection=(times>=start)&(times<=end)
        if selection.sum()<8: raise ValueError('Arrival gate needs at least eight samples.')
        t=times[selection]
        result.append(np.trapezoid(signal[selection]*np.exp(-2j*np.pi*frequency_hz*t),t))
    return np.array(result)


def coherent_mean(amplitudes):
    data=np.asarray(amplitudes)
    if data.ndim!=2 or not np.isfinite(data).all():
        raise ValueError('Expected independent samples x receivers with common phase reference.')
    return data.mean(axis=0)  # mean complex signal, NOT mean magnitude


def forward_bulk_displacement(displacement,normal_stress,frequency_hz,impedance):
    """Separate +/-x waves for exp(+iwt), using a reference P-wave impedance.

    sigma_xx=-i*w*Z*u for a +x wave, +i*w*Z*u for a -x wave.
    In porous media this is a reference-medium projection, whose residual
    backward component must still be checked; it is not an exact tensor fit.
    """
    if frequency_hz<=0 or impedance<=0:raise ValueError('Positive frequency and impedance required.')
    return .5*(np.asarray(displacement)-np.asarray(normal_stress)/(2j*np.pi*frequency_hz*impedance))


def fit_quality(fit,near_zero_np=.01):
    cumulative=abs(fit['alpha_np_m'])*fit['distance_span_m']
    return dict(cumulative_attenuation_np=float(cumulative),
                exponential_fit_accepted=bool(fit['log_amplitude_residual_max']<=near_zero_np and
                                             (cumulative<=near_zero_np or fit['r_squared']>=.95)),
                near_zero=cumulative<=near_zero_np,
                note='A solved model or fitted slope is not sufficient evidence of a stable material attenuation coefficient.')


def homogeneous_propagation(x_m,k_real_rad_m,alpha_np_m,amplitude=1.):
    if alpha_np_m<0: raise ValueError('Passive equivalent propagation requires nonnegative attenuation.')
    return amplitude*np.exp((1j*k_real_rad_m-alpha_np_m)*np.asarray(x_m))
