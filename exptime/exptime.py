"""
This file is part of the HoneyComb project.
Copyright 2015 David W. Hogg (NYU).
"""
import numpy as np

def make_one_signal(times, exptimes, A, B, omega):
    """
    `make_one_signal`

    Check the analytic integration!
    """
    t1s = times - 0.5 * exptimes
    t2s = times + 0.5 * exptimes
    return ((A / omega) * np.sin(omega * t2s) - A * np.sin(omega * t1s) -
            (B / omega) * np.cos(omega * t2s) - A * np.cos(omega * t1s))

def make_fake_data(random=True):
    """
    `make_fake_data`

    Put an integrated signal with a few sinusioidal frequencies into
    some noisy data.
    """
    # set parameters
    ntimes = 70000
    dt = 2. # sec
    n0 = 900.
    # make time vectors
    if random:
        t2s = np.cumsum(np.random.randint(0.5 * n0, high = 1.5 * n0, size=ntimes) * dt)
    else:
        t2s = np.cumsum(np.ones(ntimes) * n0 * dt)
    t1s = np.zeros_like(t2s)
    t1s[1:] = t2s[:-1]
    times = 0.5 * (t2s + t1s)
    exptimes = (t2s - t1s)
    print times[-5:], exptimes[-5:]
    # make ivars with proper exptime variation; make noise
    ivars = 1. * exptimes
    fluxes = np.random.normal(size=ntimes) / np.sqrt(ivars)
    # make signals
    omega = 0.0189 # rad / s - peak frequency in angular units
    delta_omega = 0.0003 # rad / s - large frequency difference in angular units
    fluxes += make_one_signal(times, exptimes, 0.1, 0.2, omega - 2. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.2, 0.2, omega - 1. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.2, 0.3, omega + 0. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.2, 0.2, omega + 1. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.2, 0.1, omega + 2. * delta_omega)
    return times, exptimes, fluxes, ivars

def analyze_fake_data(times, exptimes, fluxes, ivars):
    """
    `analyze fake data`

    Somehow figure out the signal-to-noise or information content of
    various frequency signals in the data.

    Likelihood ratios to start.

    This should use the superpgram; if it doesn't it is because I
    suck.
    """
    return None

if __name__ == "__main__":
    times1, exptimes1, fluxes1, ivars1 = make_fake_data()
    times2, exptimes2, fluxes2, ivars2 = make_fake_data(random=False)
    analyze_data(times1, exptimes1, fluxes1, ivars1)
