"""
This file is part of the HoneyComb project.
Copyright 2015 David W. Hogg (NYU).
"""
import numpy as np
import matplotlib.pyplot as pl

def make_one_signal(times, exptimes, A, B, omega):
    """
    `make_one_signal`

    Check the analytic integration and trig identities!

    Note:  We don't integrate, we *average*.
    """
    sineterm = 2. * np.sin(0.5 * omega * exptimes) / omega
    meanangle = omega * times
    return (A * np.cos(meanangle) * sineterm +
            B * np.sin(meanangle) * sineterm) / exptimes

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
        t2s = np.cumsum(np.random.randint(0.75 * n0, high = 1.25 * n0, size=ntimes) * dt)
    else:
        t2s = np.cumsum(np.ones(ntimes) * n0 * dt)
    t1s = np.zeros_like(t2s)
    t1s[1:] = t2s[:-1]
    times = 0.5 * (t2s + t1s)
    exptimes = (t2s - t1s)
    print "make_fake_data(): made", n0, "exposures"
    print times[-5:]
    print exptimes[-5:]
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

def plot_exptimes(times, exptimes, fluxes, prefix, title=None):
    """
    `plot_exptimes`

    Make a histogram of the exposure times we actually got.
    """
    pl.clf()
    pl.subplot(2, 1, 1)
    range = (0.001, 3600.001)
    pl.hist(exptimes, bins=100, range=range, histtype="step", color="k")
    pl.xlim(range)
    pl.xlabel("exposure time (s)")
    pl.ylabel("number")
    if title is not None:
        pl.title(title)
    pl.subplot(2, 1, 2)
    pl.plot(times, fluxes, "k.", alpha=0.5)
    pl.xlabel("time (s)")
    pl.ylabel("flux")
    pl.savefig(prefix+".png")
    return None

def analyze_fake_data(times, exptimes, fluxes, ivars):
    """
    `analyze fake data`

    Somehow figure out the signal-to-noise or information content of
    various frequency signals in the data.

    Cramer-Rao bound to start?
    """
    return None

if __name__ == "__main__":
    times1, exptimes1, fluxes1, ivars1 = make_fake_data()
    plot_exptimes(times1, exptimes1, fluxes1, "hogg", "Hogg proposal")
    times2, exptimes2, fluxes2, ivars2 = make_fake_data(random=False)
    plot_exptimes(times2, exptimes2, fluxes2, "TESS", "TESS default")
    analyze_fake_data(times1, exptimes1, fluxes1, ivars1)
