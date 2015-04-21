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
    n0 = 900
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
    ivars = 1.e10 * (exptimes / 1800.) # 1e10 at exptime = 30 min
    fluxes = np.random.normal(size=ntimes) / np.sqrt(ivars)
    # make signals
    omega = 0.0189 # rad / s - peak frequency in angular units
    delta_omega = 0.0003 # rad / s - large frequency difference in angular units
    fluxes += 1.
    fluxes += make_one_signal(times, exptimes, 0.0001, 0.0002, omega - 2. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.0002, 0.0002, omega - 1. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.0002, 0.0003, omega + 0. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.0002, 0.0002, omega + 1. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.0002, 0.0001, omega + 2. * delta_omega)
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
    pl.plot(times / 86400., fluxes, "k.")
    pl.xlim((135, 137))
    pl.xlabel("time (d)")
    pl.ylim(0.997, 1.003)
    pl.ylabel("flux")
    pl.savefig(prefix+".png")
    return None

def compute_crbs(periods, times, exptimes, ivars):
    """
    `compute_crbs`

    Compute, for each period in the `periods` list, the Cramer-Rao
    bounds on the sin() and cos() terms.

    Idiotically slow.
    """
    nperiod = len(periods)
    crbs = np.zeros((nperiod, 2))
    for ii, period in enumerate(periods):
        omega = 2 * np.pi / period
        modelA = make_one_signal(times, exptimes, 1., 0., omega)
        modelB = make_one_signal(times, exptimes, 0., 1., omega)
        crbs[ii, 0] = np.dot(modelA, ivars * modelA)
        crbs[ii, 1] = np.dot(modelB, ivars * modelB)
    return crbs

def plot_crbs(periods, crbs1, crbs2, name1=None, name2=None):
    """
    `plot_exptimes`

    Make a histogram of the exposure times we actually got.
    """
    pl.clf()
    pl.plot(periods, crbs1[:,0], "k-")
    pl.plot(periods, crbs1[:,1], "k-")
    pl.plot(periods, crbs2[:,0], "k-", alpha=0.5)
    pl.plot(periods, crbs2[:,1], "k-", alpha=0.5)
    pl.loglog()
    pl.xlabel("period (s)")
    big = max(np.max(crbs1), np.max(crbs2))
    pl.ylim(1.e-8 * big, 1.e1 * big)
    pl.ylabel("Cramer-Rao bounds")
    pl.savefig("crb.png")
    return None

if __name__ == "__main__":
    times1, exptimes1, fluxes1, ivars1 = make_fake_data()
    periods = np.exp(np.arange(0., np.log(100000.), 0.05))
    crbs1 = compute_crbs(periods, times1, exptimes1, ivars1)
    plot_exptimes(times1, exptimes1, fluxes1, "hogg", "Hogg proposal")
    times2, exptimes2, fluxes2, ivars2 = make_fake_data(random=False)
    crbs2 = compute_crbs(periods, times2, exptimes2, ivars2)
    plot_exptimes(times2, exptimes2, fluxes2, "TESS", "TESS default")
    plot_crbs(periods, crbs1, crbs2, "Hogg proposal", "TESS default")
