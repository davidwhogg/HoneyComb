"""
This file is part of the HoneyComb project.
Copyright 2015 David W. Hogg (NYU).
"""
import numpy as np
import matplotlib.pyplot as pl
import os
import cPickle as pickle
from functools import partial

# multiprocessing trix DON'T WORK
from multiprocessing import Pool
nthreads = 1
if nthreads == 1:
    pmap = map
else:
    pmap = Pool(processes=nthreads).map

# globals
np.random.seed(42)
ntimes = 1500 # number of data points
dt = 2. # exposure time of sub-exposures (s)
n0 = 900 # number of sub-exposures per exposure for S.O.P.
period = 5.55 * 60. # solar period (s)

def write_to_pickle(fn, stuff):
    filehandler = open(fn, "wb")
    pickle.dump(stuff, filehandler)
    filehandler.close()

def read_from_pickle(fn):
    filehandler = open(fn, "rb")
    stuff = pickle.load(filehandler)
    filehandler.close()
    return stuff

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
    # make time vectors
    if random:
        t2s = np.cumsum(np.random.randint(0.666667 * n0, high = 1.333333 * n0, size=ntimes) * dt)
    else:
        t2s = np.cumsum(np.ones(ntimes) * n0 * dt)
    t1s = np.zeros_like(t2s)
    t1s[1:] = t2s[:-1]
    times = 0.5 * (t2s + t1s)
    exptimes = (t2s - t1s)
    print "make_fake_data(): made", len(times), "exposures"
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
    pl.xlim((15, 17))
    pl.xlabel("time (d)")
    pl.ylim(0.997, 1.003)
    pl.ylabel("flux")
    pl.savefig(prefix+".png")
    return None

class Compute_one_crb:
    """
    `Compute_one_crb`

    Compute, for one period in the `periods` list, the Cramer-Rao
    bound on the sin() and cos() terms.

    It's a class to permit pickling.

    Idiotically slow.
    """
    def __init__(self, times, exptimes, ivars):
        self.times = times
        self.exptimes = exptimes
        self.ivars = ivars

    def __call__(self, P):
        omega = 2. * np.pi / P
        modelA = make_one_signal(self.times, self.exptimes, 1., 0., omega)
        modelB = make_one_signal(self.times, self.exptimes, 0., 1., omega)
        return np.dot(modelA, self.ivars * modelA), np.dot(modelB, self.ivars * modelB)

def compute_crbs(periods, times, exptimes, ivars):
    """
    `compute_crbs`

    Just a `map()`.
    """
    foo = Compute_one_crb(times, exptimes, ivars)
    return np.array(pmap(foo, periods))

class Compute_one_alias:
    """
    `Compute_one_alias`

    Compute, for one period in the `periods` list, the "cosine
    distance" between the given period onto all other periods.

    It's a class to permit pickling.

    Idiotically slow.
    """
    def __init__(self, period0, times, exptimes, ivars):
        self.times = times
        self.exptimes = exptimes
        self.ivars = ivars
        omega0 = 2. * np.pi / period0
        mA = make_one_signal(times, exptimes, 1., 0., omega0)
        mB = make_one_signal(times, exptimes, 0., 1., omega0)
        self.mA = mA / np.sqrt(np.dot(mA, ivars * mA))
        self.mB = mB / np.sqrt(np.dot(mB, ivars * mB))

    def __call__(self, P):
        omega = 2. * np.pi / P
        modelA = make_one_signal(self.times, self.exptimes, 1., 0., omega)
        modelB = make_one_signal(self.times, self.exptimes, 0., 1., omega)
        ivmodelA = self.ivars * modelA
        ivmodelB = self.ivars * modelB
        normA2 = np.dot(modelA, ivmodelA)
        normB2 = np.dot(modelB, ivmodelB)
        return (np.sqrt(np.dot(self.mA, ivmodelA) ** 2 / normA2 +
                        np.dot(self.mA, ivmodelB) ** 2 / normB2),
                np.sqrt(np.dot(self.mB, ivmodelA) ** 2 / normA2 +
                        np.dot(self.mB, ivmodelB) ** 2 / normB2))

def compute_aliases(period0, periods, times, exptimes, ivars):
    """
    `compute_aliases`

    Just a `map()`.
    """
    foo = Compute_one_alias(period0, times, exptimes, ivars)
    return np.array(pmap(foo, periods))

def plot_stuff(periods, crbs1, crbs2, name1, name2, ylabel, prefix, logy):
    """
    `plot_stuff`
    """
    vline_periods = np.array([2. * n0 * dt,
                              n0 * dt,
                              2. * dt,
                              dt])
    for ii in range(2):
        if ii == 0:
            x0 = period
            xs = periods
            vxs = vline_periods
        if ii == 1:
            x0 = 1. / period
            xs = 1. / periods
            vxs = 1. / vline_periods
        pl.clf()
        pl.plot(xs, crbs1, "g-", alpha=0.5, label=name1)
        pl.plot(xs, crbs2, "k-", alpha=0.25, label=name2)
        for xx in vxs:
            pl.axvline(xx, color="r", alpha=0.5)
        # pl.axvline(x0, color="b", alpha=0.5)
        pl.legend(loc=2)
        if ii == 0:
            if logy:
                pl.loglog()
            else:
                pl.semilogx()
            pl.xlabel("period (s)")
        if ii == 1:
            if logy:
                pl.semilogy()
            pl.xlabel("frequency (Hz)")
        pl.ylabel(ylabel)
        big = max(np.max(crbs1), np.max(crbs2))
        if logy:
            pl.ylim(1.e-8 * big, 1.e1 * big)
        else:
            pl.ylim(0., 1.)
        pl.xlim(min(xs), max(xs))
        pl.savefig("%s%1d.png" % (prefix, ii))
    return None

def plot_crbs(periods, crbs1, crbs2, name1, name2):
    return plot_stuff(periods,
                      np.sum(crbs1, axis=1),
                      np.sum(crbs2, axis=1),
                      name1, name2, "Cramer-Rao bounds", "crb", True)

def plot_aliases(periods, aliases1, aliases2, name1, name2):
    return plot_stuff(periods,
                      np.max(aliases1, axis=1),
                      np.max(aliases2, axis=1),
                      name1, name2, "Worst cosine distance", "alias", False)

if __name__ == "__main__":
    picklefn = "./exptime.pkl"
    if not os.path.exists(picklefn):
        dlnp = 0.05 / ntimes # many of our best people died to bring us MAGIC number "0.05"
        periods = np.exp(np.arange(np.log(3.0e0) + 0.5 * dlnp, np.log(3.0e4), dlnp)) # (s)
        print "periods", len(periods)
        times1, exptimes1, fluxes1, ivars1 = make_fake_data()
        aliases1 = compute_aliases(period, periods, times1, exptimes1, ivars1)
        crbs1 = compute_crbs(periods, times1, exptimes1, ivars1)
        times2, exptimes2, fluxes2, ivars2 = make_fake_data(random=False)
        crbs2 = compute_crbs(periods, times2, exptimes2, ivars2)
        aliases2 = compute_aliases(period, periods, times2, exptimes2, ivars2)
        write_to_pickle(picklefn, (periods,
                                   times1, exptimes1, fluxes1, ivars1, crbs1, aliases1,
                                   times2, exptimes2, fluxes2, ivars2, crbs2, aliases2))
    else:
        (periods,
         times1, exptimes1, fluxes1, ivars1, crbs1, aliases1,
         times2, exptimes2, fluxes2, ivars2, crbs2, aliases2) = read_from_pickle(picklefn)
    plot_exptimes(times1, exptimes1, fluxes1, "hogg", "Hogg proposal")
    plot_exptimes(times2, exptimes2, fluxes2, "TESS", "TESS default")
    plot_crbs(periods, crbs1, crbs2, "Hogg proposal", "TESS default")
    plot_aliases(periods, aliases1, aliases2, "Hogg proposal", "TESS default")
