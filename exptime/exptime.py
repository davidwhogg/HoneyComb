"""
This file is part of the HoneyComb project.
Copyright 2015 David W. Hogg (NYU).
"""
import os
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as pl
import matplotlib.transforms as transforms
pl.rc("text", usetex=True)
pl.rc("font", family="serif")

# multiprocessing trix DON'T WORK
from multiprocessing import Pool
nthreads = 1
if nthreads == 1:
    pmap = map
else:
    pmap = Pool(processes=nthreads).map

# globals
fontsize = 12
np.random.seed(42)
ntimes = 1500 # number of data points
dt = 2. # exposure time of sub-exposures (s)
n0 = 900 # number of sub-exposures per exposure for S.O.P.
period = 5.55 * 60. # solar period (s)

def write_to_pickle(fn, stuff):
    print "write_to_pickle(): writing", fn
    filehandler = open(fn, "wb")
    pickle.dump(stuff, filehandler)
    filehandler.close()
    print "...done"

def read_from_pickle(fn):
    print "read_from_pickle(): reading", fn
    filehandler = open(fn, "rb")
    stuff = pickle.load(filehandler)
    filehandler.close()
    print "...done"
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
    # make ivars with proper exptime variation; make noise
    ivars = 1.e10 * (exptimes / 1800.) # 1e10 at exptime = 30 min
    fluxes = np.random.normal(size=ntimes) / np.sqrt(ivars)
    # make signals
    omega = 2. * np.pi / period # rad / s - peak frequency in angular units
    delta_omega = 0.0003 # rad / s - large frequency difference in angular units
    fluxes += 1.
    fluxes += make_one_signal(times, exptimes, 0.00002, 0.00004, omega - 2. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.00003, 0.00004, omega - 1. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.00004, 0.00006, omega + 0. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.00004, 0.00004, omega + 1. * delta_omega)
    fluxes += make_one_signal(times, exptimes, 0.00004, 0.00002, omega + 2. * delta_omega)
    return times, exptimes, fluxes, ivars

def hogg_savefig(prefix):
    for suffix in ["png",]: # "pdf"]:
        fn = prefix + "." + suffix
        print "hogg_savefig(): writing", fn
        pl.savefig(fn)

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
        pl.title(title, fontsize=fontsize)
    pl.subplot(2, 1, 2)
    xlim = (15., 17.)
    # pl.plot(times / 86400., (fluxes - 1.) * 1.e6, "k.")
    for time, exptime, flux in zip(times, exptimes, fluxes):
        t1 = (time - 0.5 * exptime) / 86400.
        t2 = (time + 0.5 * exptime) / 86400.
        if t2 > xlim[0] and t1 < xlim[1]:
            dfppm = (flux - 1.) * 1.e6
            pl.plot([t1, t2], [dfppm, dfppm], "k-", lw=3.0, solid_capstyle="butt")
            pl.axvline(t2, color="k", alpha=0.5, lw=0.5)
    pl.xlim(xlim)
    pl.xlabel("time (d)")
    pl.ylim(-100., 100.)
    pl.ylabel("[(flux - mean) / mean] (ppm)")
    hogg_savefig(prefix)
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
    print "compute_crbs(): starting"
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
    print "compute_aliases(): starting"
    foo = Compute_one_alias(period0, times, exptimes, ivars)
    return np.array(pmap(foo, periods))

class Infer_one_amplitude:
    """
    `Infer_one_amplitude`

    Measure, for one period in the `periods` list, the best-fit
    squared amplitude at the given period.

    It's a class to permit pickling.

    Idiotically slow.
    """
    def __init__(self, times, exptimes, fluxes, ivars):
        self.times = times
        self.exptimes = exptimes
        self.fluxes = (fluxes - 1.) * 1.e6 # (ppm)
        self.ivars = ivars

    def __call__(self, P):
        omega = 2. * np.pi / P
        modelA = make_one_signal(self.times, self.exptimes, 1., 0., omega)
        modelB = make_one_signal(self.times, self.exptimes, 0., 1., omega)
        ones = np.ones_like(modelA)
        A = np.vstack((ones, modelA, modelB))
        ATA = np.dot(A, self.ivars[:, None] * A.T)
        ATy = np.dot(A, self.ivars * self.fluxes)
        xx = np.linalg.solve(ATA, ATy)[1:]
        return np.dot(xx, xx)

def infer_amplitudes(periods, times, exptimes, fluxes, ivars):
    print "infer_amplitudes(): starting"
    foo = Infer_one_amplitude(times, exptimes, fluxes, ivars)
    return np.array(pmap(foo, periods))

def plot_stuff(periods, crbs1, crbs2, name1, name2, ylabel, prefix, big, logy):
    """
    `plot_stuff`
    """
    vline_periods = np.array([ntimes * n0 * dt,
                              2. * n0 * dt,
                              2. * dt])
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
        pl.plot(xs, crbs1, "k-", alpha=0.75, label=name1)
        pl.plot(xs, crbs2, "g-", alpha=0.50, label=name2)
        for xx in vxs:
            pl.axvline(xx, color="k", ls="--", alpha=0.5, lw=0.5)
        pl.legend(loc=4, bbox_to_anchor=(1., 1.01), ncol=2, fontsize=fontsize,
                  borderaxespad=0.)
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
        if logy:
            ylim = (1.e-5 * big, 3.e0 * big)
        else:
            ylim = (-0.05 * big, 1.05 * big)
            pl.axhline(0., color="k", alpha=0.5, lw=0.5)
        # transform craziness to plot triangles
        ax = pl.gca()
        trans = transforms.blended_transform_factory(
            ax.transData, ax.transAxes)
        pl.plot(x0, -0.01, "k^", alpha=0.5, clip_on=False, transform=trans)
        pl.plot(x0, 1.01, "kv", alpha=0.5, clip_on=False, transform=trans)
        pl.xlim(min(xs), max(xs))
        pl.ylim(ylim)
        hogg_savefig("%s%1d" % (prefix, ii))
    return None

def plot_crbs(periods, crbs1, crbs2, name1, name2):
    big = np.max(np.sum(crbs1, axis=1))
    return plot_stuff(periods,
                      np.sum(crbs1, axis=1),
                      np.sum(crbs2, axis=1),
                      name1, name2, "Cramer-Rao bounds", "crb", big, True)

def plot_aliases(periods, aliases1, aliases2, name1, name2):
    big = 1.0
    return plot_stuff(periods,
                      np.max(aliases1, axis=1),
                      np.max(aliases2, axis=1),
                      name1, name2, "worst cosine distance", "alias", big, False)

def plot_amps(periods, amps1, amps2, name1, name2):
    big = np.max(amps1[periods > 3.e2])
    return plot_stuff(periods, amps1, amps2,
                      name1, name2, "best-fit squared amplitude", "amp", big, False)

if __name__ == "__main__":
    picklefn = "./exptime.pkl"
    if not os.path.exists(picklefn):
        #dlnp = 0.05 / ntimes # many of our best people died to bring us MAGIC number "0.05"
        #periods = np.exp(np.arange(np.log(3.0e0) + 0.5 * dlnp, np.log(3.0e4), dlnp)) # (s)
        df = 0.00005 / ntimes # (Hz)
        periods = 1. / np.arange(1./1.e5 + 0.5 * df, 1./1.e2, df) # (s)
        print "periods", len(periods)
        times1, exptimes1, fluxes1, ivars1 = make_fake_data()
        amps1 = infer_amplitudes(periods, times1, exptimes1, fluxes1, ivars1)
        crbs1 = compute_crbs(periods, times1, exptimes1, ivars1)
        aliases1 = compute_aliases(period, periods, times1, exptimes1, ivars1)
        times2, exptimes2, fluxes2, ivars2 = make_fake_data(random=False)
        amps2 = infer_amplitudes(periods, times2, exptimes2, fluxes2, ivars2)
        crbs2 = compute_crbs(periods, times2, exptimes2, ivars2)
        aliases2 = compute_aliases(period, periods, times2, exptimes2, ivars2)
        write_to_pickle(picklefn, (periods,
                                   times1, exptimes1, fluxes1, ivars1, amps1, crbs1, aliases1,
                                   times2, exptimes2, fluxes2, ivars2, amps2, crbs2, aliases2))
    else:
        (periods,
         times1, exptimes1, fluxes1, ivars1, amps1, crbs1, aliases1,
         times2, exptimes2, fluxes2, ivars2, amps2, crbs2, aliases2) = read_from_pickle(picklefn)

    print min(amps1), min(amps2)
    hoggstr = r"Hogg \textit{et al.}~proposal"
    tessstr = r"\textsl{TESS} default"
    plot_amps(periods, amps1, amps2, hoggstr, tessstr)
    plot_crbs(periods, crbs1, crbs2, hoggstr, tessstr)
    plot_aliases(periods, aliases1, aliases2, hoggstr, tessstr)
    plot_exptimes(times1, exptimes1, fluxes1, "hogg", hoggstr)
    plot_exptimes(times2, exptimes2, fluxes2, "TESS", tessstr)
