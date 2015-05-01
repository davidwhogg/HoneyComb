"""
This file is part of the HoneyComb project.
Copyright 2015 David W. Hogg (NYU).

## bugs:
* `nthreads > 1` doesn't work
"""

import os.path as path
import numpy as np
import superpgram as spg
import pylab as plt
import cPickle as pickle
import time
import glob
import pyfits
from header_time import real_footprint

nthreads = 1
if nthreads > 1:
    from multiprocessing import Pool
    pmap = Pool(nthreads).map
else:
    pmap = map

class Function:
    """
    # `Function`

    A function specification that is a constant level plus a sum of
    `J` sinusoids.
    """

    def __init__(self, a0, abks):
        self.J = 0
        self.a0 = a0
        self.abks = None
        if abks is not None:
            self.abks = np.atleast_2d(np.array(abks))
            J, foo = self.abks.shape
            assert foo == 3
            self.J = J

    def __str__(self):
        return "Function %d %f " % (self.J, self.a0) + str(self.abks)

    def evaluate(self, xs):
        ys = np.zeros_like(xs) + self.a0
        if self.J == 0:
            return ys
        for a, b, k in self.abks:
            ys += a * np.cos(k * xs) + b * np.sin(k * xs)
        return ys

    def definite_integral(self, d, u):
        I = (u - d) * self.a0
        if self.J == 0:
            return I
        for a, b, k in self.abks:
            I += (a / k) * np.sin(k * u) - (b / k) * np.cos(k * u)
            I -= (a / k) * np.sin(k * d) - (b / k) * np.cos(k * d)
        return I

class LightCurveFootprint:
    """
    # `LightCurveFootprint`

    A list of `Projection` objects that describe a Kepler-like
    lightcurve's set of exposure times.  Note the importance of the
    conversion to barycentric time.
    """

    def __init__(self):
        delta_t = 30. / 60. / 24. # days
        # make regular grid in days
        t_centers = np.arange(0.5 * delta_t, 4.1 * 365.25, delta_t) # 4.1 years
        t_exp = 0.95 * delta_t
        # distort time grid onto Barycentric Time
        self.starts = self._time_distort(t_centers - 0.5 * t_exp)
        self.stops = self._time_distort(t_centers + 0.5 * t_exp)
        self.centers = 0.5 * (self.stops + self.starts)

    def _time_distort(self, xs):
        return xs + (3.2 / 60. / 24.) * np.cos(2. * np.pi * xs / 371.) # 3.2 light-minutes in days

    def project(self, f):
        ys = np.zeros(len(self.starts))
        for n, (d, u) in enumerate(zip(self.starts, self.stops)):
            ys[n] = f.definite_integral(d, u) / (u - d)
        return ys

def make_fake_data():
    """
    # `make_fake_data()`

    make Kepler-like data
    """
    np.random.seed(42)
    k0 = (2 * np.pi / 5.132342) * 60 * 24 # 5-min period -> frequency in days
    dk = 0.002 * k0
    ks = k0 + dk * np.arange(-1., 1.1, 1.)
    abks = np.zeros((len(ks), 3))
    abks[:, 0] = 1.e-4 * np.random.normal(size=len(ks))
    abks[:, 1] = 1.e-4 * np.random.normal(size=len(ks))
    abks[:, 2] = ks
    fun = Function(1.0, abks)
    lcf = LightCurveFootprint()
    data = lcf.project(fun) # project sinusoids
    ivar = 1.e10 + np.zeros(len(data)) # incredibly low noise
    data += (1.0 / np.sqrt(ivar)) * np.random.normal(size=len(data)) # add noise
    return lcf, data, ivar, fun

def fit_one_sinusoid(k, lcf, data, ivar):
    """
    # `fit_one_sinusoid`

    fit a constant plus a sine and cosine at one k value.
    """
    f0 = Function(1.0, None)
    fa = Function(0.0, [[1., 0., k], ])
    fb = Function(0.0, [[0., 1., k], ])
    A = np.array([lcf.project(f0), lcf.project(fa), lcf.project(fb)])
    ATA = np.dot(A, ivar[:,None] * A.T)
    ATb = np.dot(A, ivar * data)
    pars = np.linalg.solve(ATA, ATb)
    resid = data - np.dot(A.T, pars)
    chi2 = np.dot(resid, ivar * resid)
    a0, a, b = pars
    return a0, a, b, chi2

if __name__ == "__main__":

    # make fake data
    strt = time.time()
    print "load real data:", time.time() - strt

    # load test star long cadence fits files
    D = "/Users/angusr/.kplr/data/lightcurves/008006161"
    kid = "8006161"  # a kepler target with lc and sc data chosen at `random'
    fnames = []
    qs = range(17)
    x, data, yerr = [], [], []
    for q in qs:
        try:
            fnames.append(glob.glob("%s/kplr%s-*_llc.fits"
                          % (D, q, kid.zfill(9)))[0])
            for fname in fnames:
                hdulist = pyfits.open(fname)
                tbdata = hdulist[1].data
                x.extend(tbdata["TIME"])
                y = tbdata["PDCSAP_FLUX"]
                med = np.median(y)
                data.extend(y/med)
                abs_yerr = tbdata["PDCSAP_FLUX_ERR"]
                yerr.extend(abs_yerr/med)
        except:
            print "no file found: Q%s, kplr%s" % (q, kid.zfill(9))
            pass
    x, data, yerr = np.array(x), np.array(data), np.array(yerr)
    l = np.isfinite(x) * np.isfinite(data) * np.isfinite(yerr)
    x, data, yerr = x[l], data[l], yerr[l]
    ivar = 1./yerr**2

    # convert BJDs to UTCs
    x = np.array(x) + 2454833
    starts, stops, centers = real_footprint(x)

    # perform inferences in Fourier space
    truths = np.genfromtxt("3427720_freqs.txt").T
    truths *= 1.842*2*np.pi # convert to rads/day
    dk = 2. / (4.1 * 365.25) # frequency resolution for testing
    for ri, truth in enumerate(truths):
        print ri, truth
        testks = np.arange(-400.5, 401., 1.) * dk + truth
        strt = time.time()

        # convert to float64
        data = np.array([d.astype('float64') for d in data])
        ivar = np.array([i.astype('float64') for i in ivar])

        amp2s = spg.superpgram(starts, stops, data, ivar, testks)
        print "computed super-resolution periodogram:", time.time() - strt

        # plot data
        strt = time.time()
        plt.clf()
        plt.plot(centers, data, "k.", ms=0.75)
        plt.xlabel("time [day]")
        plt.ylabel("intensity")
        plt.savefig("foo%s.png" % ri)

        # plot fourier tests
        plt.clf()
        plt.step(testks, np.log10(amp2s), color="k", where="mid")
        for k in truths:
            plt.axvline(k, alpha=0.5)
        plt.xlabel("wave number [uHz]")
        plt.ylabel("log10 squared amplitude of best-fit sinusoid")
        plt.xlim(np.min(testks), np.max(testks))
        plt.savefig("ruthbar%s.png" % ri)
        print "made plots:", time.time() - strt
