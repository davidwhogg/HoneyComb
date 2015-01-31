"""
This file is part of the HoneyComb project.
Copyright 2015 David W. Hogg (NYU).

## bugs:
* `nthreads > 1` doesn't work
"""

import os.path as path
import numpy as np
import pylab as plt
import cPickle as pickle

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

class Projection:
    """
    # `Projection`

    A local linear projection operator that is a linear combination of
    M top-hat definite integrals.
    """

    def __init__(self, duws):
        self.M = 0
        self.duws = None
        if duws is not None:
            self.duws = np.atleast_2d(np.array(duws))
            M, foo = self.duws.shape
            assert foo == 3
            self.M = M

    def __str__(self):
        return "Projection %d " % (self.M, ) + str(self.duws)

    def project(self, f):
        if self.M == 0:
            return 0.
        I = 0.
        for d, u, w in self.duws:
            I += w * f.definite_integral(d, u)
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
        ds = t_centers - 0.5 * t_exp
        us = t_centers + 0.5 * t_exp
        # distort time grid onto Barycentric Time
        ds = self._time_distort(ds)
        us = self._time_distort(us)
        self.x_centers = 0.5 * (ds + us)
        self.Ps = []
        for d, u in zip(us, ds):
            self.Ps.append(Projection([[d, u, 1./(u - d)], ]))

    def _time_distort(self, xs):
        return xs + (7. / 60. / 24.) * np.cos(2. * np.pi * xs / 371.) # 7 light-minutes in days

    def project(self, f):
        ys = np.zeros(len(self.Ps))
        for n, P in enumerate(self.Ps):
            ys[n] = P.project(f)
        return ys

def make_fake_data():
    """
    # `make_fake_data()`

    make Kepler-like data
    """
    np.random.seed(42)
    k0 = (2 * np.pi / 5.132342) * 60 * 24 # 5-min period -> frequency in days
    dk = 0.02 * k0
    ks = k0 + dk * np.arange(-5., 5.1, 1.)
    abks = np.zeros((len(ks), 3))
    abks[:, 0] = 0.1 * np.random.normal(size=len(ks))
    abks[:, 1] = 0.1 * np.random.normal(size=len(ks))
    abks[:, 2] = ks
    fun = Function(1.0, abks)
    lcf = LightCurveFootprint()
    data = lcf.project(fun) # project sinusoids
    ivar = 1.e20 + np.zeros(len(data)) # incredibly low noise
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

def fit_sinusoids(ks, lcf, data, ivar):
    def _fitty(k):
        print k
        a0, a, b, chi2 = fit_one_sinusoid(k, lcf, data, ivar)
        return a * a + b * b
    return np.array(pmap(_fitty, ks))

if __name__ == "__main__":
    picklefn = "nyquist.pkl"
    try:
        input = open(picklefn, "rb")
        lcf, data, ivar, truth, testks, amp2s = pickle.load(input)
        input.close()
    except:

        # make fake data
        lcf, data, ivar, truth = make_fake_data()

        # perform inferences in Fourier space
        dk = 1. / (4.1 * 365.25) # frequency resolution for testing
        testks = np.arange(-200.5, 101., 1.) * dk + (truth.abks[0])[2]
        amp2s = fit_sinusoids(testks, lcf, data, ivar)

        # save output
        output = open(picklefn, "wb")
        pickle.dump((lcf, data, ivar, truth, testks, amp2s), output)
        output.close()

    # plot data
    plt.clf()
    plt.plot(lcf.x_centers, data, "k.", ms=0.75)
    plt.xlabel("time [day]")
    plt.ylabel("intensity")
    plt.savefig("foo.png")

    # plot fourier tests
    plt.clf()
    plt.step(testks, np.log10(amp2s), color="k")
    for a, b, k in truth.abks:
        plt.axvline(k, alpha=0.5)
    plt.xlabel("wave number [rad per day]")
    plt.ylabel("log10 squared amplitude of best-fit sinusoid")
    plt.xlim(np.min(testks), np.max(testks))
    plt.savefig("bar.png")
