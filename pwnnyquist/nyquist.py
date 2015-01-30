"""
This file is part of the HoneyComb project.
Copyright 2015 David W. Hogg (NYU).
"""

import numpy as np
import pylab as plt

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
            print J, foo
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
    lightcurve's set of exposure times.
    """

    def __init__(self):
        delta_t = 30. / 60. / 24. # days
        # make regular grid in days
        t_centers = np.arange(0.5 * delta_t, 100., delta_t)
        t_exp = 0.95 * delta_t
        ds = t_centers - 0.5 * t_exp
        us = t_centers + 0.5 * t_exp
        ds = self._time_distort(ds)
        us = self._time_distort(us)
        self.x_centers = 0.5 * (ds + us)
        self.Ps = []
        for d, u in zip(us, ds):
            self.Ps.append(Projection([[d, u, 1./(u - d)], ]))

    def _time_distort(self, xs):
        return xs + (8. / 60. / 24.) * np.cos(2. * np.pi * xs / 371.) # 8 light-minutes in days

    def project(self, f):
        ys = np.zeros(len(self.Ps))
        for n, P in enumerate(self.Ps):
            ys[n] = P.project(f)
        return ys

if __name__ == "__main__":
    np.random.seed(42)
    k0 = (2 * np.pi / 5.132342) * 60 * 24 # 5-min period -> frequency in days
    dk = 0.01 * k0
    ks = k0 + dk * np.arange(-5., 5.1, 1.)
    abks = np.zeros((len(ks), 3))
    abks[:, 0] = 0.1 * np.random.normal(size=len(ks))
    abks[:, 1] = 0.1 * np.random.normal(size=len(ks))
    abks[:, 2] = ks
    f = Function(1.0, abks)
    lcf = LightCurveFootprint()
    plt.clf()
    # xplot = np.arange(0., 1000., 0.0001)
    # plt.plot(xplot, f.evaluate(xplot), "k-", alpha=0.25)
    plt.plot(lcf.x_centers, lcf.project(f), "k.", ms=0.75)
    plt.savefig("foo.png")
    print "Hello World!"
