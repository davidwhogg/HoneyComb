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

    def __init__(self, J, a0, abklist):
        self.J = int(J)
        assert self.J == J
        self.a0 = a0
        self.abklist = None
        if self.J > 0:
            self.abklist = np.atleast_2d(np.array(abklist))
        if self.J > 0:
            Jtest, foo = self.abklist.shape
            print J, foo
            assert foo == 3
            assert self.J == Jtest

    def __str__(self):
        return "Function %d %f " % (self.J, self.a0) + str(self.abklist)

    def evaluate(self, xs):
        ys = np.zeros_like(xs) + self.a0
        if self.J == 0:
            return ys
        for a, b, k in self.abklist:
            ys += a * np.cos(k * xs) + b * np.sin(k * xs)
        return ys

    def definite_integral(self, d, u):
        I = (u - d) * self.a0
        if self.J == 0:
            return I
        for a, b, k in self.abklist:
            I += (a / k) * np.sin(k * u) - (b / k) * np.cos(k * u)
            I -= (a / k) * np.sin(k * d) - (b / k) * np.cos(k * d)
        return I

class Projection:
    """
    # `Projection`

    A local linear projection operator that is a linear combination of
    M top-hat definite integrals.
    """

    def __init__(self, M, duwlist):
        self.M = int(M)
        assert self.M == M
        self.duwlist = None
        if self.M > 0:
            self.duwlist = np.atleast_2d(np.array(duwlist))
        if self.M > 0:
            Mtest, foo = self.duwlist.shape
            assert foo == 3
            assert self.M == Mtest

    def __str__(self):
        return "Projection %d " % (self.M, ) + str(self.duwlist)

    def project(self, f):
        if self.M == 0:
            return 0.
        I = 0.
        for d, u, w in self.duwlist:
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
        t_centers = np.arange(0.5 * delta_t, 400., delta_t)
        t_exp = 0.95 * delta_t
        ds = t_centers - 0.5 * t_exp
        us = t_centers + 0.5 * t_exp
        ds = self._time_distort(ds)
        us = self._time_distort(us)
        self.x_centers = 0.5 * (ds + us)
        self.Ps = []
        for d, u in zip(us, ds):
            self.Ps.append(Projection(1, [[d, u, 1./(u - d)], ]))

    def _time_distort(self, xs):
        return xs + (8 / 60. / 24.) * np.cos(2. * np.pi * xs / 371.) # days

    def project(self, f):
        ys = np.zeros(len(self.Ps))
        for n, P in enumerate(self.Ps):
            ys[n] = P.project(f)
        return ys

if __name__ == "__main__":
    f = Function(2, 2.0, [[0.01, 0.0, 1.0], [0.02, 0.02, 50.2343344]])
    lcf = LightCurveFootprint()
    plt.clf()
    xplot = np.arange(0., 400., 0.0001)
    plt.plot(xplot, f.evaluate(xplot), "k-", alpha=0.25)
    plt.plot(lcf.x_centers, lcf.project(f), "ko", mfc="none")
    plt.savefig("foo.png")
    print "Hello World!"
