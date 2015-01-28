"""
This file is part of the HoneyComb project.
Copyright 2015 David W. Hogg (NYU).
"""

import numpy as np
import pylab as plt

class Function:

    def __init__(self, J, a0, abklist):
        self.J = int(J)
        self.a0 = a0
        self.abklist = None
        if self.J > 0:
            self.abklist = np.atleast_2d(np.array(abklist))
        if self.J > 0:
            Jtest, foo = self.abklist.shape
            print J, foo
            assert foo == 3
            assert self.J == Jtest

    def evaluate(self, x):
        y = np.zeros_like(x) + self.a0
        if self.J > 0:
            for a, b, k in self.abklist:
                y += a * np.cos(k * x) + b * np.sin(k * x)
        return y

    def definite_integral(self, d, u):
        I = (u - d) * self.a0
        if self.J > 0:
            for a, b, k in self.abklist:
                I += (a / k) * np.sin(k * u) - (b / k) * np.cos(k * u)
                I -= (a / k) * np.sin(k * d) - (b / k) * np.cos(k * d)
        return I

if __name__ == "__main__":
    f = Function(2, 2.0, [[0.05, 0.0, 1.0], [0.02, 0.02, 13.7]])
    f0 = Function(0, 1.0, [])
    plt.clf()
    xplot = np.arange(0., 10., 0.001)
    plt.plot(xplot, f.evaluate(xplot), "k-")
    xtest = np.random.uniform(0., 10., size=(5))
    dxtest = 0.5
    ytest = np.zeros_like(xtest)
    for i, x in enumerate(xtest):
        d, u = x - 0.5 * dxtest, x + 0.5 * dxtest
        ytest[i] = f.definite_integral(d, u) / f0.definite_integral(d, u)
        plt.plot(xtest[i], ytest[i], "rx")
        plt.plot([d, u], [ytest[i], ytest[i]], "r-")
    plt.savefig("foo.png")
    print "Hello World!"
