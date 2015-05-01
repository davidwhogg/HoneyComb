import numpy as np
import matplotlib.pyplot as plt
from fft import load_data

def random_exptime(x, y, ivar, n):
    N = len(x)
    r = np.random.uniform(size=N-1)  # generate N-1 random numbers (this
    # is the largest number of splits you could have)
    j = N/float(n) - 1  # the mean number of points per bin
    k = 100* j / (N - 1)  # the mean number of ppb as a percentage of the
    # total number of points. e.g. 20 bins and 200 points = 10%
    print k
    index = r < np.percentile(r, k)  # select splits in the kth percentile
    print sum(index), j
    print np.diff(x[index])
#     inds = np.arange(len(x))[index]
#     yerr = 1./np.sqrt(ivar)
#     _x, _y, _yerr = [], [], []
#     for i in range(1, len(inds)):
#         _x.append(np.mean(x[inds[i-1]:inds[i]]))
#         _y.append(np.mean(y[inds[i-1]:inds[i]]))
#         yerr = np.array(yerr)
#         _yerr.append(np.sqrt(np.sum(yerr[inds[i-1]:inds[i]]**2))\
#                      / len(yerr[inds[i-1]:inds[i]]))
    numerator = np.cumsum(x[index])
    print np.shape(numerator)
    print numerator, "num"
    print np.ones_like(x)[index]
    denominator = np.cumsum(np.ones_like(x)[index])
    print denominator, "den"
#     assert denominator > 0
    print np.shape(numerator)
    return numerator / denominator

def my_rand(x, y, yerr, n):
    inds = np.random.randint(0, len(x), size=n)
    _x, _y, _yerr = [], [], []
    for i in range(len(inds)-1):
        _x.append(np.mean(x[inds[i]:inds[i+1]]))
        _y.append(np.mean(y[inds[i]:inds[i+1]]))
        _yerr.append(np.sqrt(np.sum(yerr[inds[i]:inds[i+1]]**2))\
                     / len(yerr[inds[i]:inds[i+1]]))
    _x.extend(np.mean(x[inds[-1]:]))
    _y.extend(np.mean(y[inds[-1]:]))
    _yerr.extend(np.sqrt(np.sum(yerr[inds[-1]:]**2))\
                 / len(yerr[inds[-1]:]))
    return _x, _y, _yerr

if __name__ == "__main__":
    name = "6442183"
    print "loading data..."
    D = "/Users/angusr/.kplr/data/lightcurves/%s" % str(name).zfill(9)
    x, y, yerr, ivar = load_data("%s" % name, D, sc=True, nquarters=1)
    l = 100
    x, y, yerr, ivar = x[:l], y[:l], yerr[:l], ivar[:l]
    plt.clf()
    plt.plot(x, np.ones_like(x), "k.")
    x, y, yerr = my_rand(x, y, yerr, 3)
    plt.plot(x, np.zeros_like(x), "r.")
    plt.ylim(-1, 2)
    plt.show()
#     plt.clf()
#     plt.hist(np.diff(xb))
#     print np.diff(xb)
#     plt.show()
