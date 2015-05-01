import numpy as np
import matplotlib.pyplot as plt
from fft_old import corrected_data
from newbintests import gap_split, remainder, binning

def gap_split(x, y, yerr):
    ''' Split data at the gaps, return list of arrays.
    the index of diff corresponds to the data point just before the gap
    '''
    diff = np.diff(x)
    l = diff > 60  # go above 180 seconds for fewer gaps
    inds = np.arange(len(x))[l]  # find the indices at which the gaps occur
    xs, ys, yerrs = [], [], []
    xs.append(x[:inds[0]])
    ys.append(y[:inds[0]])
    yerrs.append(yerr[:inds[0]])
    for i in range(1, len(inds)):
        xs.append(x[inds[i-1]:inds[i]])
        ys.append(y[inds[i-1]:inds[i]])
        yerrs.append(yerr[inds[i-1]:inds[i]])
    xs.append(x[inds[-1]:])
    ys.append(y[inds[-1]:])
    yerrs.append(yerr[inds[-1]:])
    return xs, ys, yerrs

def remainder(x, y, yerr, n):
    """
    Truncates the timeseries, modulo n
    """
    m = len(x) % n
    if m > 0: return x[:-m], y[:-m], yerr[:-m]
    else: return x, y, yerr

def rebin(x, y, yerr, n):
    """
    Bins the data into groups of n data points
    """
    x1, y1, yerr1 = remainder(x, y, yerr, n)
    nb = len(x1)/n
    xb, yb, yerrb = np.zeros(nb), np.zeros(nb), np.zeros(nb)
    for i in range(n):
        xb += x1[i::n]
        yb += y1[i::n]
        yerrb += (yerr1[i::n])**2
    return xb/n, yb/n, np.sqrt(yerrb)/n

if __name__ == "__main__":
    # load corrected data
    # D = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"
    # x, y, yerr, ivar = corrected_data("8006161", D)  # median to mean
    # l = 100
    # x, y, yerr, ivar = x[:l], y[:l], yerr[:l], ivar[:l]

    # make up data
    x = np.arange(0, 3000, 58.84704018)
    x = np.concatenate((x[:10], x[11:]))
    x = np.concatenate((x[:20], x[21:]))
    y = np.random.randn(len(x))
    yerr = np.ones_like(y)*.1
    ivar = 1./yerr**2

    xs, ys, yerrs = gap_split(x, y, yerr)
    print np.shape(xs[0])
    print np.shape(ys[0])
    print np.shape(yerrs[0])

    xb, yb, yerrb = rebin(x, y, yerr, ivar, 3)

    plt.clf()
    plt.errorbar(x, y, yerr=yerr, fmt="k.")
    plt.errorbar(xb, yb, yerr=yerrb, fmt="r.")
    plt.errorbar(xs[0], ys[0], yerr=yerrs[0], fmt="b.")
    plt.errorbar(xs[1], ys[1], yerr=yerrs[1], fmt="g.")
    plt.errorbar(xs[2], ys[2], yerr=yerrs[2], fmt="m.")
    plt.show()
