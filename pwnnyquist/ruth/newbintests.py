import numpy as np
import matplotlib.pyplot as plt
import glob
import nufft
from gatspy.periodic import LombScargle
import superpgram as spg
from header_time import real_footprint_sc, real_footprint
from fft import corrected_data, load_data, triangular_smoothing
import h5py
import time
from push_the_cadence import arbitrary_footprint
import sys

# split data at the gaps, return list of arrays.
# the index of diff corresponds to the data point just before the gap
def gap_split(x, y, yerr):
    diff = np.diff(x)
    l = diff > 180  # go above 180 for fewer gaps
    inds = np.arange(len(x))[l]
    x_split = [(x[inds[i-1]:inds[i]]) for i in range(1, len(inds))]
    y_split = [(y[inds[i-1]:inds[i]]) for i in range(1, len(inds))]
    yerr_split = [(yerr[inds[i-1]:inds[i]]) for i in range(1, len(inds))]
    print type(x_split), len(x_split), np.shape(x_split)
    assert 0
    return x_split, y_split, yerr_split

# take in an array and truncate remainder
def remainder(x, y, yerr, n):
    m = len(x) % n
    if m > 0:
        return x[:-m], y[:-m], yerr[:-m]
    else:
        return x, y, yerr

def binning(x, y, yerr, ivar, n):
    x = x[n:n]
    y = binit(y, n)
    yerr = binit(yerr, n, quad=True)
    ivar = 1./yerr**2
    return x, y, yerr, ivar

def binit(x, n, quad=False):  # n = number of points to bin
    xarray = np.zeros(len(x)/n)
    for i, j  in enumerate(np.arange(len(x)/n) * n):
        if quad:
            xarray[i] = np.sqrt(np.sum(x[j:j+n]**2))/float(n)
        else:
            xarray[i] = np.mean(x[j:j+n])
    return xarray

# calc spg, fft and ls. option to generate and save new results or load
# precomputed
def compare(fs, nbins, load=False):
    name = "6442183"  # lower nu max poster child
    print "loading data..."
    D = "/Users/angusr/.kplr/data/lightcurves/%s" % str(name).zfill(9)
    x, y, yerr, ivar = load_data("%s" % name, D, sc=True)
    x, y, yerr = gap_split(x, y, yerr)
    print len(fs)
    ws = np.pi*2*fs
    dt = 6.81119940564e-4 * 24 * 3600 * nbins # short cadence itrvl (seconds)
    starts, stops, centres = arbitrary_footprint(x, dt)

    if load:
        print "loading spg"
        with h5py.File("spg.h5", "r") as f:
            fs = f["spg"][0, :]
            amp2s = f["spg"][1, :]
    else:
        print "computing superpgram..."
        amp2s = spg.superpgram(starts, stops, y, ivar, ws)
        print "saving spg"
        f = h5py.File("spg%s.h5" % str(nbins).zfill(2), "w")
        data = f.create_dataset("spg", (2, len(fs)))
        data[0, :] = fs
        data[1, :] = amp2s
        f.close()
    plt.clf()
    plt.subplot(3, 1, 1)
    plt.plot(fs*1e6-1160, amp2s, "k")
    plt.ylabel("$\mathrm{Spgram}$")
    plt.xticks(visible=False)
    plt.subplot(3, 1, 2)

if __name__ == "__main__":
    fs = np.arange(1159, 1161, 1e-4) * 1e-6
    nbins = int(sys.argv[1])
    compare(fs, nbins)
