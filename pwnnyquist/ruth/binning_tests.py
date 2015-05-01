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
def compare(fs, nbins, load=False, test=False):
    name = "6442183"  # lower nu max poster child
    print "loading data..."
    D = "/Users/angusr/.kplr/data/lightcurves/%s" % str(name).zfill(9)
    x, y, yerr, ivar = load_data("%s" % name, D, sc=True)
    if test:
        n = 100
        x, y, yerr, ivar = x[:n], y[:n], yerr[:n], ivar[:n]
    print len(fs)
    ws = np.pi*2*fs
    dt = 6.81119940564e-4 * 24 * 3600 * nbins # short cadence itrvl (seconds)
#     dt = 6.81119940564e-4 * 24 * 3600  # short cadence itrvl (seconds)
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

    if load:
        print "loading fft"
        with h5py.File("fft.h5", "r") as f:
            fs = f["fft"][0, :]
            fft = f["fft"][1, :]
    else:
        print "fft..."
        fft = nufft.nufft3(stops, y, ws)
        fft = np.imag(fft)**2 + np.real(fft)**2
        print "saving fft"
        f = h5py.File("fft%s.h5" % str(nbins).zfill(2), "w")
        data = f.create_dataset("fft", (2, len(fs)))
        data[0, :] = fs
        data[1, :] = fft
        f.close()
    plt.plot(fs*1e6-1160, fft, "k")
    plt.ylabel("$\mathrm{FFT}$")
    plt.xticks(visible=False)
    plt.subplot(3, 1, 3)

    if load:
        print "loading ls"
        with h5py.File("ls.h5", "r") as f:
            fs = f["ls"][0, :]
            pgram = f["ls"][1, :]
    else:
        print "lombScargle..."
        model = LombScargle().fit(x, y, yerr)
        ps = 1/fs
        pgram = model.periodogram(ps)
        print "saving ls"
        f = h5py.File("ls%s.h5" % str(nbins).zfill(2), "w")
        data = f.create_dataset("ls", (2, len(fs)))
        data[0, :] = fs
        data[1, :] = pgram
        f.close()
    plt.plot(fs*1e6-1160, pgram/max(pgram), "k")
    plt.ylabel("$\mathrm{LombScargle}$")
    plt.xlabel("$\\nu-1160~\mathrm{(}\mu\mathrm{Hz)}$")
    plt.subplots_adjust(hspace=0)
    plt.savefig("fft_ls_spg_%s" % str(nbins).zfill(2))

if __name__ == "__main__":
    fs = np.arange(1159, 1161, 1e-4) * 1e-6
    nbins = int(sys.argv[1])
    compare(fs, nbins)
