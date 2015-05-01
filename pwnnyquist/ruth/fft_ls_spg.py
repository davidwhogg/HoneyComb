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

# calc spg, fft and ls. option to generate and save new results or load
# precomputed
def compare(fs, save=False, load=True, test=False):

    name = "6442183"  # lower nu max poster child

    print "loading data..."
    D = "/Users/angusr/.kplr/data/lightcurves/%s" % name.zfill(9)
    x, y, yerr, ivar = load_data("%s" % name, D, sc=True)
    if test:
        n = 10000
        x, y, yerr, ivar = x[:n], y[:n], yerr[:n], ivar[:n]

    print len(fs)
    ws = np.pi*2*fs

#     starts, stops, centres = real_footprint_sc(x)
    if load:
        print "loading spg"
        with h5py.File("spg.h5", "r") as f:
            fs = f["spg"][0, :]
            amp2s = f["spg"][1, :]
    else:
        print "computing superpgram..."
        amp2s = spg.superpgram(starts, stops, y, ivar, ws)
    if save:
        print "saving spg"
        f = h5py.File("spg.h5", "w")
        data = f.create_dataset("spg", (2, len(fs)))
        data[0, :] = fs
        data[1, :] = amp2s
        f.close()
    plt.clf()
    plt.subplot(3, 1, 1)
    # plt.plot(fs*1e6, amp2s, "k")
    plt.plot(fs*1e6-1160, amp2s, "k")
    plt.ylabel("$\mathrm{Spgram}$")
    # plt.xlim(1159, 1161)
    plt.xticks(visible=False)
#     plt.yticks(visible=False)

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
    if save:
        print "saving fft"
        f = h5py.File("fft.h5", "w")
        data = f.create_dataset("fft", (2, len(fs)))
        data[0, :] = fs
        data[1, :] = fft
        f.close()
    # plt.plot(fs*1e6, fft, "k")
    plt.plot(fs*1e6-1160, fft, "k")
    # plt.xlim(min(fs*1e6), max(fs*1e6))
    # plt.xlim(1159, 1161)
    plt.ylabel("$\mathrm{FFT}$")
    plt.xticks(visible=False)
#     plt.yticks(visible=False)

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
    if save:
        print "saving ls"
        f = h5py.File("ls.h5", "w")
        data = f.create_dataset("ls", (2, len(fs)))
        data[0, :] = fs
        data[1, :] = pgram
        f.close()
    # plt.plot(fs*1e6, pgram, "k")
    plt.plot(fs*1e6-1160, pgram/max(pgram), "k")

#     plt.xlim(1159, 1161)
    plt.ylabel("$\mathrm{LombScargle}$")
    plt.xlabel("$\\nu-1160~\mathrm{(}\mu\mathrm{Hz)}$")
    # plt.xlabel("$\\nu~\mathrm{(}\mu\mathrm{Hz)}$")
#     plt.yticks(visible=False)
    plt.subplots_adjust(hspace=0)

    # plt.savefig("%sfft_ls_spg" % name)
    plt.savefig("test")

if __name__ == "__main__":

#         fs = np.arange(700, 1500, 1e-2) * 1e-6
    fs = np.arange(1159, 1161, 1e-4) * 1e-6
    compare(fs, load=True, save=False)
