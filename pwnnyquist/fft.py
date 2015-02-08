import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
import nufft
from params import plot_params
reb = plot_params()
import superpgram as spg
from colours import plot_colours
c = plot_colours()
import scipy.signal as sps
import superpgram as spg
from header_time import real_footprint, real_footprint_sc
import george
from george.kernels import ExpSquaredKernel
import glob

# load KASOC light curves
def corrected_data(kid):
    D = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"
    fnames = glob.glob("%s/kplr%s*.dat" % (D, kid.zfill(9)))
    x, y, yerr = [], [], []
    for fname in fnames:
        data = np.genfromtxt(fname, skip_header=9).T
        x.extend(data[0])
        y.extend(data[5])
        yerr.extend(data[6])
    x, y, yerr = np.array(x), np.array(y), np.array(yerr)
    l = np.isfinite(x) * np.isfinite(y) * np.isfinite(yerr)
    x, y, yerr = x[l], y[l], yerr[l]
    med = np.median(y)
    mean = np.mean(y)
    y /= med
    y -= 1
    yerr /= med
#     y -= mean
    x *= 24*3600  # convert to seconds
    ivar = 1./yerr**2
    y = np.array([i.astype("float64") for i in y])
    ivar = np.array([i.astype("float64") for i in ivar])
    return x, y, yerr, ivar

def load_data(kid, sc):
    # load sc data
    D = "/Users/angusr/.kplr/data/lightcurves/%s" % kid.zfill(9)
    if sc == True:
        fnames = glob.glob("%s/kplr%s-*_slc.fits" % (D, kid.zfill(9)))
    else:
        fnames = glob.glob("%s/kplr%s-*_llc.fits" % (D, kid.zfill(9)))
    q, x, y, yerr = [], [], [], []
    for fname in fnames[:3]:
        print fname
        hdulist = pyfits.open(fname)
        tbdata = hdulist[1].data
        x.extend(tbdata["TIME"])
        data = tbdata["PDCSAP_FLUX"]
        q.extend(tbdata["SAP_QUALITY"])
#         med = np.median(data[np.isfinite(data)])
#         y.extend(data/med-1)
        mean = np.mean(data[np.isfinite(data)])
        y.extend(data-mean)
        abs_yerr = tbdata["PDCSAP_FLUX_ERR"]
#         yerr.extend(abs_yerr/med)
        yerr.extend(abs_yerr)
    x, y, yerr, q = np.array(x), np.array(y), np.array(yerr), np.array(q)
    l = np.isfinite(x) * np.isfinite(y) * np.isfinite(yerr) * (q==0)
    x, y, yerr = x[l], y[l], yerr[l]
    x *= 24*3600  # convert to seconds
    ivar = 1./yerr**2
    y = np.array([i.astype("float64") for i in y])
    ivar = np.array([i.astype("float64") for i in ivar])
    return x, y, yerr, ivar

def freqs(kid, fs):
    try:
        truths = np.genfromtxt("%s_freqs.txt" % kid).T  # uHz
        print min(truths*1e-6), max(truths*1e-6)
    except:
        truths = []
    ws = 2*np.pi*fs
    return ws, fs, truths

if __name__ == "__main__":

#     kid = "8006161"
    kid = "5515314"
    # kid = "7103006"
    # kid = "3427720"
    # kid = "3632418"
#     kid = "7341231"

#     # pwnnyquist short cadence data
# #     x, y, yerr, ivar = corrected_data(kid)
#     x, y, yerr, ivar = load_data(kid, sc=True)
# #     fs = np.linspace(0.0002, 0.0005, 1000)  # Hz
# #     fs = np.linspace(0.00265, 0.00285, 10000)  # Hz
    fs = np.linspace(0.00001, 0.00025, 10000)  # Hz
#     ws, fs, truths = freqs(kid, fs)
#     starts, stops, centres = real_footprint_sc(x)
#     amp2s = spg.superpgram(starts, stops, y, ivar, ws)
#
#     # plot light curve
#     plt.clf()
#     plt.subplot(3, 1, 1)
#     plt.plot(x, y, "k.", alpha=.1)
#
#     # plot sc superpgram
#     plt.subplot(3, 1, 2)
#     plt.plot(fs, amp2s, "k")
#     if len(truths):
#         for truth in truths:
#             plt.axvline(truth*1e-6, color=c.blue, linestyle="--")
#     plt.xlim(min(fs), max(fs))
#
    # pwnnyquist long cadence data
    x, y, yerr, ivar = load_data(kid, sc=False)
    ws, fs, truths = freqs(kid, fs)
    starts, stops, centres = real_footprint(x)
    plt.clf()
    plt.subplot(3, 1, 1)
    plt.errorbar(x, y, yerr=yerr, fmt="k.", capsize=0)
#     plt.savefig("test")
    amp2s = spg.superpgram(starts, stops, y, ivar, ws)

    # plot lc superpgram
    plt.subplot(3, 1, 3)
    plt.plot(fs, amp2s, "k")
    if len(truths):
        for truth in truths:
            plt.axvline(truth*1e-6, color=c.blue, linestyle="--")
    plt.xlim(min(fs), max(fs))
    plt.savefig("fft")
