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
from header_time import real_footprint
import george
from george.kernels import ExpSquaredKernel

def load_data(kid, sc):

    # load sc data
    D = "/Users/angusr/.kplr/data/lightcurves/%s" % kid.zfill(9)
    if sc == True:
        fnames = glob.glob("%s/kplr%s-*_slc.fits" % (D, kid.zfill(9)))
    else:
        fnames = glob.glob("%s/kplr%s-*_llc.fits" % (D, kid.zfill(9)))

    q, x, y, yerr = [], [], [], []
    for fname in fnames:
        hdulist = pyfits.open(fname)
        tbdata = hdulist[1].data
        x.extend(tbdata["TIME"])
        data = tbdata["PDCSAP_FLUX"]
        q.extend(tbdata["SAP_QUALITY"])
        med = np.median(data)
        y.extend(data/med-1)
#         mean = np.mean(data)
#         y.extend(data-mean)
        abs_yerr = tbdata["PDCSAP_FLUX_ERR"]
        yerr.extend(abs_yerr/med)
#         yerr.extend(abs_yerr)
    x, y, yerr, q = np.array(x), np.array(y), np.array(yerr), np.array(q)
    l = np.isfinite(x) * np.isfinite(y) * np.isfinite(yerr) * (q==0)
    x, y, yerr = x[l], y[l], yerr[l]
    x, y, yerr = sigma_clipping(x, y, yerr, sc)

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

#     kid = "8006161"  # a kepler target with lc and sc data chosen at `random'
    # kid = "7103006"
    # kid = "3427720"
#     kid = "3632418"
    kid = "7341231"
    x, y, yerr, ivar = load_data(kid, sc=True)
    fs = np.linspace(0.0002, 0.0005, 1000)  # Hz
    ws, fs, truths = freqs(kid, fs)

    # fft
    pgram = nufft.nufft3(x, y, ws)
    p = np.sqrt(np.real(pgram)**2 + np.imag(pgram)**2)

    # pwnnyquist
    starts, stops, centres = real_footprint(x)
    testks = ws
    amp2s = spg.superpgram(starts, stops, y, ivar, testks)

    # plot
    plt.clf()
    plt.subplot(3, 1, 1)
    plt.plot(x, y, "k.", alpha=.1)
    plt.subplot(3, 1, 2)
    plt.plot(fs, p, "k")
    if len(truths):
        for truth in truths:
            plt.axvline(truth*1e-6, color=c.blue, linestyle="--")
    plt.xlim(min(fs), max(fs))
    plt.subplot(3, 1, 3)

    x, y, yerr, ivar = load_data(kid, sc=False)
    ws, fs, truths = freqs(kid, fs)

    # fft
#     pgram = nufft.nufft3(x, y, ws)
#     p = np.sqrt(np.real(pgram)**2 + np.imag(pgram)**2)

    # pwnnyquist
    starts, stops, centres = real_footprint(x)
    testks = ws
    amp2s = spg.superpgram(starts, stops, y, ivar, testks)

    plt.plot(fs, amp2s, "k")
    if len(truths):
        for truth in truths:
            plt.axvline(truth*1e-6, color=c.blue, linestyle="--")
    plt.xlim(min(fs), max(fs))


#     if len(truths):
#         for truth in truths:
#             plt.axvline(truth*1e-6, color=c.blue, linestyle="--")
#     plt.plot(fs, amp2s, "k")
#     plt.xlim(min(fs), max(fs))
    plt.savefig("fft")
