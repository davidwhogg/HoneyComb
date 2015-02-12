import numpy as np
import matplotlib.pyplot as plt
import pyfits
import superpgram as spg
from header_time import real_footprint, real_footprint_sc
import glob

plotpar = {'axes.labelsize': 12,
           'xtick.labelsize': 12,
           'ytick.labelsize': 8,
           'text.usetex': True}
plt.rcParams.update(plotpar)

# load KASOC light curves
def corrected_data(kid, D):
    fnames = glob.glob("%s/kplr%s*.dat" % (D, kid.zfill(9)))
    x, y, yerr = [], [], []
    for fname in fnames[:3]:
        data = np.genfromtxt(fname, skip_header=9).T
        x.extend(data[0])
        y.extend(data[5])
        yerr.extend(data[6])
    x, y, yerr = np.array(x), np.array(y), np.array(yerr)
    l = np.isfinite(x) * np.isfinite(y) * np.isfinite(yerr)
    x, y, yerr = x[l], y[l], yerr[l]
#     med = np.median(y)
#     mean = np.mean(y)
#     y /= med
#     y -= 1
#     yerr /= med
#     y -= mean
    x *= 24*3600  # convert to seconds
    ivar = 1./yerr**2
    y = np.array([i.astype("float64") for i in y])
    ivar = np.array([i.astype("float64") for i in ivar])
    return x, y, yerr, ivar

def load_data(kid, D, sc):
    # load sc data
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
        med = np.median(data[np.isfinite(data)])
        y.extend(data/med-1)
#         mean = np.mean(data[np.isfinite(data)])
#         y.extend(data-mean)
        abs_yerr = tbdata["PDCSAP_FLUX_ERR"]
        yerr.extend(abs_yerr/med)
#         yerr.extend(abs_yerr)
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

def triangular_smoothing(x, y, yerr):
    win = 1.*24*3600  # one day window
    dt = 6.81119940564e-4  # interval between observations (days) (sc)
    n = int(1./dt)/2
    weights = (np.array(range(n)+[n]+range(n)[::-1])+1)/float(n)
    weights -= max(weights) % 1.
    nw = len(weights)
    newy = np.zeros_like(y)
#     plt.clf()
#     plt.plot(x, y/np.mean(y), "k.")
    for i in range(len(x)):
        w = np.zeros_like(x)
        l = (x[i]-(.5*win) < x) * (x < x[i]+(.5*win))
        if i > nw/2 and i < len(x)-nw/2:
            w[i-nw/2:i+nw/2+1] = weights
            newy[i] = y[i]/np.average(y, weights=w)
    l = newy!=0
#     plt.plot(x[l], newy[l], "r.")
#     plt.savefig("test")
    return x[l], newy[l], yerr[l]

if __name__ == "__main__":

    kid = "8006161"
#     kid = "5515314"
#     kid = "7103006"
#     kid = "3427720"
#     kid = "3632418"
#     kid = "7341231"
    kepler_lc_dir = "/Users/angusr/.kplr/data/lightcurves/%s" % kid.zfill(9)
    corrected_lc_dir = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"

    # pwnnyquist short cadence data
    x, y, yerr, ivar = corrected_data(kid, corrected_lc_dir)

#     x, y, yerr, ivar = load_data(kid, kepler_lc_dir, sc=True)
#     x, y, yerr = triangular_smoothing(x, y, yerr)
    mean = np.mean(y)
    y -= mean
    plt.clf()
    plt.errorbar(x, y, yerr=yerr, fmt="k.", alpha=.2)
    plt.savefig("test")

#     fs = np.linspace(0.00001, 0.00025, 10000)  # Hz 5515314
    fs = np.arange(0.0035895, 0.00359065, 4e-9)  # Hz 8006161
#     fs = np.arange(0.003, 0.004, 4e-7)  # Hz 8006161
#     fs = np.arange(0.0025, 0.003, 4e-8)  # Hz 5515314
    print len(fs)
    ws, fs, truths = freqs(kid, fs)
    print np.sort(truths)

    # plot sc superpgram
    print ("calculating superpgram for short cadence data")
    starts, stops, centres = real_footprint_sc(x)
    ivar = 1./yerr**2
    print np.median(yerr)
    print 1./np.sqrt(np.median(ivar))
    print np.median(y[1:]-y[:-1])
    print np.median(np.abs(np.diff(y[:10])))
    amp2s = spg.superpgram(starts, stops, y, ivar, ws)
    plt.clf()
    plt.subplot(3, 1, 1)
    uHz = 3589
    plt.plot(fs*1e6-uHz, amp2s, "k", alpha=.7, zorder=1)
    if len(truths):
        for truth in truths:
            plt.axvline(truth-uHz, color="r", linestyle="-.", zorder=0)
    plt.xlim(min(fs)*1e6-uHz, max(fs)*1e6-uHz)
#     ax = plt.gca()
#     ax.set_yticklabels([])
    plt.ylabel("$\mathrm{Short~cadence}$")
    plt.xlabel("$Frequency - %s ~(\mu Hz)$" % uHz)
#     plt.axvline(min(fs)*1e6-uHz)
#     plt.axvline(min(fs)*1e6-uHz + 8e-9)
#     plt.axvline(0.6, color=".5")
#     plt.axvline(0.6 + 8e-3, color=".5")

    # pwnnyquist long cadence data
    print ("calculating superpgram for long cadence data")
    x, y, yerr, ivar = load_data(kid, kepler_lc_dir, sc=False)
    starts, stops, centres = real_footprint(x)
    amp2s = spg.superpgram(starts, stops, y, ivar, ws)
    plt.subplot(3, 1, 2)
    plt.plot(fs*1e6-uHz, amp2s, "k", alpha=.7, zorder=1)
    if len(truths):
        for truth in truths:
            plt.axvline(truth-uHz, color="r", linestyle="-.", zorder=0)
    plt.xlim(min(fs)*1e6-uHz, max(fs)*1e6-uHz)
#     ax = plt.gca()
#     ax.set_yticklabels([])
    plt.ylim(0, 0.0000005)
    plt.ylabel("$\mathrm{Long~cadence}$")
    plt.xlabel("$Frequency - %s ~(\mu Hz)$" % uHz)
#     plt.axvline(0.6, color=".5")
#     plt.axvline(0.6 + 8e-3, color=".5")

    # plot light curve
    plt.subplot(3, 1, 3)
    plt.plot(x/24./3600, y, "k", alpha=.5)
    plt.xlabel("$\mathrm{Time~(days)}$")
    plt.ylabel("$\mathrm{Normalised~flux}$")
    plt.xlim(min(x/24/3600), max(x/24/3600))
    plt.subplots_adjust(hspace=.4)
    plt.savefig("%sfft" % kid)

#     import os
#     os.system('say "over"')
