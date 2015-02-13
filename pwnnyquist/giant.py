import numpy as np
import matplotlib.pyplot as plt
import pyfits
import nufft
import superpgram as spg
from header_time import real_footprint, real_footprint_sc
import glob
from gatspy.periodic import LombScargle

plotpar = {'axes.labelsize': 8,
           'xtick.labelsize': 8,
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

# def triangular_smoothing(x, y, yerr):

if __name__ == "__main__":

    kid = "8006161"
#     kid = "5515314"
#     kid = "7103006"
#     kid = "3427720"
#     kid = "3632418"
#     kid = "7341231"
    kepler_lc_dir = "/Users/angusr/.kplr/data/lightcurves/%s" % kid.zfill(9)
    corrected_lc_dir = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"

#     # pwnnyquist short cadence data
    x, y, yerr, ivar = corrected_data(kid, corrected_lc_dir)
# #     x, y, yerr, ivar = load_data(kid, kepler_lc_dir, sc=True)

    fs0 = .0001375  # 5515
    fs0 = .003589  # 8006161
#     fs = np.arange(0.0001372, 0.00013752, 1e-10)  # Hz 5515314
    fs = np.arange(fs0-.000001, fs0+0.000001, 1e-9)  # Hz 8006161
    fs = np.arange(fs0-.00001, fs0+0.00001, 1e-8)  # Hz 8006161
    fs = np.arange(fs0-.0005, fs0+0.0005, 1e-6)  # Hz 8006161
    print len(fs)
    ws, fs, truths = freqs(kid, fs)

#     # plot sc superpgram
#     print ("calculating superpgram for short cadence data")
#     starts, stops, centres = real_footprint_sc(x)
#     print len(ws)
#     amp2s = spg.superpgram(starts, stops, y, ivar, ws)
#     plt.clf()
#     plt.subplot(3, 1, 1)
#     plt.plot((fs-fs0)*1e6, amp2s, "k", alpha=.7)
#     if len(truths):
#         for truth in truths:
#             print truth
#             plt.axvline(truth*1e-6, color="r", linestyle="-.")
#     plt.xlim(min((fs-fs0)*1e6), max((fs-fs0)*1e6))
#     plt.ylabel("$\mathrm{Short~cadence}$")
#     plt.xlabel("$\mathrm{Frequency-%s~(uHz)}$" % (fs0*1e6))

#     # fft long cadence data
#     plt.subplot(3, 1, 2)
#     print ("calculating  fft for short cadence data")
#     print len(ws)
#     pgram = nufft.nufft3(stops, y, ws)
#     plt.plot((fs-fs0)*1e6, pgram, "k", alpha=.7)
#     if len(truths):
#         for truth in truths:
#             plt.axvline(truth*1e-6, color="r", linestyle="-.")
#     plt.xlim(min((fs-fs0)*1e6), max((fs-fs0)*1e6))
#     plt.ylabel("$\mathrm{FFT}$")
#     plt.xlabel("$\mathrm{Frequency-%s~(uHz)}$" % (fs0*1e6))

#     # lombscargle short cadence data
#     print ("calculating Lomb-Scargle for short cadence data")
#     model = LombScargle().fit(x, y, yerr)
#     period = 1. / fs
#     power = model.periodogram(period)
#     plt.subplot(3, 1, 3)
#     plt.plot(((fs)-fs0)*1e6, power, "k", alpha=.7)
#     if len(truths):
#         for truth in truths:
#             plt.axvline(truth*1e-6, color="r", linestyle="-.")
#     plt.xlim(min((fs-fs0)*1e6), max((fs-fs0)*1e6))
#     plt.ylabel("$\mathrm{sc~ls}$")
#     plt.xlabel("$\mathrm{Frequency-%s~(uHz)}$" % (fs0*1e6))
#     plt.subplots_adjust(hspace=.3)

    # pwnnyquist long cadence data
    plt.clf()
    plt.subplot(3, 1, 1)
    print ("calculating superpgram for long cadence data")
    x, y, yerr, ivar = load_data(kid, kepler_lc_dir, sc=False)
    starts, stops, centres = real_footprint(x)
    amp2s = spg.superpgram(starts, stops, y, ivar, ws)
    plt.plot((fs-fs0)*1e6, amp2s, "k", alpha=.7)
    if len(truths):
        for truth in truths:
            plt.axvline(truth*1e-6, color="r", linestyle="-.")
    plt.ylabel("$\mathrm{SuperPgram}$")
    plt.xlabel("$\mathrm{Frequency-%s~(uHz)}$" % (fs0*1e6))

    # fft long cadence data
    plt.subplot(3, 1, 2)
    print ("calculating  fft for long cadence data")
    pgram = nufft.nufft3(stops, y, ws)
    plt.plot((fs-fs0)*1e6, pgram, "k", alpha=.7)
    if len(truths):
        for truth in truths:
            plt.axvline(truth*1e-6, color="r", linestyle="-.")
    plt.ylabel("$\mathrm{FFT}$")
    plt.xlabel("$\mathrm{Frequency-%s~(uHz)}$" % (fs0*1e6))

    # LombScargle long cadence data
    plt.subplot(3, 1, 3)
    print ("calculating LombScargle for long cadence data")
    periods = 1./fs
    model = LombScargle().fit(x, y, yerr)
    power = model.periodogram(periods)
    plt.plot((fs-fs0)*1e6, power, "k", alpha=.7)
    if len(truths):
        for truth in truths:
            plt.axvline(truth*1e-6, color="r", linestyle="-.")
    plt.ylabel("$\mathrm{LS}$")
    plt.xlabel("$\mathrm{Frequency-%s~(uHz)}$" % (fs0*1e6))
    plt.subplots_adjust(hspace=.4)

#     # pwnnyquist long cadence data
#     plt.clf()
#     plt.subplot(3, 1, 1)
#     print ("calculating superpgram for long cadence data")
#     x, y, yerr, ivar = load_data(kid, kepler_lc_dir, sc=False)
#     starts, stops, centres = real_footprint(x)
#     amp2s = spg.superpgram(starts, stops, y, ivar, ws)
#     plt.plot((fs-fs0)*1e6, amp2s, "k", alpha=.7)
#     if len(truths):
#         for truth in truths:
#             plt.axvline(truth*1e-6, color=".3", linestyle="-.")
#     plt.ylabel("$\mathrm{SuperPgram}$")
#     plt.xlabel("$\mathrm{Frequency-135~(uHz)}$")
#
#     # fft long cadence data
#     plt.subplot(3, 1, 2)
#     print ("calculating  fft for long cadence data")
#     pgram = nufft.nufft3(stops, y, ws)
#     plt.plot((fs-fs0)*1e6, pgram, "k", alpha=.7)
#     if len(truths):
#         for truth in truths:
#             plt.axvline(truth*1e-6, color=".3", linestyle="-.")
#     plt.ylabel("$\mathrm{FFT}$")
#     plt.xlabel("$\mathrm{Frequency-135~(uHz)}$")
#
#     # LombScargle long cadence data
#     plt.subplot(3, 1, 3)
#     print ("calculating LombScargle for long cadence data")
#     periods = 1./fs
#     model = LombScargle().fit(x, y, yerr)
#     power = model.periodogram(periods)
#     plt.plot((fs-fs0)*1e6, power, "k", alpha=.7)
#     if len(truths):
#         for truth in truths:
#             plt.axvline(truth*1e-6, color=".3", linestyle="-.")
#     plt.ylabel("$\mathrm{LS}$")
#     plt.xlabel("$\mathrm{Frequency-135~(uHz)}$")
#     plt.subplots_adjust(hspace=.4)

    # plot light curve
#     plt.subplot(4, 1, 4)
# #     plt.subplot(3, 1, 3)
#     plt.plot(x/24./3600, y, "k.", alpha=.5)
#     plt.xlabel("$\mathrm{Time~(days)}$")
#     plt.ylabel("$\mathrm{Normalised~flux}$")
#     plt.xlim(min(x/24/3600), max(x/24/3600))
    plt.savefig("giant2")
