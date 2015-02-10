import numpy as np
import matplotlib.pyplot as plt
from header_time import real_footprint_sc
from fft import corrected_data
import superpgram as spg
import emcee
from fft import load_data

plotpar = {'axes.labelsize': 12,
           'xtick.labelsize': 12,
           'ytick.labelsize': 12,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def soup(kid, DIR, c):

    # load data
    if c == "lc":
        x, y, yerr, ivar = load_data(kid, DIR, sc=False)
    elif c == "sc":
        x, y, yerr, ivar = corrected_data(kid, corrected_lc_dir)

    mean = np.mean(y)
    y -= mean
    fs = np.arange(0.003, 0.004, 4e-7)  # Hz 8006161
    ws = 2*np.pi*fs

    print ("calculating superpgram for %s data" % c)
    starts, stops, centres = real_footprint_sc(x)
    ivar = 1./yerr**2
    y = np.array([i.astype("float64") for i in y])
    ivar = np.array([i.astype("float64") for i in ivar])
    amp2s = spg.superpgram(starts, stops, y, ivar, ws)

    # plot superpgram
    plt.clf()
    uHz = 3589
    plt.plot(fs*1e6-uHz, amp2s, "k", alpha=.7, zorder=1)
    try:
        truths = np.genfromtxt("%s_freqs.txt" % kid).T  # uHz
        for truth in truths:
            plt.axvline(truth-uHz, color="r", linestyle="-.", zorder=0)
    except: pass
    plt.xlim(min(fs)*1e6-uHz, max(fs)*1e6-uHz)
    plt.ylabel("$\mathrm{Amp}^2$")
    plt.xlabel("$Frequency - %s ~(\mu Hz)$" % uHz)
    plt.savefig("%s_%s" % (kid, c))
    return fs, amp2s

def autocorr(fs, pgram, dnu, c):
    fs *= 1e6
    dnu *= 1e6
    acor =  emcee.autocorr.function(pgram)
    df = fs[1] - fs[0]
    print "calculating acf..."
    plt.clf()
    plt.plot(np.arange(len(acor))*df, acor, ".3", zorder=1)
    print dnu
    plt.axvline(dnu, color="r", linestyle="-.", zorder=0)
#     plt.ylim(-.05, .1)
    plt.ylim(-0.005, 0.01)
    plt.xlim(100, 200)
    plt.xlabel("$\mathrm{Delta~nu~(uHz)}$")
    plt.ylabel("$\mathrm{Correlation}$")
    plt.savefig("%sacf_%s" % (kid, c))

if __name__ == "__main__":

    # ask the question: "is there extra power at this frequency separation?"

    # calculate freqs at which you need to compute spectrum
    # try just one dnu
    dnu = 149.4*1e-6 # Hz
    fs1 = np.arange(0.003, 0.004, 4e-7)  # Hz 8006161
    ws1 = 2*np.pi*fs1
    fs2 = np.arange(0.003, 0.004, 4e-7)+dnu  # Hz 8006161
    ws2 = 2*np.pi*fs2

    # pwnnyquist short cadence data
    kid = "8006161"
    corrected_dir = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"
    kepler_dir = "/Users/angusr/.kplr/data/lightcurves/%s" % kid.zfill(9)
    fs, amp2s = soup(kid, kepler_dir)

#     x, y, yerr, ivar = corrected_data(kid, corrected_lc_dir)
#     mean = np.mean(y)
#     y -= mean
#     starts, stops, centres = real_footprint_sc(x)
#     ivar = 1./yerr**2

#     print ("calculating superpgram for fs1")
#     amp2s1 = spg.superpgram(starts, stops, y, ivar, ws1)
#     autocorr(fs1, amp2s1, dnu)

    autocorr(fs, amp2s, dnu, c)

#     print ("calculating superpgram for fs")
#     amp2s2 = spg.superpgram(starts, stops, y, ivar, ws2)
#
#     plt.clf()
#     plt.plot(fs1, amp2s1+amp2s2, "k")
#     plt.savefig("test")
