import numpy as np
import matplotlib.pyplot as plt
from header_time import real_footprint_sc, real_footprint, bjd2utc
from fft import corrected_data
import superpgram as spg
import emcee
from astero import astero

plotpar = {'axes.labelsize': 12,
           'xtick.labelsize': 12,
           'ytick.labelsize': 12,
           'text.usetex': True}
plt.rcParams.update(plotpar)

# run superpgram on a star
def soup(x, y, ivar, fs, nm, kid, c, plot=False):

    # subtract mean
    mean = np.mean(y)
    y -= mean

    # format data
    y = np.array([_y.astype("float64") for _y in y])
    ivar = np.array([_ivar.astype("float64") for _ivar in ivar])

    ws = 2*np.pi*fs
    amp2s = spg.superpgram(starts, stops, y, ivar, ws)

    # plot superpgram
    if plot == True:
        print "plotting"
        plt.clf()
        plt.plot(fs*1e6, amp2s, "k", alpha=.7, zorder=1)
        print fs[0], fs*1e6
        try:
            truths = np.genfromtxt("%s_freqs.txt" % kid).T  # uHz
            print "frequency file found"
            print truths
            print max(truths) - min(truths)
            for truth in truths:
                plt.axvline(truth, color="r", linestyle="-.", zorder=0)
        except:
            print "no frequency file"
            pass
#         plt.xlim(min(fs)*1e6-uHz, max(fs)*1e6-uHz)
        plt.ylabel("$\mathrm{Amp}^2$")
        plt.xlabel("$Frequency (\mu Hz)$")
        plt.ylim(0, max(amp2s))
        plt.savefig("%spgram_%s" % (str(int(kid)), c))

    np.savetxt("results/%spgram_%s.txt" % (str(int(kid)), c),
               np.transpose((fs, amp2s)))
    return amp2s

# compute the autocorrelation function of the pgram
def autocorr(kid, dnu, fmin, fmax, c, plot=False):

    fs, pgram = np.genfromtxt("results/%spgram_%s.txt" % (str(int(kid)), c)).T
    fs *= 1e6  # convert Hz to uHz
    dnu *= 1e6

    l = (fmin < fs) * (fs < fmax)
    fs, pgram = fs[l], pgram[l]

    acor =  emcee.autocorr.function(pgram)
    df = fs[1] - fs[0]

    print "calculating acf..."
    lags = np.arange(len(acor))*df
    print "delta_nu = ", dnu
    if plot == True:
        plt.clf()
        plt.plot(lags, acor, ".3", zorder=1)
        plt.axvline(dnu, color="r", linestyle="-.", zorder=0)
        plt.xlim(dnu-.5*dnu, dnu+.5*dnu)
        plt.ylim(0, .2)
        plt.xlabel("$\mathrm{Delta~nu~(uHz)}$")
        plt.ylabel("$\mathrm{Correlation}$")
        plt.savefig("%sacf_%s" % (str(int(kid)), c))
    return lags, acor

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_best_acf_section(fmins, fmaxs, dnu, kid):
    print "computing autocorrelation function"
    col = np.linspace(.05, 1, len(fmins))[::-1]
    plt.clf()
    print fmins
    print fmaxs
    for i in range(len(fmins)):
        lags, acor = autocorr(kid, dnu*1e-6, fmins[i], fmaxs[::-1][i], "sc")
        plt.plot(lags, acor, color="%s" % col[i])
    plt.plot(lags, acor, color="%s" % col[i])
    plt.axvline(dnu, color="r", linestyle="-.", zorder=0)
    plt.xlim(dnu-.5*dnu, dnu+.5*dnu)
    plt.ylim(0, .2)
    plt.xlabel("$\mathrm{Delta~nu~(uHz)}$")
    plt.ylabel("$\mathrm{Correlation}$")
    plt.savefig("%sacf_%s" % (str(int(kid)), "combo"))

if __name__ == "__main__":

    D = "/Users/angusr/Python/HoneyComb/pwnnyquist"
    DIR = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"
    kid = "8006161"
    nm = 3481
    dnu = 149.3
#     fs = np.arange(nm - dnu*5, nm + dnu*5, 4e-2) * 1e-6
    fs = np.arange(100, 4000, 4e-1) * 1e-6
    fmins = np.arange(2000, 3500, 500)
    fmaxs = np.arange(4000, 5500, 500)

    # Giant poster child
#     kid = "3424541"
#     nm = 712
#     dnu = 41.3
#     fs = np.arange(300, 1200, 4e-1) * 1e-6
#     fmins = np.arange(100, 600, 100)
#     fmaxs = np.arange(700, 1200, 100)

    # load KASOC light curve
    x, y, yerr, ivar = corrected_data(str(int(kid)), DIR)
    starts, stops, centres = real_footprint_sc(x)

    print dnu, "dnu"
    print nm
    print len(fs), "frequencies"
    print "computing supergram..."
#     try:
    fs, amp2s = np.genfromtxt("results/8006161pgram_sc.txt").T
    print "periodogram file found"
    print "plotting"
    plt.clf()
    plt.plot(fs*1e6, amp2s, "k", alpha=.7, zorder=1)
    print fs[0], fs*1e6
    truths = np.genfromtxt("%s_freqs.txt" % kid).T  # uHz
    print "frequency file found"
    print truths
    print max(truths) - min(truths)
    for truth in truths:
        plt.axvline(truth, color="r", linestyle="-.", zorder=0)
#      plt.xlim(100, max(fs)*1e6-uHz)
    plt.ylabel("$\mathrm{Amp}^2$")
    plt.xlabel("$Frequency (\mu Hz)$")
    plt.ylim(0, max(amp2s[200:]))
    plt.savefig("8006161_test")
#     except:
#         print "no periodogram file"
#         amp2s = soup(x, y, ivar, fs, nm*1e-6, kid, "sc", plot=True)

    print "computing autocorrelation function"
    # compute over a certain section of the pgram
    plt.clf()
    lags, acor = autocorr(kid, dnu*1e-6, fmins[0], fmaxs[-1], "sc")
    plt.plot(lags, acor, "k")
    plt.axvline(dnu, color="r", linestyle="-.", zorder=0)
    plt.xlim(dnu-.5*dnu, dnu+.5*dnu)
    plt.ylim(0, .2)
    plt.xlabel("$\mathrm{Delta~nu~(uHz)}$")
    plt.ylabel("$\mathrm{Correlation}$")
    plt.savefig("%sacf_%s" % (str(int(kid)), "sc"))

#     find_best_acf_section(fmins, fmaxs, dnu, kid)
