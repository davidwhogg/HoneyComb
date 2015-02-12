import numpy as np
import matplotlib.pyplot as plt
from header_time import real_footprint_sc, real_footprint, bjd2utc
from fft import corrected_data, load_data
import superpgram as spg
import emcee

plotpar = {'axes.labelsize': 12,
           'xtick.labelsize': 12,
           'ytick.labelsize': 12,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def soup(kid, nm, DIR, c, KDIR, plot=False):

    # load data (x is in seconds)
    if c == "lc":
        x, y, yerr, ivar = load_data(kid, KDIR, sc=False)
    elif c == "sc":
        try:
            x, y, yerr, ivar = corrected_data(kid, DIR)
        except:
            x, y, yerr, ivar = load_data(kid, DIR, sc=True)

    # subtract mean
    mean = np.mean(y)
    y -= mean
#     fs = np.arange(0.003, 0.004, 4e-7)  # Hz 8006161
    fs = np.arange(nm-.2*nm, nm+.2*nm, 4e-8)
#     fs = np.arange(nm-.0001, nm+.0001, 4e-8)
#     fs = np.arange(nm-.1*nm, nm+.1*nm, 4e-7)
    print len(fs), "frequencies"
    ws = 2*np.pi*fs

    print ("calculating superpgram for %s data" % c)
    if c == "sc":
        starts, stops, centres = real_footprint_sc(x)
    elif c == "lc":
        starts, stops, centres = real_footprint(x)
    ivar = 1./yerr**2
    y = np.array([i.astype("float64") for i in y])
    ivar = np.array([i.astype("float64") for i in ivar])

#     plt.clf()
#     plt.plot(x, y, "k.")
#     plt.savefig("%s_data" % kid)

    amp2s = spg.superpgram(starts, stops, y, ivar, ws)

    # plot superpgram
    if plot == True:
        plt.clf()
        uHz = nm*1e6
        plt.plot(fs*1e6-uHz, amp2s, "k", alpha=.7, zorder=1)
        try:
            truths = np.genfromtxt("%s_freqs.txt" % kid).T  # uHz
            print "frequency file found"
            print truths
            for truth in truths:
                plt.axvline(truth-uHz, color="r", linestyle="-.", zorder=0)
        except:
            print "no frequency file"
            pass
        plt.xlim(min(fs)*1e6-uHz, max(fs)*1e6-uHz)
        plt.ylabel("$\mathrm{Amp}^2$")
        plt.xlabel("$Frequency - %s ~(\mu Hz)$" % uHz)
        plt.savefig("%spgram_%s" % (kid, c))
    return fs, amp2s

def autocorr(kid, fs, pgram, dnu, c, plot=False):
    fs *= 1e6  # convert Hz to uHz
    dnu *= 1e6
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
        plt.xlabel("$\mathrm{Delta~nu~(uHz)}$")
        plt.ylabel("$\mathrm{Correlation}$")
        plt.savefig("%sacf_%s" % (kid, c))
    return lags, acor

if __name__ == "__main__":

    D = "/Users/angusr/Python/HoneyComb/pwnnyquist"
    kids, nm, dnu = np.genfromtxt("%s/target_list.txt" % D, skip_header=1).T

    for i in range(len(kids)):
        print "kid = ", str(int(kids[i])), "nm = ", nm[i], "dnu = ", dnu[i]

        c = "sc"  # long or short cadence data?
        DIR = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"
        KDIR = "/Users/angusr/.kplr/data/lightcurves/%s" \
                % str(int(kids[i])).zfill(9)

        try:
            # compute supergram
            fs, amp2s = soup(str(int(kids[i])), nm[i]*1e-6, DIR, c, KDIR)

            # compute autocorrelation
            autocorr(str(int(kids[i])), fs, amp2s, dnu[i]*1e-6, c)
        except:
            "LinAlgError"
            print "No data found"
