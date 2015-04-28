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

# run superpgram on a star
def soup(x, y, ivar, fs, nm, kid, plot=False):

    amp2s = spg.superpgram(starts, stops, y, ivar, ws)

    # plot superpgram
    if plot == True:
        print "plotting"
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
    return amp2s

# compute the autocorrelation function of the pgram
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

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

if __name__ == "__main__":

    D = "/Users/angusr/Python/HoneyComb/pwnnyquist"
#     D = "/n/home11/rangus/Python/HoneyComb/pwnnyquist"
    kids, nm, dnu = np.genfromtxt("%s/target_list.txt" % D, skip_header=1).T

    for i, kid in enumerate(kids):
        print "kid = ", str(int(kids[i])), "nm = ", nm[i], "dnu = ", dnu[i]

        c = "sc"  # long or short cadence data?
        DIR = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"
        KDIR = "/Users/angusr/.kplr/data/lightcurves/%s" \
                % str(int(kids[i])).zfill(9)
#         KDIR = "/n/home11/rangus/.kplr/data/lightcurves/%s" \
#                 % str(int(kids[i])).zfill(9)
#         DIR = KDIR

        # load data (x is in seconds)
        if c == "lc":
            x, y, yerr, ivar = load_data(str(int(kid)), KDIR, sc=False)
            starts, stops, centres = real_footprint(x)
        elif c == "sc":
            x, y, yerr, ivar = load_data(str(int(kid)), KDIR, sc=True)
            starts, stops, centres = real_footprint_sc(x)

        # subtract mean
        mean = np.mean(y)
        y -= mean

        # format data
        ivar = 1./yerr**2
        y = np.array([_y.astype("float64") for _y in y])
        ivar = np.array([_ivar.astype("float64") for _ivar in ivar])

        # compute supergram
        fs = np.arange(nm[i]-.5*nm[i], nm[i]+.5*nm[i], 4e-8)
        print len(fs), "frequencies"
        ws = 2*np.pi*fs
        amp2s = soup(x, y, ivar, ws, nm[i]*1e-6, str(int(kids[i])),
                         plot=True)

        # compute autocorrelation
        autocorr(str(int(kids[i])), fs, amp2s, dnu[i]*1e-6, c, plot=True)
        assert 0
