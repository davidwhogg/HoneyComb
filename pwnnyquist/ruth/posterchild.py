import numpy as np
import matplotlib.pyplot as plt
from fft_old import load_data, corrected_data
import superpgram as spg
from header_time import real_footprint_sc
from plotstuff import params
reb = params()
from ruth_binning import gap_split, remainder, rebin

def binned_fp(t, n):
    dt = 6.81119940564e-4*24*3600*n # interval between observations (seconds)
    stops = t
    starts = t - dt
    centres = t - .5*dt
    return starts, stops, centres, dt

def inject_sine_wave(x, y, f):
    return y+10000*np.sin(2*np.pi*f*x)

def binned_spg(kid, fs, nq, n, bins=True,inject=True):
    """
    Calculate superpgram of KASOC corrected data for kid, over the frequencies,
    fs, with nq months of data and bins containing n data points.
    """
    # load corrected data
    D = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"
    x, y, yerr, ivar = corrected_data(kid, D, 1)  # median to mean

    print "calculating spg over ", len(fs), "frequencies"
    freqs = np.genfromtxt("data/%s_freqs.txt" % kid).T
    truth = freqs[(freqs > min(fs)*1e6) * (freqs < max(fs)*1e6)]
    print "should see the following modes: ", truth

    # inject a sine wave
    y = inject_sine_wave(x, y, truth[0]*1e-6)

    if bins:
        # split data at gaps
        xs, ys, yerrs = gap_split(x, y, yerr)

        # rebin the data (still split at gaps)
        xt = [rebin(xs[i], ys[i], yerrs[i], n)[0] for i in range(len(xs))]
        yt = [rebin(xs[i], ys[i], yerrs[i], n)[1] for i in range(len(xs))]
        yerrt = [rebin(xs[i], ys[i], yerrs[i], n)[2] for i in range(len(xs))]

        # reassemble the data
        xx = np.array([i for j in xt for i in j])
        yy = np.array([i for j in yt for i in j])
        yyerr = np.array([i for j in yerrt for i in j])
    else: xx, yy, yyerr = x, y, yerr

    print x[1] - x[0]
    print xx[1] - xx[0]
    starts, stops, centres, dt = binned_fp(xx, n)
    print dt
    print "Nyquist = ", (2./dt) * 1e6, "uHz"

    # run the spg
    ivar = 1./yyerr**2
    print "running spg..."
    amp2s = spg.superpgram(starts, stops, yy, ivar, fs*2*np.pi)

    # plot
    print "plotting..."
    plt.clf()
    plt.plot(fs*1e6, amp2s, "k")
    for i in truth:
        plt.axvline(i, color="r", ls="--")
    plt.xlabel("$\mathrm{Frequency~[} \\mu \mathrm{Hz]}$")
    if bins:
        plt.savefig("pc_%s" % n)
    else:
        plt.savefig("pc")

if __name__ == "__main__":
    kid = "8006161"
    fs = np.arange(3540, 3640, 1e-1)*1e-6
    binned_spg(kid, fs, 1, 3)
