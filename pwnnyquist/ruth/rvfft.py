import numpy as np
import matplotlib.pyplot as plt
from params import plot_params
reb = plot_params()
from colours import plot_colours
c = plot_colours()
import superpgram as spg

def load_rvs():
    x, y, yerr = np.genfromtxt("K00975.HARPSN.vels").T
    ivar = 1./yerr**2
    return x, y, yerr, ivar

def rv_footprint(x):
    data = np.genfromtxt("KOI975.HARPSN.info.txt", skip_header=1).T
    exptime = data[15]
    stops = x
    starts = x - exptime
    centres = x - .5*exptime
    return starts, stops, centres

if __name__ == "__main__":

    x, y, yerr, ivar = load_rvs()
    x -= x[0]
    starts, stops, centres = rv_footprint(x)
    y -= np.mean(y)

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.errorbar(x, y, yerr=yerr, **reb)

    plt.subplot(2, 1, 2)
    fs = np.linspace(0.00001, 0.003, 1000)  # Hz
    ws = 2*np.pi*fs
    y = np.array([i.astype("float64") for i in y])
    ivar = np.array([i.astype("float64") for i in ivar])

    y = y.copy(order='C')
    starts = starts.copy(order='C')
    stops = stops.copy(order='C')
    ivar = ivar.copy(order='C')
    ws = ws.copy(order='C')

    amp2s = spg.superpgram(starts, stops, y, ivar, ws)
    plt.plot(fs, amp2s, "k")

    truths = np.genfromtxt("3632418_freqs.txt").T  # uHz
    for truth in truths:
        plt.axvline(truth*1e-6, color=c.blue, linestyle="--")
    plt.savefig("rv")
