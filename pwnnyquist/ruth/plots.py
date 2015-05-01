import numpy as np
import matplotlib.pyplot as plt
import glob
from params import plot_params
reb = plot_params()
from astero import astero
data = astero()
import h5py
import glob

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def make_plots(D, kids, dnu, nm, c):

    nacf = 30  # 1/2 times the number of lags (just plot the acf around dn)
    acf_sum = np.zeros(2*nacf)

    for i, kid in enumerate(kids):

        # load the peridogram
#         with h5py.File("%s/results/%spgram_%s.h5"
#                        % (D, str(int(kid))), c) as f:
        with h5py.File("%s/results/%spgram.h5" % (D, str(int(kid)))) as f:
            try:
                fs = f["pgram"][:, 0]
                amp2s = f["pgram"][:, 1]
            except:
                fs = np.zeros(1000)
                amp2s = np.zeros(1000)

        # plot the periodogram
        plt.clf()
        uHz = nm[i]
        plt.plot((fs*1e6), amp2s, "k", alpha=.7, zorder=1)
        plt.axvline(uHz, color="r", linestyle="--")
        plt.xlim(min(fs)*1e6, max(fs)*1e6)
        plt.ylabel("$\mathrm{Amp}^2$")
        plt.xlabel("$Frequency~(\mu Hz)~%s$" % uHz)
        plt.savefig("%s/figs/%spgram_%s" % (D, str(int(kid)), c))

    # read acf results
#         with h5py.File("%s/results/%sacf_%s.h5" % (D, str(int(kid))), c) as f:
        with h5py.File("%s/results/%sacf.h5" % (D, str(int(kid)))) as f:
            try:
                lags = f["ACF"][:, 0]
                acf = f["ACF"][:, 1]
            except:
                lags = np.zeros_like(acf_sum)
                acf = np.zeros_like(acf_sum)

        #  plot the acf
        print kid, dnu[i], nm[i]
        print len(acf), len(fs)

        df = (fs[1] - fs[0])*1e6
        lags = np.arange(len(acf))*df
        plt.clf()
        plt.plot(lags, acf, "k")
        plt.axvline(dnu[i], color="r", linestyle="--")
#         plt.xlim(dnu[i]-.5*dnu[i], dnu[i]+.5*dnu[i])
        plt.xlabel("$\mathrm{Delta~nu~(uHz)~%s}$" % dnu[i])
        plt.ylabel("$\mathrm{Correlation}$")
        plt.savefig("%s/figs/%sacf_%s" % (D, str(int(kid)), c))

        # add the acf to the total
        nearest_dnu = find_nearest(lags, dnu[i])
        k = np.where(lags==nearest_dnu)[0]
        print k
        acf_sum += acf[k-nacf:k+nacf]

    plt.clf()
    plt.plot(acf_sum, "k")
    plt.savefig("acf_sum_%s" % c)
    np.savetxt("acf_sum_%s.txt" % c, np.transpose((acf_sum)))

if __name__ == "__main__":

    D = "/Users/angusr/Python/HoneyComb/pwnnyquist"
    kids = data.iKID
    dnu = data.idnu
    nm = data.inu_max

    nm_min, nm_max = 0, 300
    l = (nm < nm_max) * (nm > nm_min)
    kids, nm, dnu = kids[l], nm[l], dnu[l]
    print len(kids)

    make_plots(D, kids, dnu, nm, "sc")
