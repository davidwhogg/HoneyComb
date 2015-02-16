import numpy as np
import matplotlib.pyplot as plt
from delta_nu import soup, autocorr
from astero import astero
import h5py
import sys

# search for delta nu in the long cadence data of Chaplin et al.
def dn_search(nmin, nmax, c):

    # astero stellar param data directory
    D = "/n/home11/rangus/Python/HoneyComb/pwnnyquist/data"
#     D = "/Users/angusr/Python/Gyro/data"

    data = astero()
    kids = data.iKID
    nm = data.inu_max
    dnu = data.idnu
    l = (nm < nmax) * (nm > nmin)
    kids, nm, dnu = kids[l], nm[l], dnu[l]

    for i in range(len(kids)):
        print "kid = ", str(int(kids[i])), "nm = ", nm[i], "dnu = ", dnu[i]

        # light curve directory
        DIR = "/n/home11/rangus/.kplr/data/lightcurves/%s" \
                % str(int(kids[i])).zfill(9)
#         DIR = "/Users/angusr/angusr/data2/all_Qs"

        c = "lc"  # long or short cadence data?
        KDIR = DIR

        # compute supergram
        fs, amp2s = soup(str(int(kids[i])), nm[i]*1e-6, DIR, c, KDIR)
        f = h5py.File("%spgram_%s.h5" % (str(int(kids[i])), c), "w")
        data = f.create_dataset("pgram", (len(fs), 2))
        data[:, 0] = np.array(fs)
        data[:, 1] = np.array(amp2s)
        f.close()

        # compute autocorrelation
        lags, acf = autocorr(str(int(kids[i])), fs, amp2s, dnu[i]*1e-6, c)
        f = h5py.File("%sacf_%s.h5" % (str(int(kids[i])), c), "w")
        data = f.create_dataset("ACF", (len(lags), 2))
        data[:, 0] = np.array(lags)
        data[:, 1] = np.array(acf)
        f.close()

if __name__ == "__main__":

    nmin = float(sys.argv[1])
    nmax = float(sys.argv[2])
    c = str(sys.argv[3])
#     nmin, nmax, c = 0, 5000, "lc"
    dn_search(nmin, nmax, c)
