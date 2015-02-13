import numpy as np
import matplotlib.pyplot as plt
from delta_nu import soup, autocorr
from astero import astero
data = astero()
import h5py

# search for delta nu in the long cadence data of Chaplin et al.
D = "/n/home11/rangus/Python/HoneyComb/pwnnyquist/data"
# D = "/Users/angusr/Python/Gyro/data"
kids = data.iKID
nm = data.inu_max
dnu = data.idnu

for i in range(len(kids)):
    print "kid = ", str(int(kids[i])), "nm = ", nm[i], "dnu = ", dnu[i]

    c = "lc"  # long or short cadence data?
    DIR = "/n/home11/rangus/.kplr/data/lightcurves/%s" \
            % str(int(kids[i])).zfill(9)
#     DIR = "/Users/angusr/angusr/data2/all_Qs"
    KDIR = DIR

#     try:
    # compute supergram
    fs, amp2s = soup(str(int(kids[i])), nm[i]*1e-6, DIR, c, KDIR, plot=True)
    f = h5py.File("%spgram.h5" % str(int(kids[i])), "w")
    data = f.create_dataset("pgram", (len(fs), 2))
    data[:, 0] = np.array(fs)
    data[:, 1] = np.array(amp2s)
    f.close()

    # compute autocorrelation
    lags, acf = autocorr(str(int(kids[i])), fs, amp2s, dnu[i]*1e-6, c,
                         plot=True)
    f = h5py.File("%sacf.h5" % str(int(kids[i])), "w")
    data = f.create_dataset("ACF", (len(lags), 2))
    data[:, 0] = np.array(lags)
    data[:, 1] = np.array(acf)
    f.close()
#     except:
#         "LinAlgError"
#         print "No data found"
#         assert 0
