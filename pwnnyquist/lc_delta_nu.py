import numpy as np
import matplotlib.pyplot as plt
from delta_nu import soup, autocorr
from astero import astero
data = astero()

# search for delta nu in the long cadence data of Chaplin et al.
D = "/n/home11/rangus/Python/HoneyComb/pwnnyquist/data"
kids = data.iKID
nm = data.inu_max
dnu = data.idnu

for i in range(len(kids)):
    print "kid = ", str(int(kids[i])), "nm = ", nm[i], "dnu = ", dnu[i]

    c = "lc"  # long or short cadence data?
    DIR = "/n/home11/rangus/.kplr/data/lightcurves/%s" \
            % str(int(kids[i])).zfill(9)
    KDIR = DIR

    try:
        # compute supergram
        fs, amp2s = soup(str(int(kids[i])), nm[i]*1e-6, DIR, c, KDIR)
        np.savetxt("%sspg.txt" np.transpose((fs, amp2s)))

        # compute autocorrelation
        lags, acf = autocorr(str(int(kids[i])), fs, amp2s, dnu[i]*1e-6, c)
        np.savetxt("%sacf.txt" np.transpose((lags, acf)))
        assert 0
    except:
        "LinAlgError"
        print "No data found"
        assert 0
