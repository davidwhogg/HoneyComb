import numpy as np
import matplotlib.pyplot as plt
from delta_nu import soup, autocorr
from astero import astero
data = astero()

# search for delta nu in the long cadence data of Chaplin et al.
D = "/Users/angusr/Python/Gyro/data"
kids = data.iKID
nm = data.inu_max
dnu = data.idnu

for i in range(len(kids)):
    print "kid = ", str(int(kids[i])), "nm = ", nm[i], "dnu = ", dnu[i]

    c = "lc"  # long or short cadence data?
    DIR = "/Users/angusr/angusr/data2/all_Qs"
    KDIR = DIR

    try:
        # compute supergram
        fs, amp2s = soup(str(int(kids[i])), nm[i]*1e-6, DIR, c, KDIR)

        # compute autocorrelation
        autocorr(str(int(kids[i])), fs, amp2s, dnu[i]*1e-6, c)
        assert 0
    except:
        "LinAlgError"
        print "No data found"
