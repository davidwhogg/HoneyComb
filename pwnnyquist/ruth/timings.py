import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
from fft import load_data, corrected_data
from header_time import bjd2utc

kid = "8006161"
dnu = 149.4*1e-6 # Hz
c = "lc"  # long or short cadence data?
DIR = "/Users/angusr/Python/HoneyComb/pwnnyquist/KASOC"
if c == "lc":
    DIR = "/Users/angusr/.kplr/data/lightcurves/%s" % kid.zfill(9)

# load data
if c == "lc":
    x, y, yerr, ivar = load_data(kid, DIR, sc=False)
elif c == "sc":
    x, y, yerr, ivar = corrected_data(kid, DIR)

# plot bjds
plt.clf()
plt.plot(np.diff(x)[np.diff(x) < 2000])
plt.savefig("bjd")

utc = bjd2utc(x/24/3600)
utc *= 24*3600
plt.clf()
plt.plot(np.diff(utc)[np.diff(utc) < 2000])
plt.savefig("utc")
