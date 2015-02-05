import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
import nufft
from params import plot_params
reb = plot_params()
import superpgram as spg
from colours import plot_colours
c = plot_colours()
import scipy.signal as sps
import superpgram as spg
from header_time import real_footprint

# load sc data
kid = "8006161"  # a kepler target with lc and sc data chosen at `random'
# kid = "7103006"
# kid = "3427720"
D = "/Users/angusr/.kplr/data/lightcurves/%s" % kid.zfill(9)
fnames = glob.glob("%s/kplr%s-*_slc.fits" % (D, kid.zfill(9)))

q, x, y, yerr = [], [], [], []
for fname in fnames:
    hdulist = pyfits.open(fname)
    tbdata = hdulist[1].data
    x.extend(tbdata["TIME"])
    data = tbdata["PDCSAP_FLUX"]
    q.extend(tbdata["SAP_QUALITY"])
    med = np.median(data)
    y.extend(data/med-1)
    abs_yerr = tbdata["PDCSAP_FLUX_ERR"]
    yerr.extend(abs_yerr/med)
x, y, yerr, q = np.array(x), np.array(y), np.array(yerr), np.array(q)
l = np.isfinite(x) * np.isfinite(y) * np.isfinite(yerr) * (q==0)
x, y, yerr = x[l], y[l], yerr[l]

x *= 24*3600  # convert to seconds

truths = np.genfromtxt("%s_freqs.txt" % kid).T  # uHz
fs = np.linspace(0.0025, 0.004, 1000)  # Hz
print min(truths*1e-6), max(truths*1e-6)
ws = 2*np.pi*fs
pgram = nufft.nufft3(x, y, ws)

plt.clf()
plt.subplot(3, 1, 1)
plt.plot(x, y, "k.", alpha=.1)

plt.subplot(3, 1, 2)
p = np.sqrt(np.real(pgram)**2 + np.imag(pgram)**2)
plt.plot(fs, p, "k")
for truth in truths:
    plt.axvline(truth*1e-6, color=c.blue, linestyle="--")
plt.xlim(min(fs), max(fs))

ivar = 1./yerr**2
y = np.array([i.astype("float64") for i in y])
ivar = np.array([i.astype("float64") for i in ivar])
starts, stops, centres = real_footprint(x)
testks = ws
amp2s = spg.superpgram(starts, stops, y, ivar, testks)
plt.subplot(3, 1, 3)
for truth in truths:
    plt.axvline(truth*1e-6, color=c.blue, linestyle="--")
plt.plot(fs, amp2s, "k")
plt.xlim(min(fs), max(fs))
plt.savefig("fft")
