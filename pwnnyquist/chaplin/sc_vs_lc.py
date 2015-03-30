import numpy as np
import matplotlib.pyplot as plt
import nufft
from fft import load_data
import superpgram as spg
from header_time import real_footprint, real_footprint_sc

name = "6442183"
print "loading data..."
D = "/Users/angusr/.kplr/data/lightcurves/%s" % name.zfill(9)
x, y, yerr, ivar = load_data("%s" % name, D, sc=True)
x *= 24 * 3600  # convert to seconds

fs = np.arange(500, 1800) * 1e-6
print len(fs)
n = 10000
x, y = x[:n], y[:n]
print len(x)

print "calculating..."
# starts, stops, centres = real_footprint_sc(stops)
# amp2s = spg.
fft = nufft.nufft3(x, y, fs*2*np.pi)
fft = np.imag(fft)**2 + np.real(fft)**2

print "plotting..."
plt.clf()
plt.plot(fs, fft)
print "saving..."
plt.savefig("testing")
