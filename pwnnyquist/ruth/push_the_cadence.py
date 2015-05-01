import numpy as np
import matplotlib.pyplot as plt
import glob
from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel
import superpgram as spg
from header_time import real_footprint, real_footprint_sc
from fft import load_data
import h5py
from find_delta_nu import find_delta_nu as deltan
from astero import astero
import nufft

def arbitrary_footprint(t, dt):
    """
    # `arbitrary_footprint`

    Takes real Kepler time values (in seconds) for a certain target and returns
    estimates of the starts, stops and centres, given an arbitrary dt
    takes dt and t in seconds
    """
    stops = t
    starts = t - dt
    centres = t - .5*dt
    return starts, stops, centres

# calc spg for a range of cadences.
# option to generate and save new results or load precomputed spg
def cadence_tests(name, fs, load=False, test=True):

    print "loading data..."
    D = "/Users/angusr/.kplr/data/lightcurves/%s" % name.zfill(9)
    x, y, yerr, ivar = load_data("%s" % name, D, sc=True)
    x *= 24 * 3600  # convert to seconds
    if test:
        n = 1000
        x, y, yerr, ivar = x[:n], y[:n], yerr[:n], ivar[:n]

    print len(fs), "frequencies"
    ws = np.pi*2*fs

    dt = 6.81119940564e-4 * 24 * 3600 # interval between observations (seconds)
    cadences = np.arange(6)  # integer multiples of sc sampling freq

    plt.clf()
    for i in cadences:
        starts, stops, centres = arbitrary_footprint(x, dt*(i+1))
        if load:
            print "loading spg"
            with h5py.File("spg%s.h5" % cadences[i], "r") as f:
                fs = f["spg"][0, :]
                amp2s = f["spg"][1, :]
        else:
            print "computing superpgram..."
#             amp2s = spg.superpgram(starts, stops, y, ivar, ws)
            amp2s = nufft.nufft3(stops, y, ws)
            print "saving spg"
#             f = h5py.File("spg%s.h5" % cadences[i], "w")
#             data = f.create_dataset("spg", (2, len(fs)))
#             data[0, :] = fs
#             data[1, :] = amp2s
#             f.close()
        plt.plot(fs, amp2s, linewidth=cadences[i]+1)
        plt.savefig("fft")
        assert 0
    print "saving plot"
    plt.savefig("cadence_tests")

if __name__ == "__main__":
    name = "6442183"  # lower nu max poster child
#     fs = np.arange(1000, 1400, 1e-1) * 1e-6
    fs = np.arange(1159, 1161, 1e-2) * 1e-6
    cadence_tests(name, fs)
