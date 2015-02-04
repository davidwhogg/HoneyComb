import numpy as np
import pyfits
import glob
import matplotlib.pyplot as plt

def bjd2utc(t):
    """
    # `real_footprint`

    Takes real Kepler (BJD) time values for a certain target.
    Returns the spacecraft-UTC start, stop and centre times.
    A is an array of coefficients for a sinusoid + linear trend, fit to the
    timing data of 491 asteroseismic targets that are randomly distributed on
    the CCD.
    """

    A = np.genfromtxt("A.txt").T
    w = 2*np.pi/372
    dt = 0.02043359821692

    return t + A[0]*np.sin(w*t) + A[1]*np.cos(w*t) + A[2]*t + A[3]

if __name__ == "__main__":

    # Test on a real target
    D = "/Users/angusr/angusr/data2"
    kid = "7341231"  # a kepler target with lc and sc data chosen at `random'
    fnames = []
    qs = range(17)
    x = []
    for q in qs:
        fnames.append(glob.glob("%s/Q%s_public/kplr%s-*_llc.fits"
                      % (D, q, kid.zfill(9)))[0])

    # load test fits file
    for fname in fnames:
        hdulist = pyfits.open(fname)
        tbdata = hdulist[1].data
        x.extend(tbdata["TIME"])

    # convert BJDs to UTCs
    x = np.array(x) + 2454833
    utc = bjd2utc(x)

    # plot correction
    plt.clf()
    plt.plot(x, utc, "k.")
    plt.xlabel("BJD")
    plt.ylabel("BJD-UTC (days)")
    plt.savefig("demonstrate")
