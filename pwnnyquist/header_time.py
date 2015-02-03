import numpy as np
import pyfits
from astropy.time import Time
import glob
import matplotlib.pyplot as plt

def get_times(fnames):
    """
    # `get_times`

    Take a list of fits file names (in order of ascending quarter
    and read the relevent information from the header.
    """

    t_bjd, diff = [], []
    for fname in fnames:
        hdulist = pyfits.open(fname)
        hdr = hdulist[1].header
        # read times in header
        bjdrefi = hdr["BJDREFI"]  # integer part of time = 2454833
        tstart = hdr["TSTART"]  # observation start time in BJD-BJDREF
        tstop = hdr["TSTOP"]  # observation stop time in BJD-BJDREF
        date_obs = hdr["DATE-OBS"]  # TSTART as UTC calendar date
        date_end = hdr["DATE-END"]  # TSTOP as UTC calendar date
        dt = hdr["TIMEDEL"]
        # format dates
        iso_tstrt = " ".join((date_obs[:10], date_obs[11:23]))
        iso_tstp = " ".join((date_end[:10], date_end[11:23]))
        tstrt = Time(iso_tstrt, format="iso", scale="utc").jd
        tstp = Time(iso_tstp, format="iso", scale="utc").jd
        t_bjd.append(tstart)
        t_bjd.append(tstop)
        diff.append(tstrt - bjdrefi - tstart)  # utc - bjd
        diff.append(tstp - bjdrefi - tstop)  # utc - bjd
    return np.array(t_bjd), np.array(diff), dt

def fit_sine_and_linear(x, y, w, make_plot=False):
    """
    # `fit_sine_and_linear`

    Fit a sum of sines + cosines plus a linear trend
    """

    M = np.ones((len(x), 4))
    M[:, 0] = np.sin(w*x)
    M[:, 1] = np.cos(w*x)
    M[:, 2] = x
    A = np.linalg.solve(np.dot(M.T, M), np.dot(M.T, y))
    ys = A[0]*np.sin(w*x) + A[1]*np.cos(w*x) + A[2]*x + A[3]

    # plot difference between UTC and BJD
    if make_plot == True:
        plt.clf()
        plt.plot(x, y, "k.")
        xs = np.linspace(x[0], x[-1], 1000)
        ys = A[0]*np.sin(w*xs) + A[1]*np.cos(w*xs) + A[2]*xs + A[3]
        plt.plot(xs, ys)
        plt.xlabel("BJD - 2454833")
        plt.ylabel("BJD - UTC (days)")
        plt.savefig("test")

    return ys, A

def bjd2utc(bjd, fnames):
    """
    # `bjd2utc`

    Convert BJD time to spacecraft UTC time. Fits a sinusoid with a period of
    372(?) days and linear trend to the difference between BJD and UTC time
    """

    w = 2*np.pi/372
    t_bjd, diff, dt = get_times(fnames)
    y, A = fit_sine_and_linear(t_bjd, diff, w, make_plot=True)
    return bjd + A[0]*np.sin(w*bjd) + A[1]*np.cos(w*bjd) + A[2]*bjd + A[3]

def real_footprint(t):
    """
    # `real_footprint`

    Takes real Kepler time values for a certain target.
    Returns the spacecraft-UTC start, stop and centre times.
    A is an array of coefficients for a sinusoid + linear trend, fit to the
    timing data of 491 asteroseismic targets that are randomly distributed on
    the CCD.
    """

    A = np.genfromtxt("A.txt").T
    w = 2*np.pi/372
    dt = 0.02043359821692
    stops = t + A[0]*np.sin(w*t) + A[1]*np.cos(w*t) + A[2]*t + A[3]
    starts = stops - dt
    centers = stops - .5*dt
    return starts, stops, centers

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
    starts, stops, centers = real_footprint(x)

    # plot correction
    plt.clf()
    plt.plot(x, x-stops, "k.")
    plt.xlabel("BJD")
    plt.ylabel("BJD-UTC (days)")
    plt.savefig("demonstrate")
