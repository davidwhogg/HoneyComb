import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
from astropy.time import Time

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


def plot_locations(fnames):
    RA, dec = [], []
    for fname in fnames:
        hdulist = pyfits.open(fname)
        hdr = hdulist[1].header
        RA.append(hdr["RA_OBJ"])
        dec.append(hdr["DEC_OBJ"])
    plt.clf()
    plt.plot(RA, dec, '.')
    plt.savefig("CCD_positions")
    return np.array(RA), np.array(dec)

def mo_data(D, kids):

    # read timing data for all asteroseismic stars and save to text file
    qs = range(17)
    t_bjds, diffs = [], []
    for i, kid in enumerate(kids):
        print kid, i, "of", len(kids)
        fnames = []
        x = []
        for q in qs:
            try:
                fnames.append(glob.glob("%s/Q%s_public/kplr%s-*_llc.fits"
                              % (D, q, kid.zfill(9)))[0])
                t_bjd, diff, dt = get_times(fnames)
                t_bjds.extend(t_bjd)
                diffs.extend(diff)
            except:
                print "File not found: Q%s, kplr%s" % (q, kid)
    t_bjds, diffs = np.array(t_bjds), np.array(diffs)
    np.savetxt("timing_data.txt", np.transpose((t_bjds, diffs)))

    plt.clf()
    plt.plot(t_bjds, diffs, "k.")
    plt.savefig("mo_data")

if __name__ == "__main__":

    # load a bunch of Kepler stars and plot where they fall on the CCD
    D = "/Users/angusr/angusr/data2/Q16_public"
    fnames = glob.glob("%s/kplr*llc.fits" % D)
    RA, dec = plot_locations(fnames)

    # load lots of timing data for stars
    D = "/Users/angusr/angusr/data2"
    kids = np.genfromtxt("kids.txt")
    kids = [str(int(k)) for k in kids]
#     mo_data(D, kids)

    # fit to the timing data
    t_bjd, diff = np.genfromtxt("timing_data.txt").T
    ys, A = fit_sine_and_linear(t_bjd, diff, 2*np.pi/372., make_plot=True)
    np.savetxt("A.txt", np.transpose((A)))
