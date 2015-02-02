import numpy as np
import pyfits
from astropy.time import Time
import glob
import matplotlib.pyplot as plt
from sin_tests import fit_sine, show_sine

def get_times(D, kid):
    qs = range(17)
    t_bjd, diff = [], []
    for q in qs:
        # load fits file
        fname = glob.glob("%s/Q%s_public/kplr%s-*_llc.fits" % (D, q, kid.zfill(9)))
        hdulist = pyfits.open(fname[0])
        hdr = hdulist[1].header

        # read times in header
        bjdrefi = hdr["BJDREFI"]  # integer part of time = 2454833
        tstart = hdr["TSTART"]  # observation start time in BJD-BJDREF
        tstop = hdr["TSTOP"]  # observation stop time in BJD-BJDREF
        date_obs = hdr["DATE-OBS"]  # TSTART as UTC calendar date
        date_end = hdr["DATE-END"]  # TSTOP as UTC calendar date

        # format dates
        iso_tstrt = " ".join((date_obs[:10], date_obs[11:23]))
        iso_tstp = " ".join((date_end[:10], date_end[11:23]))
        tstrt = Time(iso_tstrt, format="iso", scale="utc").jd
        tstp = Time(iso_tstp, format="iso", scale="utc").jd

        t_bjd.append(tstart)
        t_bjd.append(tstop)
        diff.append(tstrt - bjdrefi - tstart)  # utc - bjd
        diff.append(tstp - bjdrefi - tstop)  # utc - bjd
    return np.array(t_bjd), np.array(diff)

def fit_sine_and_linear(x, y, w):
    M = np.ones((len(x), 4))
    M[:, 0] = np.sin(w*x)
    M[:, 1] = np.cos(w*x)
    M[:, 2] = x
    A = np.linalg.solve(np.dot(M.T, M), np.dot(M.T, y))
    ys = A[0]*np.sin(w*x) + A[1]*np.cos(w*x) + A[2]*x + A[3]
    return ys, A

# this is the main function that takes BJDs and converts to UTCs
def bjd2utc(bjd):
    A = np.genfromtxt("t_coeffs.txt")
    w = 2*np.pi/372
    return bjd + A[0]*np.sin(w*bjd) + A[1]*np.cos(w*bjd) + A[2]*bjd + A[3]

if __name__ == "__main__":

    D = "/Users/angusr/angusr/data2"
    kid = "7341231"

    w = 2*np.pi*1./372
    t_bjd, diff = get_times(D, kid)  # find all recorded times in headers
    y, A = fit_sine_and_linear(t_bjd, diff, w)  # fit sin+cos + linear trend
    np.savetxt("t_coeffs.txt", np.transpose((A)))

    # plot difference between UTC and BJD
    plt.clf()
    plt.plot(t_bjd, diff, "k.")
    xs = np.linspace(t_bjd[0], t_bjd[-1], 1000)
    ys = A[0]*np.sin(w*xs) + A[1]*np.cos(w*xs) + A[2]*xs + A[3]
    plt.plot(xs, ys)
    plt.xlabel("BJD - 2454833")
    plt.ylabel("BJD - UTC (days)")
    plt.savefig("test")

    # load test fits file
    hdulist = pyfits.open("%s/Q16_public/kplr%s-2013098041711_llc.fits"
                          % (D, kid.zfill(9)))
    tbdata = hdulist[1].data
    x = tbdata["TIME"]

    # convert BJDs to UTCs
    utc = bjd2utc(x)

    # plot correction
    plt.clf()
    plt.plot(x, x-utc, "k.")
    plt.xlabel("BJD")
    plt.ylabel("BJD-UTC (days)")
    plt.savefig("demonstrate_q16")
