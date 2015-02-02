import numpy as np
import pyfits
from astropy.time import Time
import glob
import matplotlib.pyplot as plt
from sin_tests import fit_sine, show_sine

# hdr["BJDREFF"]  # frac of day in BJD reference date
# hdr["TELAPSE"]  # [d] TSTOP - TSTART
# hdr["LIVETIME"]  # [d] TELAPSE multiplied by DEADC
# hdr["DEADC"]  # deadtime correction
# hdr["TIMEPIXR"]  # 0.5 / bin time beginning=0 middle=0.5 end=1
# hdr["TIERRELA"]  # [d] relative time error
# hdr["TIERABSO"]  # [d] absolute time error
# hdr["INT_TIME"]  # [s] photon accumulation time per frame
# hdr["READTIME"]  # [s] readout time per frame
# hdr["FRAMETIM"]  # [s] frame time (INT_TIME + READTIME)
# hdr["NUM_FRM"]  # number of frames per time stamp
# hdr["TIMEDEL"]  # [d] time resolution of data
# hdr["DEADAPP"]  # deadtime applied
# hdr["NREADOUT"]  # number of read per cadence

D = "/Users/angusr/angusr/data2"
kid = "7341231"
qs = range(16)

st_utc, sp_utc, st_bjd, sp_bjd = [], [], [], []
st_diff, sp_diff = [], []
for q in qs:
    fname = glob.glob("%s/Q%s_public/kplr%s-*_llc.fits" % (D, q, kid.zfill(9)))
    print fname
    hdulist = pyfits.open(fname[0])
    hdr = hdulist[1].header

    bjdrefi = hdr["BJDREFI"]  # integer part of time = 2454833
    tstart = hdr["TSTART"]  # observation start time in BJD-BJDREF
    tstop = hdr["TSTOP"]  # observation stop time in BJD-BJDREF
    lc_start = hdr["LC_START"]  # mid point of first cadence in MJD
    lc_end = hdr["LC_END"]  # mid point of last cadence in MJD
    date_obs = hdr["DATE-OBS"]  # TSTART as UTC calendar date
    date_end = hdr["DATE-END"]  # TSTOP as UTC calendar date

    iso_tstrt = " ".join((date_obs[:10], date_obs[11:23]))
    iso_tstp = " ".join((date_end[:10], date_end[11:23]))
    tstrt = Time(iso_tstrt, format="iso", scale="utc").jd
    tstp = Time(iso_tstp, format="iso", scale="utc").jd

    print tstrt - bjdrefi, "start time in UTC"
    print tstart, "start time in BJD"
    print "diff = ", (tstrt - bjdrefi - tstart)*24*3600, "s", "\n"
    print tstp - bjdrefi, "stop time in UTC"
    print tstop, "stop time in BJD"
    print "diff = ", (tstp - bjdrefi - tstop)*24*3600, "s", "\n"

    st_utc.append(tstrt - bjdrefi)
    sp_utc.append(tstp - bjdrefi)
    st_bjd.append(tstart)
    sp_bjd.append(tstop)
    st_diff.append(tstrt - bjdrefi - tstart)
    sp_diff.append(tstp - bjdrefi - tstop)

# def fit_sine(x, y, w):
#     M = np.ones((len(x), 2+1))
#     M[:, 0] = np.sin(w*x)
#     M[:, 1] = np.cos(w*x)
#     A = np.linalg.solve(np.dot(M.T, M), np.dot(M.T, y))
#     ys = A[0]*np.sin(w*x) + A[1]*np.cos(w*x) + A[2]
#     return ys, A

def fit_sine(x, y, w):
#     M = np.vstack((np.sin(w*x), np.cos(w*x), np.ones_like(x), x))
    M = np.ones((len(x), 4))
    M[:, 0] = np.sin(w*x)
    M[:, 1] = np.cos(w*x)
    M[:, 2] = x
    A = np.linalg.solve(np.dot(M.T, M), np.dot(M.T, y))
    ys = A[0]*np.sin(w*x) + A[1]*np.cos(w*x) + A[2]*x + A[3]
    return ys, A

st_diff = np.array(st_diff)
st_bjd = np.array(st_bjd)

w = 2*np.pi*1./372
y, A = fit_sine(st_bjd, st_diff, w)
xs = np.linspace(st_bjd[0], st_bjd[-1], 1000)
ys = A[0]*np.sin(w*xs) + A[1]*np.cos(w*xs) + A[-2]

# Calculate least-square values
A = np.vstack((np.ones_like(st_bjd), st_bjd)).T
c, m = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, st_diff))

plt.clf()
plt.plot(st_bjd, st_diff, "k.")
plt.plot(st_bjd, m*st_bjd+c)
plt.plot(xs, ys+m*ys+c)
plt.plot(sp_bjd, sp_diff, "k.")
plt.xlabel("BJD")
plt.ylabel("BJD - UTC (days)")
plt.savefig("test")
