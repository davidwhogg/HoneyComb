import numpy as np
import pyfits
from astropy.time import Time

D = "/Users/angusr/angusr/data2/Q16_public"
hdulist = pyfits.open("%s/kplr007341231-2013098041711_llc.fits" % D)
hdr = hdulist[1].header

timeref = hdr["TIMEREF"]  # barycentric correction = solarsystem
tassign = hdr["TASSIGN"]  # where the time is assigned = spacecraft
timesys = hdr["TIMESYS"]  # time system = BJD
bjdrefi = hdr["BJDREFI"]  # integer part of time = 2454833
bjdreff = hdr["BJDREFF"]  # frac of day in BJD reference date
timeunit = hdr["TIMEUNIT"]  # time unit for TIME, TSTART and TSTOP
telapse = hdr["TELAPSE"]  # [d] TSTOP - TSTART
livetime = hdr["LIVETIME"]  # [d] TELAPSE multiplied by DEADC
tstart = hdr["TSTART"]  # observation start time in BJD-BJDREF
tstop = hdr["TSTOP"]  # observation stop time in BJD-BJDREF
lc_start = hdr["LC_START"]  # mid point of first cadence in MJD
lc_end = hdr["LC_END"]  # mid point of last cadence in MJD
deadc = hdr["DEADC"]  # deadtime correction
timepixr = hdr["TIMEPIXR"]  # 0.5 / bin time beginning=0 middle=0.5 end=1
tierrela = hdr["TIERRELA"]  # [d] relative time error
tierabso = hdr["TIERABSO"]  # [d] absolute time error
int_time = hdr["INT_TIME"]  # [s] photon accumulation time per frame
readtime = hdr["READTIME"]  # [s] readout time per frame
frametim = hdr["FRAMETIM"]  # [s] frame time (INT_TIME + READTIME)
num_frm = hdr["NUM_FRM"]  # number of frames per time stamp
timedel = hdr["TIMEDEL"]  # [d] time resolution of data
date_obs = hdr["DATE-OBS"]  # TSTART as UTC calendar date
date_end = hdr["DATE-END"]  # TSTOP as UTC calendar date
deadapp = hdr["DEADAPP"]  # deadtime applied
nreadout = hdr["NREADOUT"]  # number of read per cadence

iso_tstrt = " ".join((date_obs[:10], date_obs[11:23]))
iso_tstp = " ".join((date_end[:10], date_end[11:23]))
tstrt = Time(iso_tstrt, format="iso", scale="utc").jd
tstp = Time(iso_tstp, format="iso", scale="utc").jd
print tstrt - bjdrefi, "start time in UTC"
print tstart, "start time in BJD"
print "diff = ", ((tstrt-bjdrefi)-tstart)*24*3600, "s", "\n"
print tstp - bjdrefi, "stop time in UTC"
print tstop, "stop time in BJD"
print "diff = ", ((tstp-bjdrefi)-tstop)*24*3600, "s", "\n"

lc_strt = Time(lc_start, format="mjd").jd
print lc_strt - bjdrefi, "midpoint of 1st cadence on s/c"

# tbdata = hdulist[1].data
# x = tbdata["TIME"]  # BJD
# print (lc_strt-bjdrefi-x[0])*24*3600, "difference between midpoint (BJD) and x[0]"
# print (tstart-x[0])*24*3600, "difference between start point (BJD) and x[0]"
