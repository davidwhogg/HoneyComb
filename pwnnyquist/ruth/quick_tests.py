import numpy as np
import glob
import pyfits

D = "/Users/angusr/angusr/data2"
kid = "3427720"  # a sc target from Metcalfe et al.

fnames = []
qs = range(17)
for q in qs:
    fnames.append(glob.glob("%s/Q%s_public/kplr%s-*_llc.fits"
                  % (D, q, kid.zfill(9)))[0])

for fname in fnames:
    hdulist = pyfits.open(fname)
    hdr = hdulist[1].header
