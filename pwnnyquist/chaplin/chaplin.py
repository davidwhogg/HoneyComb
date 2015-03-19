# script to measure nu_max and delta_nu for all chaplin stars

import numpy as np
import matplotlib.pyplot as plt
import glob
from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel
import kplr
from useful_things import color2teff
from scaling_relations import delta_nu, nu_max
import superpgram as spg
from header_time import real_footprint, real_footprint_sc
from fft import load_data
import h5py
import nufft
from find_delta_nu import find_delta_nu as deltan
from astero import astero
import os

# using the kic colours, use isochrones to produce posterior samples
# of the physical stellar parameters
def colour2physical(kid, plot=False):

    # fetch photometry
    client = kplr.API()
    star = client.star(kid)
    mags = {'H': (star.kic_hmag, 0.08),
            'J': (star.kic_jmag, 0.08),
            'K': (star.kic_kmag, 0.08),
            'g':(star.kic_gmag, 0.08),
            'r':(star.kic_rmag, 0.08)}

    # find the mass and radius using isochrones
    dar = Dartmouth_Isochrone()
    smod_phot = StarModel(dar, **mags)
    result = smod_phot.fit_mcmc()

    # save samples and triangle plot
    if plot:
        smod_phot.triangle_plots("test", format="pdf")
    smod_phot.save_hdf("samples_%s.h5" % str(int(kid)))
    max_like = smod_phot.maxlike()
    print max_like
    np.savetxt("maxlike.txt", max_like)

# load samples and find astero parameters
def samples2astero(fname, max_like=False):
    print fname

    if max_like:
        max_like = np.genfromtxt(fname).T
    else:
        with h5py.File(fname) as f:
            dataset = f["/samples"][...]

    # find the nu_max and delta_nu
#     nm = nu_max(mass, radius, teff)
#     dn = delta_nu(mass, radius)

# compute a periodogram around that area using short cadence
def make_pgram(x, y, yerr, ivar, nm, dn, fft=False):

    # compute spg or fft
    fs = np.arange(nm-(10*dn), nm+(10*dn), 1e-2) * 1e-6
    if fft:
        comp_pgram = nufft.nufft3(x, y, fs*2*np.pi)
        pgram = np.imag(comp_pgram)**2 + np.real(comp_pgram)**2
    else:
        starts, stops, centres = real_footprint_sc(x)
        pgram = spg.superpgram(starts, stops, y, ivar, fs*2*np.pi)

    # save spg
    np.savetxt("pgram_%s.txt" % str(int(kid)), np.transpose((fs, pgram)))
    plt.clf()
    plt.plot(fs, pgram)
    plt.savefig("pgram_%s" % str(int(kid)))
    return fs, pgram

# add up the signal for stars of the same type... how? try adding both the
# matrix and the collapsed
def add_up(kids):
    fname = "pow_%s.txt" % str(int(kids[0]))
    if os.path.isfile(fname):
        pos, power = np.genfromtxt("pow_%s.txt" % str(int(kids[0]))).T
        for i, kid in enumerate(kids[1:]):
            kid = str(int(kid))
            data = np.genfromtxt("pow_%s.txt" % kid).T
            power += data[-1]
    plt.clf()
    plt.plot(pos, power, "k")
    plt.savefig("test")

# calculates the pgram and dnu, etc
def calc_pgram_and_dnu(kid, dn, nm):
        D = "/Users/angusr/.kplr/data/lightcurves/%s" % kid.zfill(9)
        x, y, yerr, ivar = load_data(kid, D, sc=True)
        if len(x):  # only run if you find data
            print "data found for ", kid
            print "dnu = ", dn, "nm = ", nm
            fs, pgram = make_pgram(x, y, yerr, ivar, nm, dn, kid)

            # run delta nu
            width = .1
            fs, pgram = np.genfromtxt("pgram_%s.txt" % str(int(kid))).T
            mx, my, lags, pos, acor, power = deltan(fs, pgram, str(int(kid)),
                                                    width, sub=100,
                                                    truths=[dn, nm],
                                                    smooth=False)
            np.savetxt("acor_%s.txt" % kid, np.transpose((lags, acor)))
            np.savetxt("pow_%s.txt" % kid, np.transpose((pos, power)))

            plt.clf()
            plt.subplot(2, 1, 1)
            plt.plot(lags, acor, "k")
            plt.axvline(nm, color="r", linestyle="--")
            plt.subplot(2, 1, 2)
            plt.plot(pos, power, "k")
            plt.savefig("test_%s" % kid)

#     add_up(kids[l])

if __name__ == "__main__":

#     kid = "6442183"
#     dn = 65.
#     nm = 1160.
#     calc_pgram_and_dnu(kid, dn, nm)

    data = astero()
    kids = data.iKID
    dnu, nm = data.idnu, data.inu_max

    # select stars with similar nu_max
    l = (500 < nm) * (nm < 600)
    print len(kids[l])

#     fname = "samples_%s.h5" % kid
#     fname = "maxlike.txt"

#     colour2physical(kid)
#     samples2astero(fname)

    for i, kid in enumerate(kids[l][:3]):
        kid = str(int(kid))
        calc_pgram_and_dnu(kid, dnu[l][i], nm[l][i])
