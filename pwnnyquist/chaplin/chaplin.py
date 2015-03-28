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
import scipy.interpolate as spi

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
def make_pgram(kid, x, y, yerr, ivar, nm, dn, fft=False):

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

def interp(grid, x, y):
    f = spi.interp1d(x, y)
    ynew = f(grid)
    return ynew

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
        fs, pgram = make_pgram(kid, x, y, yerr, ivar, nm, dn, kid)

        # run delta nu
        width = .2  # fs range is 20*dn, so this is about 4*dn
        fs, pgram = np.genfromtxt("pgram_%s.txt" % str(int(kid))).T
        mx, my, lags, pos, acor, power = deltan(fs, pgram, str(int(kid)),
                                                width, sub=100,
                                                truths=[dn, nm],
                                                smooth=False)
        np.savetxt("acor_%s.txt" % kid, np.transpose((lags, acor)))
        np.savetxt("pow_%s.txt" % kid, np.transpose((pos, power)))

        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(lags*1e6, acor, "k")
        plt.subplot(2, 1, 2)
        plt.plot(pos, power, "k")
        plt.savefig("test_%s" % kid)

    return lags, acor, pos, power

def add_up(kids, dns, nms, load=False):

    lag_list, pos_list, acor_list, power_list = [], [], [], []
    for i, kid in enumerate(kids[:3]):
        kid = str(int(kid))
        if load:
            lags, acor = np.genfromtxt("acor_%s.txt" % kid).T
            pos, power = np.genfromtxt("pow_%s.txt" % kid).T
        else:
            lags, acor, pos, power = calc_pgram_and_dnu(kid, dnu[i], nm[i])

        # save acfs and powers
        lag_list.append(lags)
        acor_list.append(acor)
        pos_list.append(pos)
        power_list.append(power)

    # find the shortest acf so that you interpolate the others onto it
#     lag_size = np.array([len(j) for j in lag_list])
#     pos_size = np.array([len(j) for j in pos_list])
#     l1 = np.where(lag_size==max(lag_size))[0]
#     l2 = np.where(pos_size==max(pos_size))[0]
#     lag_grid = lag_list[l1]
#     pos_grid = pos_list[l2]
    dl = lag_list[0][1] - lag_list[0][0]
    print dl
    max_start = max([j[0] for j in lag_list])
    min_end = min([j[0] for j in lag_list])
    lag_grid = np.arange(max_start, min_end, dl)
    print len(lag_grid), len(lag_list[0])
    assert 0

    # interpolate so that you can add up the signals
    lag_interps = np.zeros((len(lag_list), len(lag_grid)))
    pos_interps = np.zeros((len(pos_list), len(pos_grid)))
    for i, lags in enumerate(lag_list):
        l = lags > lag_grid[0]
        print min(lag_grid), min(lags[l])
        print max(lag_grid), max(lags[l])
        new_acor = interp(lag_grid, lags[l], acor_list[i][l])
        l = pos > pos_grid[0]
        print min(pos_grid), min(pos[l])
        new_power = interp(pos_grid, pos_list[i][l], power_list[i][l])
        lag_interps[i, :] = new_acor
        pos_interps[i, :] = new_power

    # make plots of summed acfs and power thingys
    plt.clf()
    plt.subplot(2, 1, 1)
    plt.plot(lag_grid, np.sum(lag_interps, axis=0))
    plt.subplot(2, 1, 2)
    plt.plot(pos_grid, np.sum(pos_interps, axis=0))
    plt.savefig("test")

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

    add_up(kids[l], dnu[l], nm[l])
