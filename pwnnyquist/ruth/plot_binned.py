import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
from params import plot_params
reb = plot_params()

fnames = glob.glob("spg??.h5")
plt.clf()
for fname in fnames:
    print fname
    with h5py.File(fname, "r") as f:
        fs = f["spg"][0, :]
        amp2s = f["spg"][1, :]
    nbins = int(fname[3:5])
#     nbins = int(fname[3:4])
    print nbins
    plt.plot(fs*1e6-1160, amp2s, label="$%s~\mathrm{points~per~bin}$" % nbins)
plt.legend()
plt.xlim(min(fs*1e6-1160), max(fs*1e6-1160))
plt.ylabel("$\mathrm{Amp}^2$")
plt.xlabel("$\\nu-1160~\mathrm{(}\mu\mathrm{Hz)}$")
plt.savefig("binning_comparison_spg")

fnames = glob.glob("fft??.h5")
plt.clf()
for fname in fnames:
    print fname
    with h5py.File(fname, "r") as f:
        fs = f["fft"][0, :]
        amp2s = f["fft"][1, :]
    nbins = int(fname[3:5])
#     nbins = int(fname[3:4])
    plt.plot(fs*1e6-1160, amp2s, label="$%s~\mathrm{points~per~bin}$" % nbins)
plt.legend()
plt.xlim(min(fs*1e6-1160), max(fs*1e6-1160))
plt.ylabel("$\mathrm{Amp}^2$")
plt.xlabel("$\\nu-1160~\mathrm{(}\mu\mathrm{Hz)}$")
plt.savefig("binning_comparison_fft")

fnames = glob.glob("ls??.h5")
plt.clf()
for fname in fnames:
    print fname
    with h5py.File(fname, "r") as f:
        fs = f["ls"][0, :]
        amp2s = f["ls"][1, :]
    nbins = int(fname[2:4])
#     nbins = int(fname[2:3])
    plt.plot(fs*1e6-1160, amp2s, label="$%s~\mathrm{points~per~bin}$" % nbins)
plt.legend()
plt.xlim(min(fs*1e6-1160), max(fs*1e6-1160))
plt.ylabel("$\mathrm{Amp}^2$")
plt.xlabel("$\\nu-1160~\mathrm{(}\mu\mathrm{Hz)}$")
plt.savefig("binning_comparison_ls")
