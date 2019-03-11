#
# ladeeocctime_all.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: September 2017
# description: analyzing processed LADEE UVS data
#

import os
import argparse
import pickle, pprint
from numpy import sqrt, pi, exp, linspace, random
import numpy as np
import scipy as sp
import pydl as dl
from scipy.optimize import curve_fit
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import spiceypy as spice
from scipy import fftpack
import utilities as ut


# converts calendar time in Julian or UTC to ephemeris seconds past J2000
def cspice_utc2et(kernelfile, stimes, ttimes):
    global ets, etx
    nt = len(stimes)
    ets = np.zeros(nt)
    etx = np.zeros(nt)
    spice_ver = spice.tkvrsn('TOOLKIT')
    spice.furnsh(kernelfile)
    kernels_loaded = spice.ktotal("ALL")
    #print(kernels_loaded)
    
    for i in range(0,nt):
        sstr = stimes[i][:23]
        tstr = ttimes[i][:23]
        et1 = spice.utc2et(sstr)
        et2 = spice.utc2et(tstr)
        ets[i] = 0.5*(et1 + et2)

    for i in range(0,nt):
        etx[i] = ets[i] - ets[0]

    print("Spice version: ", spice_ver)


# fit a Guassian to input data
def fit_data():
    nw = len(wavels)
    nt = len(ets)
    global profs, derivs, fits, profs_base
    profs = np.zeros(nw*nt).reshape(nw,nt)
    derivs = np.zeros(nw*nt).reshape(nw, nt)
    fits = np.zeros(3*nw).reshape(nw,3)
    profs_base = np.zeros(nt)
    
    #fit to gaussian
    for i in range(nw):
        profs[i] = specs[i]
        derivs[i] = np.gradient(profs[i])
        ymax = np.argmax(-derivs[i])
        guess = [1, etx[ymax], 1]
        bounds = ([0, etx[ymax]-dchan, 0], [1., etx[ymax]+dchan, 1.])
        fits[i] = np.fit_gaussian(etx, -derivs[i], guess, bounds)

    # use peak location from fit to offset etx
    goodfit = []
    for i in range(nw):
        if (fits[i][1] != 0):
            goodfit.append(fits[i][1])
    
    global ety
    ety = np.zeros(nw*nt).reshape(nw,nt)
    mt = np.median(goodfit)
    for i in range(nw):
        for j in range(nt):
            ety[i][j] = etx[j] - mt


# normalize the data
def norm_data():
    nw = len(wavels)
    nt = len(ets)
  
    for i in range(nw):
        qb = []
        for j in range(nt):
            if (ety[i][j] > br[0] and ety[i][j] < br[1]):
                qb.append(j)
        etq = np.asarray(qb)
        lq = len(etq)
        pmean = np.mean(profs[i][etq[0]:etq[lq-1]])
        profs[i] = profs[i] / pmean


# plot results
def plot_results():
    
    nw = len(wavels)
    nt = len(ety[0])

    xr = [nt-180, nt-60]
    yr1 = [229, 822]
    
    f, ax1 = plt.subplots(1)
    f.suptitle(filenameonly, fontsize=12)
    
    wavelengths = np.linspace(yr1[0], yr1[1], 20)
    times = np.linspace(ety[0][xr[0]], ety[0][xr[1]], 20)

    # wavelength of around 400 nm
    widx = 300
    profs_avg = profs[widx, xr[0]:xr[1]]

    nuprofs = profs[:, xr[0]:xr[1]] - profs_avg
    ax1.imshow(nuprofs, shape=(nw,nt), cmap='Greys', extent=[ety[0][xr[0]], ety[0][xr[1]], yr1[0], yr1[1]], origin='lower', interpolation='bicubic', aspect='auto')
    ax1.set_title("Trans. Difference")
    ax1.set_xlabel('Time (seconds from transit)', fontsize=8)
    ax1.set_ylabel('Wavelengths (nm)', fontsize=8)
    ax1.tick_params(labelsize=6)
    ax1.set_xticks(times)
    ax1.set_yticks(wavelengths)
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(.5))
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(100))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(25))

    if (saveplot):
        plotfile = os.path.splitext(filenameonly)[0] + "_all.pdf"
        plt.savefig(plotfile, figsize=(9,9), format='pdf')
    else:
        plt.show()


# load processed uvs data
def readuvs_data(filename):
    global specs, wavels, stimes, ttimes, tellat, tellon, sunlat, sunlon, telalt, solgraz
    
    if (os.path.exists(filename) == False):
        print("File does not exist!")
        
    print("File " + filename + " exists!")
    dataload = pickle.load(open(filename, "rb"))
    
    oldspec = dataload['specs']
    specs = oldspec.T
    wavels = dataload['wavels']
    ttimes = dataload['ttimes']
    stimes = dataload['stimes']
    tellat = dataload['tellat']
    tellon = dataload['tellon']
    sunlat = dataload['sunlat']
    sunlon = dataload['sunlon']
    telalt = dataload['telalt']
    solgraz = dataload['solgraze']


# some of the wavelengths have spikes in the data.
# check for this and ignore that data
def clean_data():
    ns = len(specs)
    # only applying check for spikes for first 300 wavelengths
    ns = 300#len(wavels)
    nt = len(etx)
    diff = np.zeros(ns*nt).reshape(ns, nt)
    smooth = np.zeros(ns*nt).reshape(ns, nt)
    
    for i in range(ns):
        diff[i] = specs[i] - dl.smooth(specs[i], 3)

    threshold = 0.03
    for i in range(ns):
        spectmp = []
        ttmp = []
        for j in range(nt):
            if (diff[i][j] < threshold):
                spectmp.append(specs[i][j])
                ttmp.append(j)
        trange = np.arange(0,len(etx))
        specs[i] = np.interp(trange, np.asarray(ttmp), np.asarray(spectmp))


def main():
# parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-f", "--filename", help="data file")
    parser.add_argument("-k", "--kernels", help="kernel meta file")
    parser.add_argument("-s", "--saveplot", help="save plot file", action="store_true")
    args = parser.parse_args()
    
# define some global variables
    global bc, bw, dchan, br

    global saveplot
    saveplot = args.saveplot
    
    bc = [325, 375, 425, 475, 525, 575, 625, 675, 725, 775]
    bw = 25
    dchan = 50
    br = [-20, -10]
    
# get file name (no directory)
    filename = args.filename
    global filenameonly;
    filenameonly = os.path.basename(filename)
    #filenameonly = os.path.splitext(filenameonly)[0]

# read processed uvs data
    readuvs_data(filename)
    
# compute ephemeris times using spice
    kernelfile = args.kernels
    cspice_utc2et(kernelfile, stimes, ttimes)
    
# clean data
    clean_data()

# fit data (determine fit peak using derivative)
    fit_data()
    
# normalize data
    norm_data()

# plot results
    plot_results()


if __name__ == "__main__":
    main()


