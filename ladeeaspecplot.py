#
# ladeeaspecplot.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: December 2017
# description: analyzing processed LADEE UVS data
#

import os
import argparse
import pickle, pprint
import pydl as dl
from numpy import sqrt, pi, exp, linspace, random
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy import ndimage
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter, FormatStrFormatter)
import spiceypy as spice
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

    #print(spice_ver)


# load processed uvs data
def readuvs_data(filename):
    global specs, wavels, stimes, ttimes, tellat, tellon, sunlat, sunlon, telalt, solgraz, telaltx
    
    if (os.path.exists(filename) == False):
        print("File " + filename + " does not exist!")
        
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
    telaltx = dataload['telaltx']
    solgraz = dataload['solgraze']


# solar spectrum data
def read_solarflux(filename):
    global solflux

    if (os.path.exists(filename) == False):
        print("File " + filename + " does not exist!")

    filein = open(filename, "r")

    solflux = []
    for line in filein:
        newline = line.strip()
        rowdata = newline.split()
        rowdata[0] = float(rowdata[0]) * 1.0e03
        rowdata[1] = float(rowdata[1])
        solflux.append(rowdata)
    #print(solflux)

    filein.close()


# some of the wavelengths have spikes in the data.
# check for this and ignore that data
def clean_data():
    ns = len(specs)
    # only applying check for spikes for first 300 wavelengths
    ns = 300
    nt = len(etx)
    diff = np.zeros(ns*nt).reshape(ns,nt)
    smooth = np.zeros(ns*nt).reshape(ns,nt)
    
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


# plot results
def plot_results():
    figdim = 9
    
    xr = [-10, 10]

    nw = len(wavels)
    widx = 0
    for i in range(nw):
        if (abs(float(wavels[i][1]) - 400.) < 1.0 and float(wavels[i][1]) > 400.):
            widx = i
            break

    yrange = np.linspace(0., 25., 6., endpoint=True)
    
    f = plt.figure(figsize=(figdim,figdim))
    ax3 = plt.subplot2grid((3, 3), (2, 0), colspan=1)
    ax2 = plt.subplot2grid((3, 3), (1, 0), colspan=1, sharex=ax3)
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=1, sharex=ax2)
    ax4 = plt.subplot2grid((3, 3), (0, 1), colspan=2, rowspan=1)
    ax5 = plt.subplot2grid((3, 3), (2, 1), colspan=2, rowspan=1, sharex=ax4)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax5.get_xticklabels(), visible=False)

    ax1.xaxis.set_tick_params(which='both', direction='inout')
    ax2.xaxis.set_tick_params(which='both', direction='inout')
    ax3.xaxis.set_tick_params(which='both', direction='in')
    ax1.yaxis.set_tick_params(which='both', direction='in')
    ax2.yaxis.set_tick_params(which='both', direction='in')
    ax3.yaxis.set_tick_params(which='both', direction='in')
    ax4.xaxis.set_tick_params(which='both', direction='in')
    ax4.yaxis.set_tick_params(which='both', direction='in')
    #ax4.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax5.xaxis.set_tick_params(which='both', direction='inout')
    ax5.yaxis.set_tick_params(which='both', direction='in')
    #ax5.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax3.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax4.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax5.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax3.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax3.yaxis.set_minor_locator(ticker.MultipleLocator(1e-04))
    #ax4.yaxis.set_minor_locator(ticker.MultipleLocator(.1e-04))
    ax5.yaxis.set_minor_locator(ticker.MultipleLocator(.1e-04))

    f.subplots_adjust(hspace=0.00)

    f.suptitle(filenameonly, fontsize=12)

    telaltr = np.asarray(telaltx).astype(np.float)
    solgrazr = np.asarray(solgraz).astype(np.float)

    tmax = np.amax(telaltr)
    smax = np.amax(solgrazr)
    tmin = np.amin(telaltr)
    smin = np.amin(solgrazr)
    yr1 = [0, 1.1*tmax]
    yr2 = [0, 1.1*smax]
    #yr3 = [0, 1.12*np.amax(specs[widx])]
    yr3 = [1.0*np.amin(specs[widx]), 0.0025*np.amax(specs[widx])]
    yr4 = [-2.0e-04, 2.0e-04]
    xr1 = [0, 500]
    xr4 = [200, 800]
    ax1.set_ylim(yr1)
    ax2.set_ylim(yr2)
    ax3.set_ylim(yr3)
    #ax4.set_ylim(yr4)
    #ax5.set_ylim(yr4)
    ax1.set_xlim(xr1)
    ax2.set_xlim(xr1)
    ax3.set_xlim(xr1)
    ax4.set_xlim(xr4)
    ax5.set_xlim(xr4)

    lines1 = ax1.plot(etx, telaltr, lw=1.0)

    # define regions to use for computing specs excess
    chan = len(etx)
    tx1 = etx[np.amin(np.where(specs[chan,:] < 0.01))]
    #print(tx1)
    ty1 = np.amax(etx, 0)
    tx = [tx1+xs_bands[0]*(ty1-tx1), tx1+xs_bands[1]*(ty1-tx1), tx1+xs_bands[2]*(ty1-tx1),tx1+xs_bands[3]*(ty1-tx1)]

    for i in range(0,2):
        ax1.axvline(tx[i], color='red', alpha=0.5, lw=1.0, ls='--')
        ax2.axvline(tx[i], color='red', alpha=0.5, lw=1.0, ls='--')
        ax3.axvline(tx[i], color='red', alpha=0.5, lw=1.0, ls='--')
    for i in range(2,4):
        ax1.axvline(tx[i], color='blue', alpha=0.5, lw=1.0, ls='--')
        ax2.axvline(tx[i], color='blue', alpha=0.5, lw=1.0, ls='--')
        ax3.axvline(tx[i], color='blue', alpha=0.5, lw=1.0, ls='--')

    ax1.set_ylabel('Telescope Boresight \nAltitude (km)', fontsize=10)

    #ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ymax = np.around(smax) + 5 - np.around(smax)%5
    nticks = np.around(smax)/5 + 2
    yrange = np.linspace(0., ymax, nticks, endpoint=True)
    ax1.set_yticks(yrange)
    #ax1.set_xticks(np.linspace(xr[0], xr[1], int(0.5*(xr[1]-xr[0])+1), endpoint=True))
    #ax1.set_xticks(np.linspace(xr[0], xr[1], int(2*(xr[1]-xr[0])+1), endpoint=True), minor=True)

    # shift first tick up so it doesn't overlap tick on other plot
    for tick in ax1.yaxis.get_majorticklabels():
        tick.set_verticalalignment("bottom")
        break

    lines2 = ax2.plot(etx, solgrazr, lw=1.0)

    ax2.set_ylabel('Solar Graze \nAltitude (km)', fontsize=10)
    #ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax2.set_yticks(yrange)

    # shift first tick up so it doesn't overlap tick on other plot
    for tick in ax2.yaxis.get_majorticklabels():
        tick.set_verticalalignment("bottom")
        break

    lines3 = ax3.plot(etx, specs[widx], lw=1.0)

    ax3.set_ylabel('Radiance at 400 nm \n(W/m$^2$/nm/sr)', fontsize=10)
    ax3.yaxis.set_label_coords(-0.27, 0.5)
    #ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

    ax3.set_xlabel('Time from Start (s)')

    # shift first tick up so it doesn't overlap tick on other plot
    for tick in ax3.yaxis.get_majorticklabels():
        tick.set_verticalalignment("bottom")
        break

    # compute specs excess
    rb = 256
    q1 = []
    q2 = []
    for i in range(0, len(etx)):
        if (etx[i] > tx[0] and etx[i] < tx[1]):
            q1.append(i)
        if (etx[i] > tx[2] and etx[i] < tx[3]):
            q2.append(i)
    spec1 = dl.rebin(specs[:,q1], (rb,1)).reshape(rb)
    spec2 = dl.rebin(specs[:,q2], (rb,1)).reshape(rb)
    wl = dl.rebin(wavels[:,1], (rb,))
    solf = np.interp(wl, np.asarray(solflux)[:,0], np.asarray(solflux)[:,1])
    #nspec = (spec2 - spec1) * np.pi / (10.*solf)
    nspec = (spec1 - spec2) * np.pi / (10.*solf)

    ax4.set_xlabel('wavelength (nm)', fontsize=10)
    ax4.set_ylabel('I/F', fontsize=10)
    ax4.yaxis.set_label_coords(-0.10, 0.5)
    #ax4.plot(wl, nspec, lw=1.0)
    ax4.semilogy(wl, nspec, lw=1.0)
    #ax4.yaxis.set_major_formatter(ScalarFormatter())
    #ax4.yaxis.set_minor_formatter(NullFormatter())
    #ax4.set_yticks(np.linspace(0.1, 1, 10))
    #ax4.set_ylim([1.0e-05, 0.2])
    [left, bottom, width, height] = [0.47, 0.11, 0.50, 0.5*0.77]
    ax4.set_position([left, bottom, width, height], which='both')

    ax5.set_ylabel('Radiance (W/m$^2$/nm/sr)', fontsize=10)
    ax5.yaxis.set_label_coords(-0.12, 0.5)
    ax5.plot(wl, spec1-spec2, lw=1.0)
    yr5 = [-1.0e-04, 4.0e-04]
    ax5.set_ylim(yr5)
    [left, bottom, width, height] = [0.47, 0.49, 0.50, 0.5*0.77]
    ax5.set_position([left, bottom, width, height], which='both')

    # shift last tick down so it doesn't overlap tick on other plot
    ticks = ax4.yaxis.get_major_ticks()
    nticks = len(ticks)
    ticknum = 0
    for tick in ax4.yaxis.get_majorticklabels():
        if (ticknum == nticks-1):
            tick.set_verticalalignment("top")
        ticknum += 1

    # shift first tick up so it doesn't overlap tick on other plot
    for tick in ax5.yaxis.get_majorticklabels():
        tick.set_verticalalignment("bottom")
        break

    if (saveplot):
        plotfile = os.path.splitext(filenameonly)[0] + "_spec.pdf"
        plt.savefig(plotfile, figsize=(figdim,figdim), format='pdf')
    else:
        plt.show()


def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-f", "--filename", help="file name")
    parser.add_argument("-k", "--kernels", help="kernel meta file")
    parser.add_argument("-s", "--solarflux", help="solar flux file")
    parser.add_argument("-p", "--saveplot", help="save plot file", action="store_true")
    parser.add_argument("-x", "--xs_bands", nargs='+', type=float, help="regions for computing specs excess")
    args = parser.parse_args()

    global saveplot, xs_bands

    # get args
    filename = args.filename
    kernelfile = args.kernels
    solarflux = args.solarflux
    saveplot = args.saveplot

    xs_bands = [0.1, 0.4, 0.6, 0.9]
    if (args.xs_bands != None):
        xs_bands = args.xs_bands

# get file name (no directory)
    global filenameonly;
    filenameonly = os.path.basename(filename)
    #filenameonly = os.path.splitext(filenameonly)[0]

    # read processed uvs data
    readuvs_data(filename)
    
    # compute ephemeris times using spice
    kernelfile = args.kernels
    cspice_utc2et(kernelfile, stimes, ttimes)

    # clean data
    #clean_data()

    # read solarflux data
    read_solarflux(solarflux)

    # plot results
    plot_results()


if __name__ == "__main__":
    main()
