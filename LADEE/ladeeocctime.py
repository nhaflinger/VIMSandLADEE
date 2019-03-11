#
# ladeeocctime.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: August 2017
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

    print(spice_ver)


# fit a Guassian to input data
def fit_data():
    nc = len(bc)
    nw = len(wavels)
    nt = len(ets)
    global profs, derivs, fits
    profs = np.zeros(nc*nt).reshape(nc,nt)
    derivs = np.zeros(nc*nt).reshape(nc, nt)
    fits = np.zeros(3*nc).reshape(nc,3)
    
    #find indices of values in wavels that lie in ranges of bc
    for i in range(nc):
        qc = []
        specl = []
        for j in range(nw):
            if (wavels[j][1] > (bc[i]-bw) and wavels[j][1] < (bc[i]+bw)):
                qc.append(j)
        for j in range(len(qc)):
            for k in range(nt):
                specl.append(specs[qc[j]][k])

        lqc = len(qc)
        nuspec = np.asarray(specl).reshape(lqc,nt)
        
        if (qc[0] != -1):
            #profs[i] = rebin(nuspec[i], nt, 1)
            profs[i] = np.mean(nuspec, axis=0)
            derivs[i] = np.gradient(profs[i])
            ymax = np.argmax(-derivs[i])
            guess = [1, etx[ymax], 1]
            bounds = ([0, etx[ymax]-dchan, 0], [1., etx[ymax]+dchan, 1.])
            fits[i] = ut.fit_gaussian(etx, -derivs[i], guess, bounds)

    #use peak location from fit to offset etx
    goodfit = []
    for i in range(nc):
        if (fits[i][1] != 0):
            goodfit.append(fits[i][1])
    
    global ety
    ety = np.zeros(nc*nt).reshape(nc,nt)
    mt = np.median(goodfit)
    for i in range(nc):
        for j in range(nt):
            ety[i][j] = etx[j] - mt


# normalize the data
def norm_data():
    nc = len(bc)
    nt = len(ets)
    for i in range(nc):
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

    colors = ['grey','magenta','violet','blue','cyan','green','orange','red','maroon','brown']
    #colors=['violet','blue','cyan','orange','red','maroon']
    #colors=['grey','green','forest green', 'yellow','red','magenta','blue','cyan','turquoise']

    figdim = 9
    
    xr = [-10, 10]
    yr1 = [0, 1.2]
    yr2 = [-0.05, 0.05]
    
    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
    f.suptitle(filenameonly, fontsize=12)

    nc = len(bc)
    for i in range(nc):
        lines = ax1.plot(ety[i], profs[i])
        ax1.set_xlim(xr)
        ax1.set_ylim(yr1)
        ax1.set_ylabel('Transmission')
        ax1.set_xticks(np.linspace(xr[0], xr[1], int(0.5*(xr[1]-xr[0])+1), endpoint=True))
        ax1.set_xticks(np.linspace(xr[0], xr[1], int(2*(xr[1]-xr[0])+1), endpoint=True), minor=True)
        ax1.text(-xr[1]*0.25, yr1[1]*1.1,  "Wavelength Range (nm)", fontsize=8, fontstyle='normal')
        for j in range(len(bc)):
            wstr = str(bc[j]-bw) + "-" + str(bc[j]+bw)
            ax1.text(-0.8*figdim+j*1.5, yr1[1]*1.05, wstr, fontsize=6, fontstyle='normal', color=colors[j])
        plt.setp(lines, color=colors[i], linewidth=1)
    
    for i in range(nc):
        normed = ax2.plot(ety[i], profs[i] - profs[2])
        ax2.set_xlim(xr)
        ax2.set_ylim(yr2)
        ax2.set_xlabel('Time (seconds from transit)')
        ax2.set_ylabel('Trans. Difference')
        ax2.set_xticks(np.linspace(xr[0], xr[1], int(0.5*(xr[1]-xr[0])+1), endpoint=True))
        ax2.set_xticks(np.linspace(xr[0], xr[1], int(2*(xr[1]-xr[0])+1), endpoint=True), minor=True)
        plt.setp(normed, color=colors[i], linewidth=1)

    f.subplots_adjust(hspace=0.1)

    if (saveplot):
        plotfile = os.path.splitext(filenameonly)[0] + "_trans.pdf"
        plt.savefig(plotfile, figsize=(figdim,figdim), format='pdf')
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


# load elevation data
def load_elevation(elevationfile):
    elevation = sp.misc.imread(elevationfile)
    (lx, ly, lz) = elevation.shape
    scalefactor = 0.5
    #print(elevation.shape)
    height = np.zeros(lx*ly).reshape(lx, ly)
    for i in range(lx):
        for j in range(ly):
            height[i][j] = elevation[i][j][0] * scalefactor
            #print(height[i][j])

    #plt.imshow(height, shape=(lx,ly))
    #plt.axis('off')
    #plt.show()


def main():
# parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-f", "--filename", help="data file")
    parser.add_argument("-k", "--kernels", help="kernel meta file")
    parser.add_argument("-e", "--elevation", help="elevation data file (image)")
    parser.add_argument("-s", "--saveplot", help="save plot file", action="store_true")
    args = parser.parse_args()
    
# define some global variables
    global bc, bw, bcw, chans, dchan, br

    global saveplot
    saveplot = args.saveplot
    
    bc = [325, 375, 425, 475, 525, 575, 625, 675, 725, 775]
    bw = 25
    chans = [150, 250, 350, 450, 550, 650, 750, 850, 950]
    dchan = 50
    br = [-20, -10]
    
# get file name (no directory)
    filename = args.filename
    global filenameonly;
    filenameonly = os.path.basename(filename)
    #filenameonly = os.path.splitext(filenameonly)[0]
    
# get elevation data file
    elevfile = args.elevation

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
    
# load elevation data
    if (elevfile):
        load_elevation(elevfile)

# plot results
    plot_results()


if __name__ == "__main__":
    main()


