#
# ladeeocctime_list.py
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

    print("Using Spice version ", spice_ver)


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
            if (ety[i][j] >= br[0] and ety[i][j] <= br[1]):
                qb.append(j)
        etq = np.asarray(qb)
        lq = len(etq)
        pmean = np.mean(profs[i][etq[0]:etq[lq-1]])
        profs[i] = profs[i] / pmean


# plot results
def plot_results(datafile):

    colors = ['grey','magenta','violet','blue','cyan','green','orange','red','maroon','brown']
    #colors=['violet','blue','cyan','orange','red','maroon']
    #colors=['grey','green','forest green', 'yellow','red','magenta','blue','cyan','turquoise']

    xr = [-10, 10]
    yr1 = [-0, 1.2]
    yr2 = [-0.035, 0.035]
    yr3 = [0, 100]
    
    #f, (ax1, ax2) = plt.subplots(2)
    f, ax2 = plt.subplots(1, sharex=True, sharey=False)
    f.suptitle(datafile, fontsize=12)
    
    nc = len(bc)
    #for i in range(nc):
        #lines = ax1.scatter(ety[i], altitude, s=1)
        #ax1.set_ylabel('Altitude (km)')
        #ax1.set_xlim(xr)

    nt = len(ets)
    pmean = np.zeros(nt)
    pmean = np.mean(profs[2:8], axis=0)
    #pmean = np.mean(profs[0:10], axis=0)

    for i in range(nc):
        #normed = ax2.plot(profs[2], profs[i] - profs[2])
        normed = ax2.plot(pmean, profs[i] - pmean)
     
        ax2.set_ylim(yr2)
        ax2.set_xlabel('Mean Transmission', fontsize=8)
        ax2.set_ylabel('Trans. Difference', fontsize=8)
        ax2.text(0.35, yr2[1]*1.10, "Wavelength Range (nm)", fontsize=8)
        
        for j in range(len(bc)):
            wstr = str(bc[j]-bw) + "-" + str(bc[j]+bw)
            ax2.text(-0.0*figdim+j*0.1, yr2[1]*1.03, wstr, fontsize=6, fontstyle='normal', color=colors[j])
        plt.setp(normed, color=colors[i], linewidth=1)
        
        #ax3 = ax2.twinx()
        #ax3.set_ylabel('Altitude (km)', fontsize=8)
        #ax3.tick_params('y')
        #ax3.scatter(profs[2], altitude, s=1)
        #ax3.set_ylim(yr3)
        
        #plt.tight_layout()

    #f.subplots_adjust(hspace=0.1)


# load processed uvs data
def readuvs_data(filename):
    global specs, wavels, stimes, ttimes, tellat, tellon, sunlat, sunlon, telaltx, solgraz
    
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
    telaltx = dataload['telaltx']
    solgraz = dataload['solgraze']

    global tlon, tlat, talt
    nt = len(tellon)
    tlon = np.zeros(nt)
    tlat = np.zeros(nt)
    talt = np.zeros(nt)
    for i in range(nt):
        tlon[i] = (float(tellon[i]) + 180.)%360 - 180
        tlat[i] = float(tellat[i])
        talt[i] = float(telaltx[i])

    global altitude
    lx = elevmap.shape[0]
    ly = elevmap.shape[1]
    #print(lx, ly)
    scalefactor = 0.5
    altitude = np.zeros(nt)
    for i in range(nt):
        lon = int(lx * 0.5 * (tlon[i] + 180.0) / 180.0 + 0.0)
        lat = int(ly * 0.5 * (tlat[i] + 90.0) / 90.0 + 0.0)
        #print(lon,lat)
        altitude[i] = elevmap[lon][lat][0] * scalefactor


# some of the wavelengths have spikes in the data.
# check for this and ignore that data
def clean_data():
    ns = len(specs)
    # only applying check for spikes for first 350 wavelengths
    #ns = 300
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


# load elevation and tracking data
def read_map(elevationfile):
    global elevmap
    
    if (os.path.exists(elevationfile) == False):
        print("Elevation map does not exist!")
        quit()
    
    elevmap = sp.misc.imread(elevationfile)
    #plt.imshow(elevmap, interpolation="nearest")


def main():
# parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-f", "--filelist", help="file names list")
    parser.add_argument("-d", "--datadir", help="data directory")
    parser.add_argument("-k", "--kernels", help="kernel meta file")
    parser.add_argument("-e", "--elevation", help="elevation data file (image)")
    parser.add_argument("-r", "--results", help="destination directory")
    parser.add_argument("-s", "--saveplot", help="save plot file", action="store_true")
    args = parser.parse_args()
    
# define some global variables
    global bc, bw, bcw, chans, dchan, br

    global saveplot, filenameonly, figdim
    
    bc = [325, 375, 425, 475, 525, 575, 625, 675, 725, 775]
    bw = 25
    chans = [150, 250, 350, 450, 550, 650, 750, 850, 950]
    dchan = 50
    br = [-20, -10]
    figdim = 9
    
# get args
    filename = args.filelist
    datadir = args.datadir
    kernelfile = args.kernels
    saveplot = args.saveplot
    elevfile = args.elevation
    results = args.results
    
    filenameonly = os.path.basename(filename)
        
    read_map(elevfile)
    
# list of tracks
    obs = []
    filein = open(filename, "r")
    for line in filein:
        obs.append(line.strip())
    
    nobs = len(obs)
    for i in range(nobs):
        datafile = datadir + "/" + obs[i] + ".pkl"
        #print(datafile)

        if (os.path.exists(datafile)):
# process uvs data
            readuvs_data(datafile)
            cspice_utc2et(kernelfile, stimes, ttimes)
            clean_data()
            fit_data()
            norm_data()
            plot_results(obs[i])

            if (saveplot):
                newdir = results + "/" + obs[i] 
                if (os.path.exists(newdir) == False):
                    os.mkdir(newdir)
                plotname = "transmission_" + obs[i] + ".pdf"
                print("Plotting results for " + plotname)
                plotfile = newdir + "/" + plotname 
                #plt.figure(i)
                plt.savefig(plotfile, figsize=(figdim,figdim), format='pdf')
                plt.close()

            else:
                plt.show()

    print("All tasks completed!")



if __name__ == "__main__":
    main()


