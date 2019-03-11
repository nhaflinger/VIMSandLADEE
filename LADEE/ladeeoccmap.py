#
# ladeeoccmap.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: September 2017
# description: displays occulation tracks on lunar map
#

import os
import argparse
import pickle, pprint
from numpy import sqrt, pi, exp, linspace, random
import numpy as np
import scipy as sp
from scipy import ndimage
import matplotlib.pyplot as plt


# read tracking data
def read_track(filename):
    global telllat, tellon, sunlat, sunlon, telalt, solgraz

    if (os.path.exists(filename) == False):
        print("Track file does not exist!")
    
    dataload = pickle.load(open(filename, "rb"))
    
    tellat = dataload['tellat']
    tellon = dataload['tellon']
    sunlat = dataload['sunlat']
    sunlon = dataload['sunlon']
    telalt = dataload['telalt']
    solgraz = dataload['solgraze']

    global slon, slat

    nt = len(sunlon)
    goodsol = []
    for i in range(nt):
        if (float(solgraz[i]) > 0.0 and float(solgraz[i]) < 0.5):
        #if (float(telalt[i]) < 1.0):
            goodsol.append(i)
    ns = len(goodsol)

    slon = np.zeros(ns)
    slat = np.zeros(ns)
    talt = np.zeros(ns)
    sgraze = np.zeros(ns)

    for i in range(ns):
        idx = goodsol[i]
        slon[i] = (float(sunlon[idx]) + 180.)%360 - 180
        slat[i] = float(sunlat[idx])
        talt[i] = float(telalt[idx])
        sgraze[i] = float(solgraz[idx])
        #print(idx, slon[i])

    return ns


# load elevation and tracking data
def read_map(elevationfile):
    global elevmap
    
    if (os.path.exists(elevationfile) == False):
        print("Elevation map does not exist!")
    
    elevmap = sp.misc.imread(elevationfile)

# plot results
def setup_plot():
    global axis1
    lx = elevmap.shape[0]
    ly = elevmap.shape[1]
    #print(lx,ly)

    scalefactor = 0.5
    height = np.zeros(lx*ly).reshape(lx, ly)
    for i in range(lx):
        for j in range(ly):
            height[i][j] = elevmap[i][j][0] * scalefactor

    latitude = np.linspace(-90, 90, 20-1)
    longitude = np.linspace(-180, 180, 20-1)

    f, axis1 = plt.subplots(1)
    axis1.imshow(elevmap, shape=(lx,ly), cmap='Greys', extent=[-180, 180, -90, 90])
    axis1.set_xticks(longitude)
    axis1.set_yticks(latitude)
    axis1.tick_params(labelsize=6)
    axis1.set_title("Occultation Tracks")
    axis1.set_xlabel("Longitude", fontsize=8)
    axis1.set_ylabel("Latitude", fontsize=8)


# lay down individual tracks
def plot_track():
    #print(slon[0], slat[0])
    ns = len(slat)
    axis1.text(slon[0], slat[ns-1]+8, shortname, fontsize=6, color='lawngreen', rotation='vertical')
    axis1.plot(slon, slat, color='lawngreen', linewidth=1.0)


def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-d", "--datadir", help="data directory")
    parser.add_argument("-e", "--elevation", help="elevation data file (image)")
    args = parser.parse_args()
    
    # get elevation data file
    elevfile = args.elevation
    read_map(elevfile)
    
    # get file name (no directory)
    datadir = args.datadir
    global filenameonly, shortname;
    filenameonly = os.path.basename(datadir)
    
    # plot results
    setup_plot()
    
    files = []
    filelist = []
    
    if (os.path.exists(datadir)):
        filelist = os.listdir(datadir)
        nfiles = len(filelist)
        for i in range(nfiles):
            fullname = datadir + "/" + filelist[i]
            shortname = os.path.splitext(filelist[i])[0][0:4]
            if (os.path.isfile(fullname)):
                #print("Plotting track for " + fullname)
                ns = read_track(fullname)
                if (ns != 0):
                    plot_track()
    else:
        print("Data directory does not exist!")
    
    plt.show()
    #plotfile = "occultation_tracks.pdf"
    #plt.savefig(plotfile, format='pdf')



if __name__ == "__main__":
    main()


