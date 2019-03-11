#
# ladeemap.py
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
    global telllat, tellon, telalt

    if (os.path.exists(filename) == False):
        print("Track file does not exist!")
    
    dataload = pickle.load(open(filename, "rb"))
    
    tellat = dataload['tellat']
    tellon = dataload['tellon']
    telalt = dataload['telalt']

    global tlon, tlat

    nt = len(tellon)
    goodtel = []
    for i in range(nt):
        #if (float(telalt[i]) <= 1.0):
        #if (float(telalt[i]) > 1.0):
        goodtel.append(i)
    ng = len(goodtel)

    tlon = np.zeros(ng)
    tlat = np.zeros(ng)
    talt = np.zeros(ng)

    for i in range(ng):
        idx = goodtel[i]
        tlon[i] = (float(tellon[idx]) + 180.)%360 - 180
        tlat[i] = float(tellat[idx])
        talt[i] = float(telalt[idx])

    return ng


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
    #print(tlon[0], tlat[0])
    ns = len(tlat)
    axis1.text(tlon[0], tlat[ns-1]+8, shortname, fontsize=6, color='lawngreen', rotation='vertical')
    axis1.plot(tlon, tlat, 'o', markersize=1, color='lawngreen')


def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-d", "--datadir", help="data directory")
    parser.add_argument("-t", "--track", help="specific track")
    parser.add_argument("-e", "--elevation", help="elevation data file (image)")
    parser.add_argument("-s", "--saveplot", help="save plot file", action="store_true")
    args = parser.parse_args()
    
    global saveplot
    saveplot = args.saveplot
    
    track = args.track.strip()
    
    # get elevation data file
    elevfile = args.elevation
    read_map(elevfile)
    
    # get file name (no directory)
    datadir = args.datadir
    global shortname;
    
    # plot results
    setup_plot()
    
    files = []
    filelist = []
    
    if (os.path.exists(datadir)):
        filelist = os.listdir(datadir)
        nfiles = len(filelist)
        for i in range(nfiles):
            fullname = datadir + "/" + filelist[i]
            shortname = os.path.splitext(filelist[i])[0][0:5]
         
            if (track == shortname):
                print("Plotting track for " + fullname)
                ns = read_track(fullname)
                if (ns != 0):
                    plot_track()
                break
            if (os.path.isfile(fullname) and track == ""):
                #print("Plotting track for " + fullname)
                ns = read_track(fullname)
                if (ns != 0):
                    plot_track()
    else:
        print("Data directory does not exist!")
    
    if (saveplot):
        plotfile = "occultation_tracks.pdf"
        plt.savefig(plotfile, format='pdf')
    else:
        plt.show()



if __name__ == "__main__":
    main()


