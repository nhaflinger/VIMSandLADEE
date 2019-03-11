#
# sort_tracks.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: September 2017
# description: sort tracks by longitude 
#

import os
import argparse
import pickle
import numpy as np



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

def write_track(filehandle, shortname, range):
    # tracks in marea
    if (tlon[0] < range[1] and tlon[0] > range[0]):
        filehandle.write(shortname+"\n")
        #print(shortname, tlon[0], tlon[ng-1])


def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-d", "--datadir", help="data directory")
    parser.add_argument("-s", "--savefile", help="save text file")
    parser.add_argument("-r", "--range", nargs='+', type=float, help="longitude range")
    args = parser.parse_args()

    filename = args.savefile
    datadir = args.datadir
    lrange = args.range
    if (lrange[0] > lrange[1]):
        print("Range values in wrong order")
        quit()

    files = []
    filelist = []

    if (os.path.exists(datadir)):
        fh = open(filename, 'w')

        filelist = os.listdir(datadir)
        nfiles = len(filelist)
        for i in range(nfiles):
            fullname = datadir + "/" + filelist[i]
            shortname = os.path.splitext(filelist[i])[0][0:5]
            if (os.path.isfile(fullname)):
                ns = read_track(fullname)
                if (ns != 0):
                    ns = write_track(fh, shortname, lrange)

        fh.close()
    else:
        print("Data directory does not exist!")



if __name__ == "__main__":
    main()


