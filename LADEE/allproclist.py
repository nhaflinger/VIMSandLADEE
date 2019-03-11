#
# allproclist.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: August 2017
# description: processing LADEE UVS data
#

import os
import argparse
import re
import pickle

import xml.etree.ElementTree as et
import numpy as np


# compute this data
# specs, wavels, stimes, ttimes, solgraze, solalt, telalt, telsol, solsol, telaltx, tellat, tellon, sunlat, sunlon, tsolaz, tsolel, tscx,tscy, tscz

def readuvs(datadir, filename):
    filedir = datadir + "/" + filename
    files = []
    filelist = []
    
    if (os.path.exists(filedir)):
        #print("Directory " + filedir + " exists!")
        filelist = os.listdir(filedir)
        nfiles = len(filelist)
        for i in range(nfiles):
            #print(filelist[i])
            fullname = filedir + "/" + filelist[i]
            if (os.path.isfile(fullname)) :
                if (".TAB" in fullname):
                    files.append(fullname)

    nf = len(files)
    global specs
    specs = np.zeros(nf*1024).reshape(nf, 1024)
    specx = np.zeros(1024)

# store radiance values in array specs
    for i in range(nf):
        filein = open(files[i], "r")
        #print("Computing specs for" + files[i])
        
        nline = 0
        for line in filein:
            specx[nline] = line.strip()
            specs[i][nline] = specx[nline]
            #print("DEBUG: " + str(i) + " " + str(nline) + " " + str(specs[i][nline]))
            nline += 1
        
        filein.close()

    filex = []
    filelist = []
    
    if (os.path.exists(filedir)):
        #print("Directory " + filedir + " exists!")
        filelist = os.listdir(filedir)
        nfiles = len(filelist)
        for i in range(nfiles):
            #print(filelist[i])
            fullname = filedir + "/" + filelist[i]
            if (os.path.isfile(fullname)) :
                if (".XML" in fullname):
                    filex.append(fullname)
        
    global stimes, ttimes, solgraze, solalt, telalt, telsol, solsol, telaltx, solgrazex, tellat, tellon, sunlat, sunlon, tsolaz, tsolel, tscx, tscy, tscz
    stimes = []
    ttimes = []
    solgraze = []
    solgrazex = []
    solalt = []
    telalt = []
    telaltx = []
    telsol = []
    solsol = []
    tellat = []
    tellon = []
    sunlat = []
    sunlon = []
    tsolaz = []
    tsolel = []
    tscx = []
    tscy = []
    tscz = []
    
# parse xml files
    for i in range(nf):
        #print("DEBUG: parsing " + filex[i])
        tree = et.iterparse(filex[i])
        for event,elem in tree:
            tag = elem.tag
            value = elem.text
            if ("start_date_time" in tag):
                stimes.append(value)
                #print(tag + ":" + stimes[i])
            if ("stop_date_time" in tag):
                ttimes.append(value)
            if ("sun_graze_altitude_above_terrain" in tag):
                solgraze.append(value)
            if ("sun_graze_altitude" in tag):
                solgrazex.append(value)
            if ("solar_viewer_graze_altitude_above_terrain" in tag):
                solalt.append(value)
            if ("telescope_graze_altitude_above_terrain" in tag):
                telaltx.append(value)
            if ("telescope_fov_graze_altitude_above_terrain" in tag):
                telalt.append(value)
            if ("telescope_solar_elongation" in tag):
                telsol.append(value)
            if ("solar_viewer_solar_elongation" in tag):
                solsol.append(value)
            if ("telescope_graze_latitude" in tag):
                tellat.append(value)
            if ("telescope_graze_longitude" in tag):
                tellon.append(value)
            if ("solar_viewer_graze_latitude" in tag):
                sunlat.append(value)
            if ("solar_viewer_graze_longitude" in tag):
                sunlon.append(value)
            if ("telescope_sun_azimuth" in tag):
                tsolaz.append(value)
            if ("telescope_sun_elevation" in tag):
                tsolel.append(value)
            if ("moon_fixed_x" in tag):
                tscx.append(value)
            if ("moon_fixed_y" in tag):
                tscy.append(value)
            if ("moon_fixed_z" in tag):
                tscz.append(value)



# store wavelength calibration data
def readuvs_wavelength(filename):
    global wavels
    wavels = np.zeros(2*1024).reshape(1024,2)
    
    filein = open(filename, "r")
    
    nline = 0
    for line in filein:
        wavedata = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        wavels[nline][0] = wavedata[0]
        wavels[nline][1] = wavedata[1]
        #print("DEBUG: " + str(wavels[nline][0]) + " " + str(wavels[nline][1]))
        nline += 1
    
    filein.close()


def main():
# parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-f", "--filelist", help="input file names")
    parser.add_argument("-d", "--datadir", help="data directory")
    parser.add_argument("-s", "--savedir", help="save directory")
    parser.add_argument("-w", "--wavelengths", help="wavelengths data")
    args = parser.parse_args()
    
# read wavelength calibration data file
    filename = args.wavelengths
    readuvs_wavelength(filename)

# read list of data files containing radiance data
    obs = []
    filename = args.filelist
    datadir = args.datadir
    filein = open(filename, "r")

    for line in filein:
        obs.append(line.strip())

    nobs = len(obs)
    for i in range(nobs):
        #print("Loading " + obs[i])
        if (obs[i] == "0075R" or obs[i] == "0486F"):
            #print("Ignoring " + obs[i])
            continue

        readuvs(datadir, obs[i])
    
# save data files
        filedir = datadir + "/" + obs[i]
        if (os.path.exists(filedir)):
            print("Dumping data for directory " + obs[i])
            datadump = {'specs':specs, 'wavels':wavels, 'stimes':stimes, 'ttimes':ttimes, 'solgraze':solgraze,'solgrazex':solgrazex, 'solalt':solalt, 'telalt':telalt, 'telsol':telsol, 'solsol':solsol,'telaltx':telaltx, 'tellat':tellat, 'tellon':tellon, 'sunlat':sunlat, 'sunlon':sunlon, 'tsolaz':tsolaz, 'tsolel':tsolel, 'tscx':tscx, 'tscy':tscy, 'tscz':tscz}

            pckl_file = open(args.savedir + "/" + obs[i] + ".pkl", "wb")
            pickle.dump(datadump, pckl_file)

    filein.close()



if __name__ == "__main__":
    main()

