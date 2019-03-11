#
# get_vims_occultations.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: October 2017
# description: download Cassini VIMS ring occultation data
#

import os
import argparse
import wget


def main():
# parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-f", "--filelist", help="text file with list of files")
    parser.add_argument("-d", "--datadir", help="PDS data directory (URL)")
    parser.add_argument("-s", "--savedir", help="destination directory")
    args = parser.parse_args()

    # read list of data files containing ring occultation data
    filename = args.filelist
    datadir = args.datadir
    savedir = args.savedir

    filein = open(filename, "r")

    for line in filein:
        url = datadir + "/" + line.strip()
        #print(url)
        filename = wget.download(url, out=savedir)

    filein.close()



if __name__ == "__main__":
    main()


