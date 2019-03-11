#
# get_cassini_kernels.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: October 2017
# description: download Cassini trajectory kernel files
#

import os
import argparse
import wget
from urllib.request import Request, urlopen, urlretrieve
from bs4 import BeautifulSoup

def read_url(url):
    urllist = []
    url = url.replace(" ","%20")
    req = Request(url)
    a = urlopen(req).read()
    soup = BeautifulSoup(a, 'html.parser')
    x = (soup.find_all('a'))
    for i in x:
        file_name = i.extract().get_text()
        url_new = url + file_name
        url_new = url_new.replace(" ","%20")
        if(file_name[-1]=='/' and file_name[0]!='.'):
            read_url(url_new)

        #print(url_new)
        urllist.append(url_new)
    return urllist


def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-e", "--fileext", help="File extension to match")
    parser.add_argument("-p", "--pattern", help="Pattern to match")
    parser.add_argument("-d", "--datadir", help="PDS data directory (URL)")
    parser.add_argument("-s", "--savedir", help="destination directory")
    args = parser.parse_args()

    # get list of files under URL
    datadir = args.datadir
    savedir = args.savedir
    fileext = args.fileext
    pattern = args.pattern

    urllist = read_url(datadir)
    #print(urllist)

    for line in urllist:
        url = line.strip()
        file_parts = url.split(".")
        file_segs = url.split("_")
        ext = len(file_parts)-1
        if (file_parts[ext] == fileext):
            if (file_segs[1] == pattern):
                print("\n")
                print(url)
                filename = wget.download(url, out=savedir)



if __name__ == "__main__":
    main()


