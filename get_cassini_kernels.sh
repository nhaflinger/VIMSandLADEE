#!/bin/csh

python get_cassini_kernels.py -d https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/ -s ../../kernels/cassini_traj -p "SCPSE" -e "bsp"
