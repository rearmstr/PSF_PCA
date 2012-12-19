#!/usr/bin/env python
#
# make sure to run "setup scipy" and "setup esutil -r ~esheldon/exports/esutil-local"
# before excute this script.
#
# run the script by "cat xy.dat | image2sky-example.py" where xy.dat list ccd, x, and y.
#
from sys import stdin
import esutil
import cPickle as pickle
# import scipy
# import numpy
wcs_example_file='/astro/u/esheldon/oh/desfiles/wcsexample/wcs-decam--24--15-i-6.pickle'
# wcs_example_file='wcs-decam--24--15-i-6.pickle'

with open(wcs_example_file) as wcsfile:
    many_wcs = pickle.load(wcsfile)

for line in stdin:
    ls = line.split()

    ccd = int(ls[0])
    x = float(ls[1])
    y = float(ls[2])

    ra,dec = many_wcs[ccd].image2sky(x,y)
    print ccd,x,y,ra,dec
