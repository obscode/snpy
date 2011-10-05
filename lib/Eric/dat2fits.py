#!/usr/bin/env python

import sys
import FITS
from Numeric import *
import string

file = sys.argv[1]
outfile = file+".fits"

f = open(file, 'r')
day = []
wave = []
flux = []
for line in f.readlines():
   fs = string.split(line)
   day.append(float(fs[0]))
   wave.append(float(fs[1]))
   flux.append(float(fs[2]))
f.close()

day = array(day);  wave = array(wave);   flux = array(flux)
# figure out the dimensions
val1 = day[0]
numwaves = len(nonzero(equal(day, day[0])))
nx = numwaves
x0 = wave[0]
dx = wave[1]-wave[0]
dy = day[numwaves] - day[0]
ny = len(day)/nx
y0 = day[0]

flux.shape = (ny,nx)

of = FITS.FITS(outfile, create=1)
of.newhead("T", -32, 2, [nx,ny])
of.putdata(flux)
print 1, x0, dx, 1, y0, dy
of.inskey("NAXIS2","CRPIX1",1)
of.inskey("CRPIX1","CRVAL1",x0)
of.inskey("CRVAL1","CDELT1",dx)
of.inskey("CDELT1","CRPIX2",1)
of.inskey("CRPIX2","CRVAL2",y0)
of.inskey("CRVAL2","CDELT2",dy)
of.close

