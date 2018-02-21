#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards
from snpy.filters import standard_mags
import scipy

CTs = {'u':(0.046, 0.017),
       'B':(0.061, 0.022),
       'g':(-0.014,0.011),
       'V':(-0.058,0.011),
       'r':(-0.016,0.015),
       'i':(-0.002,0.015),
       'J':(0.016, 0.00),
       'H':(-0.029, 0.00),
       'Ydw':(-0.042, 0.00),
       'u2':(0.034, 0.001),
       'B2':(0.085, 0.001),
       'g2':(-0.010, 0.001),
       'V2':(-0.060, 0.001),
       'r2':(-0.013, 0.001),
       'i2':(-0.002, 0.001)}


# this function converts the USNO40 magnitudes to natural system magnitudes
def USNO2LCO(up, gp, rp, ip):
   '''Convert Photometric Telescope to CSP'''
   u = up - CTs['u2'][0]*(up - gp)
   du = abs(CTs['u2'][1]*(up - gp))
   g = gp - CTs['g2'][0]*(gp - rp)
   dg = abs(CTs['g2'][1]*(gp - rp))
   r = rp - CTs['r2'][0]*(rp - ip)
   dr = abs(CTs['r2'][1]*(rp - ip))
   i = ip - CTs['i2'][0]*(rp - ip)
   di = abs(CTs['i2'][1]*(rp - ip))
   return((u,g,r,i),(du,dg,dr,di))

# Filters to process
fs = ['u2','g2','r2','i2']

# Standards to use
stds = ['bd17B','p177d_stisnic_006','p330e_stisnic_007']
std_mags = [[standard_mags['Smith']['bd17'][f[0]+'_40'] for f in fs],
            [15.128,13.771,13.316,13.170],
            [14.539,13.287,12.848,12.701]]

# Convert Standard mags to natural mags
nat_mags = [USNO2LCO(*m)[0] for m in std_mags]
d_nat_mags = [USNO2LCO(*m)[1] for m in std_mags]
zpts = {}
for f in fs:
   zpts[f] = []
for i,st in enumerate(stds):
   print "Assuming standard and natural magnitudes for %s:" % (st)
   for j,f in enumerate(fs):
      print "   %s  %6.3f    %6.3f +/- %5.3f" \
            % (f, std_mags[i][j], nat_mags[i][j], d_nat_mags[i][j])
   print "Zero-points:"
   for j,f in enumerate(fs):
      zpt = fset[f].compute_zpt(standards[st], nat_mags[i][j],
            zeropad=True)
      zpts[f].append(zpt)
      print "   %s  %7.4f +/- %6.4f" % (f, zpt, d_nat_mags[i][j])

for f in fs:
   print "std for %s:  %f" % (f, std(zpts[f]))
