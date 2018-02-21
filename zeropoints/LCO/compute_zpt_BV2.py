#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards,fset

fs = ['B2','V2']

CTs = {'u':(0.046, 0.017),
       'B':(0.061, 0.022),
       'g':(-0.014,0.011),
       'V':(-0.058,0.011),
       'r':(-0.016,0.015),
       'i':(-0.002,0.015),
       'J':(0.016, 0.00),
       'H':(-0.029, 0.00),
       'Ydw':(-0.042, 0.00),
       'u':(0.034, 0.001),
       'B':(0.085, 0.001),
       'g':(-0.010, 0.001),
       'V':(-0.060, 0.001),
       'r':(-0.013, 0.001),
       'i':(-0.002, 0.001),
       'u2':(0.034, 0.001),
       'B2':(0.085, 0.001),
       'g2':(-0.010, 0.001),
       'V2':(-0.060, 0.001),
       'r2':(-0.013, 0.001),
       'i2':(-0.002, 0.001)}

def stand2lco(B, V, i):
   '''Convert photometry from standard to CSP'''
   Bnat = B - CTs['B2'][0]*(B - V)
   dBnat = abs(CTs['B2'][1]*(B - V))
   Vnat = V - CTs['V2'][0]*(V - i)
   dVnat = abs(CTs['V2'][1]*(V - i))
   return ((Bnat,Vnat),(dBnat,dVnat))

stds = ['VegaB','p177d_stisnic_006','p330e_stisnic_007']
ip_vega = fset['i_40'].synth_mag(standards['VegaB'])
std_mags = [[0.030, 0.030, ip_vega],
            [14.138, 13.492, 13.170],
            [13.658, 13.028, 12.701]]

# Need to convert these "standard" vega magnitudes into Natrual magnitudes 
#  using Contreras color transformations:
# CAUTION:  we need a synthetic mag. for Vega in i', so this filters's 
#  zeropoint needs to  be set FIRST!!!!

nat_mags = [stand2lco(*m)[0] for m in std_mags]
d_nat_mags = [stand2lco(*m)[1] for m in std_mags]
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


