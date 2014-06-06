#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards
from snpy.filters import standard_mags
import scipy

# this function converts the USNO40 magnitudes to natural system magnitudes
def USNO2LCO(up, gp, rp, ip):
   '''Convert Photometric Telescope to CSP'''
   u = up - 0.045*(up - gp)
   g = gp + 0.015*(gp - rp)
   r = rp + 0.014*(rp - ip)
   i = ip - 0.002*(rp - ip)
   return(u,g,r,i)

# Filters to process
fs = ['u','g','r','i']

# dictionary of natural magnitudes
sm = {}

sloan_mags = standard_mags['Smith']
sloan_standards = standards['Smith']
# Convert Smith et al. magnitudes to LCO natural magnitudes
sm['u'],sm['g'],sm['r'],sm['i'] = USNO2LCO(sloan_mags['bd17']['u_40'], 
      sloan_mags['bd17']['g_40'], sloan_mags['bd17']['r_40'],
      sloan_mags['bd17']['i_40'])
print "Assuming standard and natural magnitdes:"
for f in fs:
   print f,"standard:  ",sloan_mags['bd17'][f+'_40'], ', natural:',sm[f]

# Do this three times, just to see who different SED's change the answer
print "zero-points with the original BD+174708:"
for f in fs:
   zpt = fset[f].compute_zpt(sloan_standards['bd17'], sm[f], zeropad=1)
   print f, zpt

print "zero-points with Mark's modified BD+174708:"
for f in fs:
   zpt = fset[f].compute_zpt(sloan_standards['bd17adj'], sm[f], 
         zeropad=1)
   print f, zpt

print "zero-points with Bohlin & Gilliland (2004) SED:"
for f in fs:
   zpt = fset[f].compute_zpt(sloan_standards['bd17B'], sm[f], 
         zeropad=1)
   print f, zpt
