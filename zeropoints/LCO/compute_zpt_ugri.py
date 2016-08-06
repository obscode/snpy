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
   u = up - 0.046*(up - gp)
   g = gp + 0.014*(gp - rp)
   r = rp + 0.016*(rp - ip)
   i = ip + 0.002*(rp - ip)
   return(u,g,r,i)

# Filters to process
fs = ['u','g','r','i']

# Standards to use
stds = ['bd17B','p177d_stisnic_006','p330e_stisnic_007']
std_mags = [[standard_mags['Smith']['bd17'][f+'_40'] for f in fs],
            [15.118,13.750,13.302,13.161],
            [14.548,13.285,12.843,12.704]]

# Convert Standard mags to natural mags
nat_mags = [USNO2LCO(*m) for m in std_mags]
for i,std in enumerate(stds):
   print "Assuming standard and natural magnitdes for %s:" % (std)
   for j,f in enumerate(fs):
      print "   %s  %6.3f    %6.3f" % (f, std_mags[i][j], nat_mags[i][j])
   print "Zero-points:"
   for j,f in enumerate(fs):
      zpt = fset[f].compute_zpt(standards[std], nat_mags[i][j],
            zeropad=True)
      print "   %s  %7.4f" % (f, zpt)


# Do this three times, just to see who different SED's change the answer
#print "zero-points with the original BD+174708:"
#for f in fs:
#   zpt = fset[f].compute_zpt(sloan_standards['bd17'], sm[f], zeropad=1)
#   print f, zpt
#
#print "zero-points with Mark's modified BD+174708:"
#for f in fs:
#   zpt = fset[f].compute_zpt(sloan_standards['bd17adj'], sm[f], 
#         zeropad=1)
#   print f, zpt
#print "zero-points with Bohlin & Gilliland (2004) SED:"
#for f in fs:
#   zpt = fset[f].compute_zpt(sloan_standards['bd17B'], sm[f], 
#         zeropad=1)
#   print f, zpt
#
#print "zero-points with Bohlin & Gilliland (2004) SED:"
#for f in fs:
#   zpt = fset[f].compute_zpt(sloan_standards['bd17B'], sm[f], 
#         zeropad=1)
#   print f, zpt
#
#print "zero-points with Bohlin & Gilliland (2004) SED:"
#for f in fs:
#   zpt = fset[f].compute_zpt(sloan_standards['bd17B'], sm[f], 
#         zeropad=1)
#   print f, zpt
