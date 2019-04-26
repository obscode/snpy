#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards,fset,standard_mags
import scipy

sloan_standards = standards['Smith']
sloan_mags = standard_mags['Smith']

def PT2Img(up, gp, rp, ip, zp):
   '''Convert Photometric Telescope to 2.5m imager.'''
   u = up
   g = gp + 0.060*((gp - rp) - 0.53)
   r = rp + 0.035*((rp - ip) - 0.21)
   i = ip + 0.041*((rp - ip) - 0.21)
   z = zp - 0.030*((ip - zp) - 0.09)
   return(u,g,r,i,z)

fs = ['u_s','g_s','r_s','i_s','z_s']
sm = {}

stand = 'bd+174708'
sm['u_s'],sm['g_s'],sm['r_s'],sm['i_s'],sm['z_s'] = \
      PT2Img(sloan_mags[stand]['u_40'], sloan_mags[stand]['g_40'], 
      sloan_mags[stand]['r_40'], sloan_mags[stand]['i_40'],
      sloan_mags[stand]['z_40'])

print("The zero-points for the orignal BD+174708 SED:")
for f in fs:
   zpt = fset[f].compute_zpt(sloan_standards['bd+174708'], 
         sm[f],zeropad=1)
   print(f,zpt)
print("The zero-points for Mark's modified BD+174708 SED:")
for f in fs:
   zpt = fset[f].compute_zpt(sloan_standards['bd+174708adj'], 
         sm[f],zeropad=1)
   print(f,zpt)

