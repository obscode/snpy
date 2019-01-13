#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards,fset
from snpy.filters import standard_mags
import scipy

vegaB = standards['Vega']['VegaB']
sloan_standards = standards['Smith']
sloan_mags = standard_mags['Smith']

fs = ['u_40','g_40','r_40','i_40','z_40']

# Note:  because the u'g'r'i'z' filters were defined with only 1.o airmass
#        and the standards were observed at 1.3 airmasses, we need to
#        modify the standard magnitudes:
sloan_mags['bd+174708']['u_40'] = 10.51
sloan_mags['bd+174708']['z_40'] = 9.22

# Smith et al. system is based on bd+174708, so just do that one, plus the
#   modified version that Mark came up with.
print("Here are the zero-points with the original bd+174708 spectrum:")
for f in fs:
   zpt = fset[f].compute_zpt(sloan_standards['bd+174708'], 
         sloan_mags['bd+174708'][f], zeropad=1)
   print(f, zpt)

print("here are the zero-points with Mark's adjusted bd+174708 spectrum:")
for f in fs:
   zpt = fset[f].compute_zpt(sloan_standards['bd+174708adj'], 
         sloan_mags['bd+174708adj'][f], zeropad=1)
   print(f, zpt)
