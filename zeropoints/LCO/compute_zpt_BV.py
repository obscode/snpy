#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards,fset

fs = ['B','V','V0','V1']
Bvega = 0.03
Vvega = 0.03
vegaB = standards['Vega']['VegaB']   # Bohlin version 5 (2004)

# Need to convert these "standard" vega magnitudes into Natrual magnitudes using Contreras
#   color transformations:
# CAUTION:  we need a synthetic mag. for Vega in i', so this filters's zeropoint needs to 
#   be set FIRST!!!!

ip_vega = fset['i_40'].synth_mag(vegaB.wave, vegaB.resp)

# Use the color terms to transform to natural system from Contreras et al. (2010)
B_vega = Bvega - 0.061*(Bvega - Vvega)
V_vega = Vvega + 0.058*(Vvega - ip_vega)
# I derived this color term for the 'temp' V-band by fitting synthetic photometry
V1_vega = Vvega + 0.047*(Vvega - ip_vega)

vega_mags = [B_vega, V_vega, V_vega, V1_vega]

print "Using the following Natural (shifted) magnitudes for Vega:"
print "B = ",B_vega
print "V = ",V_vega
print "i' = ",ip_vega

for i in range(len(fs)):
   print fs[i],fset[fs[i]].compute_zpt(vegaB, vega_mags[i])
