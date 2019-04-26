#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards,fset

fs = ['Bk','Bk2','Vk','Vk2']
Bvega = 0.03
Vvega = 0.03
vegaB = standards['Vega']['VegaB']   # Bohlin version 5 (2004)

# Using the color terms from Hicken et al. (2011)
# to transform to CfA3/4 natural systems

V_vega = Vvega + 0.0233*(Bvega - Vvega)
B_vega1 = V_vega + 0.9294*(Bvega - Vvega)
B_vega2 = V_vega + 0.8734*(Bvega - Vvega)

vega_mags = [B_vega1, B_vega2, V_vega, V_vega]

print("Using the following Natural (shifted) magnitudes for Vega:")
print("Bk = ",B_vega1)
print("Bk2 = ",B_vega2)
print("Vk = ",V_vega)
print("Vk2 = ",V_vega)

for i in range(len(fs)):
   print(fs[i],fset[fs[i]].compute_zpt(vegaB, vega_mags[i]))
