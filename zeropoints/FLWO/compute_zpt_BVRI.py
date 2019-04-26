#!/usr/bin/env python

# compute zeropoints for cfa3 4-shooter

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards,fset,standard_mags

fs = ['U4sh','B4sh','V4sh','R4sh','I4sh']
Uvega = standard_mags['Vega']['VegaB4']['Us']
Bvega = standard_mags['Vega']['VegaB4']['Bs']
Vvega = standard_mags['Vega']['VegaB4']['Vs']
Rvega = standard_mags['Vega']['VegaB4']['Rs']
Ivega = standard_mags['Vega']['VegaB4']['Is']
vegaB = standards['Vega']['VegaB']   # Bohlin version 5 (2004)

# Using the color terms from Hicken et al. (2009)
# to transform to CfA3 natural systems

vega_mags = {}
vega_mags['V4sh'] = Vvega + 0.0336*(Bvega - Vvega)
vega_mags['B4sh'] = vega_mags['V4sh'] + 0.8928*(Bvega - Vvega)
vega_mags['U4sh'] = vega_mags['B4sh'] + 0.9912*(Uvega - Bvega)
vega_mags['R4sh'] = vega_mags['V4sh'] - 1.0855*(Vvega - Rvega)
vega_mags['I4sh'] = vega_mags['V4sh'] - 1.0166*(Vvega - Ivega)


print("Using the following Natural magnitudes for Vega:")
for i,f in enumerate(fs):
   print("   ",f,vega_mags[f])

for i,f in enumerate(fs):
   print(f,fset[f].compute_zpt(vegaB, vega_mags[f]))
