#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards,fset

vegaB = standards['Vega']['VegaB']

fs = ['Us','Bs','Vs','Rs','Is', 'g_m','r_m','i_m','z_m']

# These come from Alex Connley:
vega_mags = [0.0190, 0.0210, 0.0230, 0.0290, 0.0230, 0.0220, 0.0281, 
             0.0241, 0.0214]

for i in range(len(fs)):
   print fs[i],fset[fs[i]].compute_zpt(vegaB, vega_mags[i])
