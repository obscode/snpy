#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards

fs = ['Y','J','H','K','Yc','Jc']
vegaB = standards['Vega']['vegaB4']

# JHK comes from Stritzinger et. al.,  Yc/Y comes from assumption.
mags = {'J':-0.001, 'H':0.000, 'K':-0.001, 'Jc':-0.001, 'Y':0.0, 'Yc':0.0}

for f in fs:
   print f,fset[f].compute_zpt(vegaB, mags[f])
