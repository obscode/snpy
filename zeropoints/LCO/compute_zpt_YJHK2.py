#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards

fs = ['Yd','Jd','Hd']

CTs = {'Jd':0.018,
       'Hd':-0.035,
       'Yd':0}

def YP98(J, K):
   '''Convert P98 J,Ks to Y'''
   JK = J-K
   return K + 2.398*JK - 4.545*JK**2 + 10.945*JK**3 - 13.321*JK**4 + 6.378*JK**5

def stand2lco(JP98,HP98,KP98):
   '''Convert photometry from standard to CSP, returns (Y,Ydw,J,Jrc2,H,K)'''
   JHP98 = (JP98-HP98)
   Y = YP98(JP98,KP98)
   J = JP98 - CTs['Jd']*JHP98
   H = HP98 - CTs['Hd']*JHP98
   return dict(Yd=Y,Jd=J,Hd=H)


stands = [standards['Vega']['VegaB4'], standards['CALSPEC']['ksi2ceti_stis_003']]
mags = {
   standards['Vega']['VegaB4']:(0.00,0.00,0.00),
   standards['CALSPEC']['ksi2ceti_stis_003']:(4.365, 4.38, 4.39) # Elias+1982
   }  

for stand in stands:
   print(stand)
   print("   standard JHK = ",mags[stand])
   smag = stand2lco(*mags[stand])
   for f in fs:
      print("      ",f,smag[f],fset[f].compute_zpt(stand, smag[f]))
