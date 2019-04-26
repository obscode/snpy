#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards

fs = ['Y','Ydw','J','Jrc2','H','K']

CTs = {'J':0.039,
       'Jrc2':0.016,
       'H':-0.029,
       'Y':0,
       'Ydw':-0.042,
       'K':0}

def YP98(J, K):
   '''Convert P98 J,Ks to Y'''
   JK = J-K
   return K + 2.393*JK - 4.473*JK**2 + 10.715*JK**3 - 13.011*JK**4 + 6.232*JK**5

def stand2lco(JP98,HP98,KP98):
   '''Convert photometry from standard to CSP, returns (Y,Ydw,J,Jrc2,H,K)'''
   JHP98 = (JP98-HP98)
   Y = YP98(JP98,KP98)
   Ydw = Y - CTs['Ydw']*JHP98
   J = JP98 - CTs['J']*JHP98
   Jrc2 = JP98 - CTs['Jrc2']*JHP98
   H = HP98 - CTs['H']*JHP98
   K = KP98
   return dict(Y=Y,Ydw=Ydw,J=J,Jrc2=Jrc2,H=H,K=K)


stands = [standards['Vega']['VegaB4'], standards['CALSPEC']['ksi2ceti']]
mags = {
   standards['Vega']['VegaB4']:(0.00,0.00,0.00),
   standards['CALSPEC']['ksi2ceti']:(4.365, 4.38, 4.39) # Elias+1982
   }  

for stand in stands:
   print(stand)
   print("   standard JHK = ",mags[stand])
   smag = stand2lco(*mags[stand])
   for f in fs:
      print("      ",f,smag[f],fset[f].compute_zpt(stand, smag[f]))
