#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import myplotlib
from numpy import *
import sys,os,string
from getSNdata2 import getdata
from scipy.interpolate import CloughTocher2DInterpolator as interp2D
#import FITS
from astropy.io import fits
import re

rcParams['font.size'] = 20

pat = re.compile(r'([ugriYJHBV])_mean2.fits')


f = open('BV_max.dat')
lines = f.readlines()
lines = map(string.strip, lines)
lines = map(string.split, lines)
names = [line[0] for line in lines]
name_sts = array([float(line[7])/30. for line in lines])

files = ['B_mean2.fits','Y_mean2.fits']
stmin = 0.2
stmax = 0.7
Nbins = 5

fig = myplotlib.PanelPlot(2,Nbins, figsize=(7,9))

for fi,file in enumerate(files):
   filt = pat.search(file).group(1)
   
   f = fits.open(file)
   h = f[0].header
   xs = h['CRVAL1'] + (arange(1,h['NAXIS1']+1) - h['CRPIX1'])*h['CDELT1']
   ys = h['CRVAL2'] + (arange(1,h['NAXIS2']+1) - h['CRPIX2'])*h['CDELT2']
   N = h['NAXIS1']*h['NAXIS2']
   X = dstack(meshgrid(xs,ys)).reshape((N,2))
   fl = interp2D(X, ravel(f[0].data))
   f.close()
   f = fits.open(file.replace('mean','std'))
   efl = interp2D(X,ravel(f[0].data))
   
   x,y,v,efluxes,mesh,xoff,xscale,yoff,yscale,mean_flux=getdata(filt, False)
   
   x = x*xscale + xoff
   y = y*yscale + yoff
   
   dst = (stmax - stmin)/Nbins
   
   ids = concatenate([nonzero(diff(y))[0],array([len(y)-1])])
   sts = y[ids]
   
   #fig = plt.figure(figsize=(6,9))
   #ax = fig.add_axes([0.1,0.1,0.8,0.8])
   for k in range(Nbins):
      stlow = stmin + k*dst
      sthigh = stlow + dst
   
      # The mean st for this bin
      tts = arange(xs[0], xs[-1]+0.5, 0.5)
      XX = zeros((tts.shape[0],2), dtype=float32)
      XX[:,0] = tts
      XX[:,1] = stlow + 0.5*dst
      ZZ = []
      eZZ = []
      mZZ = []
      ZZ.append(fl(XX))
      eZZ.append(efl(XX))
      mZZ.append(-isnan(ZZ)*-isnan(eZZ)*greater(ZZ,0))
         
      symbols = ['v','^','d','s','o']
      kk = k*2 + fi
      for i in range(len(ids)):
         st = sts[i]
         if st < stlow or st > sthigh:
            continue
         if st <= ys[0]:
            st = ys[0] + 0.0001
         if st >= ys[-1]:
            st = ys[-1] - 0.0001
         name = names[argmin(absolute(name_sts - st))]
         if i == 0:
            i1 = 0
            i2 = ids[i]+1
         else:
            i1 = ids[i-1]+1
            i2 = ids[i]+1
         
         ts = x[i1:i2]
         vs = v[i1:i2] + mean_flux
         evs = efluxes[i1:i2]
         XX[:,1] = st
         Z = fl(XX)
         eZ = efl(XX)
         #mZ = -isnan(Z)*-isnan(eZ)*greater(Z,0)
         #Z = where(mZ, Z, 1)
         #eZ = where(mZ, eZ, -1)
         ZZ.append(Z)
         
         fig.axes[kk].plot(ts, vs, symbols.pop(), color='k',
               label=name)
         
      if len(ZZ) > 1:
         print array(ZZ).shape
         ZZ = mean(array(ZZ[1:]), axis=0)
      else:
         ZZ = ZZ[0]
      fig.axes[kk].plot(tts[ZZ>0], ZZ[ZZ>0], '-', color='k')
      if not kk % 2:
         fig.axes[kk].legend(prop={'size':10}, frameon=False, numpoints=1,
               handlelength=0.75)
      else:
         fig.axes[kk].text(0.95,0.9, "$%.1f < s_{BV} < %.1f$" % \
            (stlow,sthigh), ha='right', va='top', 
            transform=fig.axes[kk].transAxes, fontsize=14)
      fig.axes[kk].yaxis.set_major_locator(MaxNLocator(4))
   
fig.xlabel('Epoch (days)', fontsize=22)
fig.ylabel('Normalized Flux', fontsize=22)
fig.set_limits()
for ax in fig.axes:
   ax.set_ylim(0.01,1.39)
fig.draw()

fig.fig.savefig('fast_BY_st.eps')
   
