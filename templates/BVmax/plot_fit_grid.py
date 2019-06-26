#!/usr/bin/env python

from matplotlib import pyplot as plt
from matplotlib import rcParams
from numpy import *
import sys,os,string
from getSNdata2 import getdata
from scipy.interpolate import CloughTocher2DInterpolator as interp2D
import FITS
import re

pat = re.compile(r'([ugriYJHBV])_mean2.fits')
f = open('BV_max.dat')
lines = f.readlines()
lines = map(string.strip, lines)
lines = map(string.split, lines)
names = [line[0] for line in lines]
name_sts = array([float(line[7])/30. for line in lines])

file = sys.argv[1]
if len(sys.argv) > 2:
   Nbins = int(sys.argv[2])
else:
   Nbins = 20
filt = pat.search(file).group(1)

f = FITS.FITS(file)
xs = f['CRVAL1'] + (arange(1,f['NAXIS1']+1) - f['CRPIX1'])*f['CDELT1']
ys = f['CRVAL2'] + (arange(1,f['NAXIS2']+1) - f['CRPIX2'])*f['CDELT2']
N = f['NAXIS1']*f['NAXIS2']
X = dstack(meshgrid(xs,ys)).reshape((N,2))
fl = interp2D(X, ravel(f.data()))
efl = interp2D(X,ravel(FITS.qread(file.replace('mean','std'))))

x,y,v,efluxes,mesh,xoff,xscale,yoff,yscale,mean_flux = getdata(filt, False)

x = x*xscale + xoff
y = y*yscale + yoff

stmin = y.min()
stmax = y.max()
dst = (stmax - stmin)/Nbins


ids = concatenate([nonzero(diff(y))[0],array([len(y)-1])])
sts = y[ids]

N = len(sts)

for k in range(Nbins):
   stlow = stmin + k*dst
   sthigh = stlow + dst

   plt.autoscale(enable=True)
   fig = plt.figure()
   ax1 = fig.add_axes([0.1,0.1,0.8,0.2])
   ax2 = fig.add_axes([0.1,0.3,0.8,0.6])

   tts = arange(xs[0], xs[-1]+0.5, 0.5)
   XX = zeros((tts.shape[0],2), dtype=float32)
   XX[:,0] = tts
   XX[:,1] = stlow + 0.5*dst
   ZZ = fl(XX)
   eZZ = efl(XX)
   mZZ = -isnan(ZZ)*-isnan(eZZ)*greater(ZZ,0)

   for i in range(len(ids)):
      st = sts[i]
      if st < stlow or st > sthigh:  continue
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
      
      X = zeros((ts.shape[0],2), dtype=float32)
      X[:,0] = ts
      X[:,1] = st
      Z = fl(X)
      eZ = efl(X)
      mZ = -isnan(Z)*-isnan(eZ)*greater(Z,0)
      Z = where(mZ, Z, 1)
      eZ = where(mZ, eZ, -1)
   
      ax2.errorbar(ts, -2.5*log10(vs), yerr=evs/vs, fmt='o', label=name)
   
      ax1.errorbar(ts, -2.5*log10(vs)+2.5*log10(Z), yerr=evs/vs, fmt='o')
   
   ax1.axhline(0.0)
   ax1.set_xlabel('Epoch')
   ax2.set_ylabel('flux')
   ax2.set_title(name + " s = %.2f" % (stlow+0.5*dst))
   ax2.plot(tts, -2.5*log10(ZZ), '-', color='red')
   plt.autoscale(enable=False)
   ax2.plot(tts, -2.5*log10(ZZ+eZZ), '--', color='red')
   ax2.plot(tts, -2.5*log10(ZZ-eZZ), '--', color='red')
   ax2.invert_yaxis()
   ax2.legend(prop={'size':14}, numpoints=1, frameon=False)
   ax1.set_xlim(ax2.get_xlim())
   ax1.set_ylim(-0.1,0.1)
   ax2.set_xticklabels([])
   plt.draw()
   plt.savefig(file.replace('.fits','')+"_%02d.eps" % (k))

