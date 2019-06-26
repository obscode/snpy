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


ids = concatenate([nonzero(diff(y))[0],array([len(y)-1])])
sts = y[ids]

N = len(sts)
n_col = int(ceil(sqrt(N)))
n_row = int(ceil(float(N)/n_col))
print N,n_col, n_row

for i in range(len(ids)):
   st = sts[i]
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
   tts = arange(xs[0], xs[-1]+0.5, 0.5)
   XX = zeros((tts.shape[0],2), dtype=float32)
   XX[:,0] = tts
   XX[:,1] = st
   ZZ = fl(XX)
   eZZ = efl(XX)
   mZZ = -isnan(ZZ)*-isnan(eZZ)*greater(ZZ,0)
   
   X = zeros((ts.shape[0],2), dtype=float32)
   X[:,0] = ts
   X[:,1] = st
   Z = fl(X)
   eZ = efl(X)
   mZ = -isnan(Z)*-isnan(eZ)*greater(Z,0)
   Z = where(mZ, Z, 1)
   eZ = where(mZ, eZ, -1)
   
   fig = plt.figure()
   ax1 = fig.add_axes([0.1,0.1,0.8,0.2])
   ax2 = fig.add_axes([0.1,0.3,0.8,0.6])
   ax2.errorbar(ts, vs, yerr=evs, fmt='o')
   ax2.plot(tts, ZZ, '-', color='red')
   ax2.plot(tts, ZZ+eZZ, '--', color='red')
   ax2.plot(tts, ZZ-eZZ, '--', color='red')
   ax1.set_xlabel('Epoch')
   ax2.set_ylabel('flux')
   ax2.set_title(name + " s = %.2f" % st)
   
   ax1.errorbar(ts, vs-Z, yerr=sqrt(evs**2 + eZ**2), fmt='o')
   ax1.axhline(0.0)
   
   plt.draw()
   plt.savefig(file.replace('.fits','')+"%.2f_%s.eps" % (st,name))

