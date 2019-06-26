#!/usr/bin/env python
'''A python module to generate 2nd generation lightcurve templates.
Author:  Chris Burns
Version:  0.2 

In this generation of templates, we still fit a smooth surface defined by
(t-Tmax, dm15/s, m-Mmax) data points from our best observed SNeIa.  The
difference is that we have determined the surface using Gaussian Processes.

The details of this are not important.  The template surface is not tabulated
in a FITS file and interpolated using the new Clauchy-Tocher interpolator,
which fits bi-variate cubic splines on the convex hulls of the data.

We still keep the same class structure as the old templates (maybe even
incorporate the old ones in for back-compatibility...

This module provides an object called template.  You create any number of these
and then simply ask it to generate a template of a given dm15 and then you can
evaluate that template at any time you like.  This is not a natural way to do
things, but it keeps it comaptible with the other dm15temp module I made for
Prieto's templates.

Example:
>>> import dm15temp2
>>> t = dm15temp2.template()
>>> t.mktemplate(1.6)
>>> flux,e_flux,mask = t.eval(arange(-10, 70, 1.0), 'B')

Here, we create an instance of the template, ask it to generate a lc template
with dm15 = 1.6, then we evaluate the flux on the interval [-10,70] in steps
of 1 day.  It returns a 3-tuple:  flux (normalized to peak of 1.0), error in
flux and a mask.  The mask is 1 where the template is well-behaved (you could,
afterall ask for a template that goes to day 1000, but its mask value will be 0,
because we have no constraining data there).

This module can also be run as a standalone script:

   ./dm15temp2 [dm15] [filter]

in which case, the template is sent to starndard out.
'''
import sys,os,string,re
import numpy as num
import FITS
from scipy.interpolate import CloughTocher2DInterpolator as interp2d
import scipy.optimize

base = os.path.dirname(globals()['__file__'])
if base == '':  base = '.'
# Cache dictionaries
mean = {}
std = {}

filepat = re.compile(r'(.*)_mean([0-9])\.fits')
def load_filter(filter, gen=2):
   '''Given the filter and generation, attempt to load the data from
   FITS file and put it in the cache.'''
   file = os.path.join(base,"%s_mean%d.fits" % (filter,gen))
   f = FITS.FITS(file)
   xs = f['CRVAL1'] + (num.arange(f['NAXIS1']) + 1 - f['CRPIX1'])*f['CDELT1']
   ys = f['CRVAL2'] + (num.arange(f['NAXIS2']) + 1 - f['CRPIX2'])*f['CDELT2']
   xx,yy = num.meshgrid(xs,ys)
   points = num.array([num.ravel(xx),num.ravel(yy)]).T
   values = f.data()
   mean[(filter,gen)] = interp2d(points, num.ravel(values))
   f.close()
   f = FITS.FITS(file.replace('mean','std'))
   values = f.data()
   std[(filter,gen)] = interp2d(points,num.ravel(values))
   f.close()

def dm152s(dm15):
   '''Convert from dm15 parameter to stretch.'''
   return((3.06-dm15)/2.04)

class stretch_template:
   #def __init__(self, flux=True):
   #   self.flux = flux   # Do we return fluxes or magnitudes?

   def _find_max(self, filter, s, gen=2):
      '''Given the filter and stretch, find where maximum occurs and value.'''
      # An initial evaluation
      ts = num.arange(-10,20,1.0)
      fs,efs,gids = self.__call__(filter,s, ts, gen,
            normalize_flux=False, normalize_time=False, sextrap=True)
      dfs = num.diff(fs[gids])
      id = num.nonzero(num.less(dfs,0))[0][0]
      brack = (ts[gids][id-1],ts[gids][id],ts[gids][id+1])
      f = lambda x: -self.__call__(filter,s,x,gen,normalize_flux=False,
            normalize_time=False)[0]
      tmax,fmax,iter,fcs = scipy.optimize.brent(f,brack=brack,full_output=True)
      return tmax,-fmax

   def __call__(self, filter, s, t, gen=2, normalize_flux=True, 
         normalize_time=True, sextrap=True):
      '''Given the filter, stretch, and evaluation times, return interpolated
      fluxes.'''
      global mean,std
      if (filter,gen) not in mean:
         load_filter(filter,gen)
      if normalize_time or normalize_flux:
         tmax,fmax = self._find_max(filter,s,gen)
      else:
         tmax = 0
         fmax = 1

      m = mean[(filter,gen)]
      v = std[(filter,gen)]
      min_st = m.points[:,1].min()
      max_st = m.points[:,1].max()
      if len(num.shape(t)) == 0:
         scalar = True
         t = num.array([t])
      else:
         scalar = False
      t = t - tmax
      if s > max_st and sextrap:
         evalpts = num.array([t*max_st/s,t*0+max_st]).T
      elif s < min_st and sextrap:
         evalpts = num.array([t*min_st/s,t*0+min_st]).T
      else:
         evalpts = num.array([t,t*0+s]).T
      flux = m(evalpts)/fmax
      vflux = v(evalpts)
      gids = -num.isnan(flux)
      flux[-gids] = 0
      vflux[-gids] = 0
      if scalar:
         return flux[0],vflux[0],gids[0]
      else:
         return flux,vflux,gids
