#!/usr/bin/env python
'''A python module to generate lightcurve templates.  1st generation, based
on ATLAS training set for filters g+r and r+i, which we will call gr and ri.
Author:  Chris Burns
Version:  0.1

This module provides an object called template.  You create any number of these
and then simply ask it to generate a template of a given dm15 and then you can
evaluate that template at any time you like.  This is not a natural way to do
things, but it keeps it comaptible with the other dm15temp module I made for
Prieto's templates.

Another template, based on stretch, is also defined now and is based on the 
observation that (B-V) color curves have a well-defined maximum.  We define
the stretch as the time of this maximum divided by 30 days.

'''
from __future__ import print_function
import sys,os,string
import numpy as num
from scipy.interpolate import bisplrep,bisplev,splev,make_interp_spline
import scipy.optimize
import pickle
import json
try:
   from astropy.io import fits as pyfits
except ImportError:
   try:
      import pyfits
   except ImportError:
      sys.stderr.write('Error:  You need pyfits to run snpy.  You can get it\n')
      sys.stderr.write('        from:  http://www.stsci.edu/resources/'+\
                       'software_hardware/pyfits/\n')
      raise ImportError

debug=0

template_bands = ['ATgr','ATri']

base = os.path.dirname(globals()['__file__'])
if base == '':  base = '.'

# bivariate splines for light-curve generators
flux = {'st':{}}
eflux = {'st':{}}

def load_data(band, param='st', gen=1):
   '''Given a band param and generation of template generator, load the
   appropriate data and create a inter2D instance.'''
   if (band,gen) not in flux[param]:
      fil = os.path.join(base, 'fits', '{}_{}_mean{}.fits'.format(param,band,gen))
      pfile = fil.replace('.fits','.pickle')
      if debug:  print("Getting data from ",fil)
      if not os.path.isfile(fil) and not os.path.isfile(pfile):
         raise IOError("Could not find surface data for band %s, gen %d" \
               % (band,gen))
      # Look for previously saved pickle that's newer than the original FITS
      # file
      if os.path.isfile(pfile) and \
            os.path.getctime(pfile) > os.path.getctime(fil):
         f = open(pfile, 'rb')
         flux[param][(band,gen)] = pickle.load(f)
         f.close()
         f = open(pfile.replace('mean','std'), 'rb')
         eflux[param][(band,gen)] = pickle.load(f)
         f.close()
         return

      # pickle file does not exist, need to generate the bivriate spline interp
      f = pyfits.open(fil)
      h = f[0].header
      fdata = f[0].data
      f2 = pyfits.open(fil.replace('mean','std'))
      edata = f2[0].data
      f.close()
      f2.close()

      xs = h['CRVAL1'] + (num.arange(1,h['NAXIS1']+1) - h['CRPIX1'])\
            *h['CDELT1']
      ys = h['CRVAL2'] + (num.arange(1,h['NAXIS2']+1) - h['CRPIX2'])\
            *h['CDELT2']
      N = h['NAXIS1']*h['NAXIS2']
      x,y = num.meshgrid(xs,ys)
      tx = xs[::2]
      ty = ys[::2]
      flux[param][(band,gen)] = bisplrep(num.ravel(x), num.ravel(y), 
                                         num.ravel(fdata), task=-1, tx=tx, ty=ty)
      eflux[param][(band,gen)] = bisplrep(num.ravel(x), num.ravel(y), 
                                          num.ravel(edata), task=-1, tx=tx, ty=ty)
      # Try to save it as a pickle, to speed things up later
      try:
         f = open(pfile, 'wb')
         pickle.dump(flux[param][(band,gen)], f)
         f.close()
         f = open(pfile.replace('mean','std'), 'wb')
         pickle.dump(eflux[param][(band,gen)], f)
         f.close()
      except:
         pass

      if debug:
         rms = num.sqrt(num.mean(num.power(bisplev(xs,ys,df[(band,gen)])-fdata.T,2)))
         mad = num.median(num.absolute(bisplev(xs,ys,df[(band,gen)])-fdata.T))
         print("rms = ",rms, 'mad = ',mad)
   else:
      return

def poly(x, x0, coefs):
   return coefs[0] + coefs[1]*(x-x0) + coefs[2]/2*(x-x0)**2

def breakpoly(x, xb, coefs, before=True):
   if before:
      return coefs[0] + coefs[1]*(x-xb) + (x<xb)*coefs[2]*(x-xb)**2
   else:
      return coefs[0] + coefs[1]*(x-xb) + (x>xb)*coefs[2]*(x-xb)**2

def finterp(band, t, p, param, gen, extrap=False):
   '''interpolate at time t and param p for param,gen combo.'''
   if (band,gen) not in flux[param]:
      # Need to generate the bosplrep
      load_data(band,param,gen)
   # Select correct coefficients
   f = flux[param][(band,gen)]
   ef = eflux[param][(band,gen)]

   Z0 = 0    # Maybe in the future, we try a mean SN

   if len(num.shape(t)) == 0:
      scalar = 1
   else:
      scalar = 0
   t = num.atleast_1d(t)
   # First the evaluation mtarix:
   Z = num.atleast_2d(bisplev(t, p, f))[:,0] + Z0
   eZ = num.atleast_2d(bisplev(t, p, ef))[:,0]
   if not extrap:
      mask = num.greater_equal(t,f[0][0])*num.less_equal(t,f[0][-1])
      mask = mask*num.greater(Z, 0)
      Z = num.where(mask, Z, 1)
      eZ = num.where(mask, eZ, -1)
   else:
      t1,t2 = get_t_lim(band, param, gen)
      mask = num.logical_not(num.isnan(Z))
      # extrapolate lower with t^2 law
      if num.sometrue(num.less(t,t1)):
         Tp = bisplev(t1, p, f, dx=1)
         T = bisplev(t1, p, f)
         eT = bisplev(t, p, ef)
         t0 = t1 - 2*T/Tp; a = T/(t1-t0)**2
         Z = num.where(num.less(t, t1), a*num.power(t-t0,2), Z)
         eZ = num.where(num.less(t, t1), eT, eZ)
         mask = mask*num.greater(Z,0)*num.greater(t, t0)
      if num.sometrue(num.greater(t, t2)):
         # extrapolate with m = a*(t-t2)+b
         Tp = bisplev(t2, p, f, dx=1)
         T = bisplev(t2, p, f)
         eT = bisplev(t2, p, ef)
         b = -2.5*num.log10(T)
         a = -2.5/num.log(10)/T*Tp
         f = num.power(10, -0.4*(a*(t-t2)+b))
         Z = num.where(num.greater(t, t2), f, Z)
         eZ = num.where(num.greater(t, t2), eT, eZ)
      Z = num.where(mask, Z, 1)
   if scalar:
      return Z[0],eZ[0],mask[0]
   else:
      return Z,eZ,mask

def get_p_lim(band, param, gen):
   if (band,gen) not in flux[param]:
      # Need to generate the bosplrep
      load_data(band,param,gen)
   x0 = flux[param][(band,gen)][1][0]
   x1 = flux[param][(band,gen)][1][-1]

   return x0+0.00001, x1-0.00001

def get_t_lim(band, param, gen):
   '''Get the limit in time-space for trusting the templates.'''
   return -15.,40.
   

def dm152s(dm15):
   '''Convert from dm15 parameter to stretch.'''
   return((3.06-dm15)/2.04)

class st_template:
   def __init__(self):
      self.st = None
      self.normalize = 1   # Do we force max of lightcurve = 0?
      self.t = None
      self.ATgr = None;  self.eATgr = None
      self.ATri = None;  self.eATri = None

   def mktemplate(self, st, dm15_int=None, dm15_colors='int', generate=0):
      '''2nd and 3rd arguments ignored.'''
      self.st = st

      if generate:
         # generate the model light-curve from -10 to 80 in 1 day increments.
         self.t = num.arange(-15,40, 1.0)
         for band in ['ATgr','ATri']:
            self.__dict__[band],self.__dict__['e'+band], mask = \
               self.eval(band, self.t)
            self.__dict__['e'+band] = num.where(mask, self.__dict__['e'+band], -1.0)

   def deltaTmax(self, band, gen=1):
      if band == 'ATgr': return 0.23
      if band == 'ATri': return 0.38
      return 0

   def teval(self, band, time, st, gen=1, extrap=False):
      f,ef,mask = finterp(band, time, st, 'st', gen, extrap=extrap)

      return(num.where(mask,-2.5*num.log10(f),-1))

   def domain(self, band, gen=1):
      return get_t_lim(band, 'st', gen)

   def eval(self, band, times, z=0, mag=1, sextrap=1, gen=1, toff=True,
         extrap=False):
      '''Evaluate the template in band [band] at epochs [times].  Optionally
      redshift by (1+[z]).  If [mag]=1, return in magnitudes, otherwise return
      in flux units.  If [sextrap]=1, extrapolate beyond the training sample
      by using a stretch.  Use [gen] to specifiy the generation of the template.
      If you want the Tmax - Tmax(B) offset applied, set [toff] to True,
      otherwise, Tmax will be at 0 for every filter.'''

      if toff:
         evt = (times - self.deltaTmax(band))/(1+z)
      else:
         evt = times/(1+z)

      if self.st <= 0:
         return (evt*0, evt*0, num.zeros(evt.shape, dtype=num.bool))

      stmin,stmax = get_p_lim(band,'st', gen)
      if sextrap and not (stmin < self.st < stmax):
         if self.st < stmin:
            s = stmin/self.st
         else:
            s = stmax/self.st
      else:
         s = 1.0

      st = self.st
      if sextrap and st <= stmin:
         st = stmin + 0.001
      if sextrap and st >= stmax:
         st = stmax - 0.001

      tmin,tmax = get_t_lim(band, 'st', gen)
      tmask = num.greater_equal(evt, tmin)*num.less_equal(evt, tmax)

      # Apply any stretch
      evt = evt*s
      evd,eevd,mask = finterp(band, evt, st, 'st', gen, extrap=extrap)
      if not extrap:
         mask = mask*tmask

      if not mag:
         return(evd, eevd, mask)
      else:
         return(-2.5*num.log10(evd), eevd/evd*1.0857, mask)
