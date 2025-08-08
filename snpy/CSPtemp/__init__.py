#!/usr/bin/env python
'''A python module to generate lightcurve templates.  2nd generation, based
on DR2
Author:  Chris Burns
Version:  0.2

This module provides an object called template.  You create any number of these
and then simply ask it to generate a template of a given dm15 and then you can
evaluate that template at any time you like.  This is not a natural way to do
things, but it keeps it comaptible with the other dm15temp module I made for
Prieto's templates.

Another template, based on stretch, is also defined now and is based on the 
observation that (B-V) color curves have a well-defined maximum.  We define
the stretch as the time of this maximum divided by 30 days.

Example:
>>> import CSPtemp
>>> t = CSPtemp.dm15_template()
>>> t.mktemplate(1.6)
>>> flux,e_flux,mask = t('B', arange(-10, 70, 1.0))

>>> t = CSPtemp.s_template()
>>> t.mktemplate(1.1)
>>> flux,e_flux,mask = t('B', arange(-10, 70, 1.0))

Here, we create an instance of the template, ask it to generate a lc template
with dm15 = 1.6, then we evaluate the flux on the interval [-10,70] in steps
of 1 day.  It returns a 3-tuple:  flux (normalized to peak of 1.0), error in
flux and a mask.  The mask is 1 where the template is well-behaved (you could,
afterall ask for a template that goes to day 1000, but its mask value will be 0,
because we have no constraining data there).

NEW:  Using GLOES is rather slow and expensive.  However, you can always 
generate the surface and fit it with a bi-variate spline, capturing the 
information, but making things far faster.  Now, we check to see if the tck
in a file tck.pickle is available.  If so, then use this instead of calling
gloes.

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

template_bands = ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']

base = os.path.dirname(globals()['__file__'])
if base == '':  base = '.'

NIR_range = [-12,10]

# bivariate splines for light-curve generators
flux = {'st':{}, 'dm15':{}}
eflux = {'st':{}, 'dm15':{}}
# Mean LC generation 3 LC generators
meanflux = None

# Spines for Tmax(\lambda) - Tmax(B) for generation 3 LC generators
T_knots = [0.3, 0.5, 0.8, 1.0, 1.2]
coefs = {
   'u':[-2.11940249,-2.32966301,-1.48594271,-1.31983376,-1.07048405],
   'g':[ 0.64812471, 0.64715437, 0.51726053, 0.33568555, 0.29665327],
   'r':[ 3.37679485, 1.94205847, 1.91287074, 2.01362703, 2.73100208],
   'i':[ 5.23293714, 1.16878210,-2.46636081,-2.94255067,-2.69343348],
   'V':[ 1.76026462, 2.15903867, 1.83847808, 1.71077089, 1.46472528],
   'Y':[ 5.84131387, 2.05615124,-3.36406763,-4.15808803,-4.58703313],
   'J':[ 3.92132999, 1.14454630,-2.96897932,-3.69300855,-3.56201202],
   'H':[ 4.16136404, 0.87414653,-3.31381139,-4.18343476,-3.76891347]}
Tmax_splines = {}
for filt in coefs:
   Tmax_splines[filt] = make_interp_spline(T_knots, coefs[filt], bc_type='clamped')

def load_data(band, param='dm15', gen=1):
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

def load_meanSN():
   '''Load in the mean SN template, for use with generation 3 templates.'''
   global meanflux
   fil = os.path.join(base, 'SN2012fr_mean.json')
   pfil = fil.replace('json','pkl')
   if os.path.isfile(pfil) and os.path.getctime(pfil) > os.path.getctime(fil):
      with open(pfil, 'rb') as fin:
         meanflux = pickle.load(fin)
      return
   with open(fil, 'r') as fin:
      meanflux = json.load(fin)
   for filt in meanflux:
      meanflux[filt]['tck'] = (num.array(meanflux[filt]['tck'][0]),
                               num.array(meanflux[filt]['tck'][1]),
                               meanflux[filt]['tck'][2])
   # now safe as pickle for faster work
   with open(pfil, 'wb') as fout:
      pickle.dump(meanflux, fout)
   return

def meanSN(filt, t, s):
   '''Given time relative to filter maximum, t, and stretch, return the stretched
   mean SN light-curve.'''
   if meanflux is None:
      load_meanSN()
   di = meanflux[filt]
   xx = t/s
   F = splev(xx,di['tck'])
   F = num.where(xx < di['t0'], di['f0']*(xx - di['te'])**2, F)
   F = num.where(xx > di['t1'], di['f1']*num.exp(-xx/di['tau']), F)
   F = num.where(xx < di['te'], 0, F) 
   return F


def poly(x, x0, coefs):
   return coefs[0] + coefs[1]*(x-x0) + coefs[2]/2*(x-x0)**2

def breakpoly(x, xb, coefs, before=True):
   if before:
      return coefs[0] + coefs[1]*(x-xb) + (x<xb)*coefs[2]*(x-xb)**2
   else:
      return coefs[0] + coefs[1]*(x-xb) + (x>xb)*coefs[2]*(x-xb)**2

def finterp(band, t, p, param, gen, extrap=False):
   '''interpolate at time t and param p for param,gen combo.'''
   global meanflux
   if (band,gen) not in flux[param]:
      # Need to generate the bosplrep
      load_data(band,param,gen)
   #if gen > 2 and meanflux is None:
   #   # Need to load up the mean SN for generation 3 and beyond
   #   load_meanSN()

   # Select correct coefficients
   f = flux[param][(band,gen)]
   ef = eflux[param][(band,gen)]

   #if gen > 2:
   #   z0 = 0
   #   # Need to construct a mean LC by stretching SN2012fr
   #   if param == 'dm15':
   #      s = dm152s(p)
   #   else:
   #      s = p
   #   # CRB
   #   Z0 = meanSN(band, t, p)
   #else:
   #   # No mean LC, so just set to 0
   #   Z0 = 0
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
      if num.any(num.less(t,t1)):
         Tp = bisplev(t1, p, f, dx=1)
         T = bisplev(t1, p, f)
         eT = bisplev(t, p, ef)
         t0 = t1 - 2*T/Tp; a = T/(t1-t0)**2
         Z = num.where(num.less(t, t1), a*num.power(t-t0,2), Z)
         eZ = num.where(num.less(t, t1), eT, eZ)
         mask = mask*num.greater(Z,0)*num.greater(t, t0)
      if num.any(num.greater(t, t2)):
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

   if gen==2:
      return x0 + 0.1*(x1-x0), x1 - 0.1*(x1-x0)
   else:
      return x0+0.00001, x1-0.00001

def get_t_lim(band, param, gen):
   '''Get the limit in time-space for trusting the templates.'''
   if (band,gen) not in flux[param]:
      # Need to generate the bosplrep
      load_data(band,param,gen)
   if band in ['Y','J','H']:
      return -10.,55.
   elif band in ['u'] and gen < 3:
         return -10.,50.
   else:
      return -10.,70.
   

def dm152s(dm15):
   '''Convert from dm15 parameter to stretch.'''
   return((3.06-dm15)/2.04)

class dm15_template:
   def __init__(self):
      self.dm15 = None
      self.normalize = 1   # Do we force max of lightcurve = 0?
      self.t = None
      self.B = None;  self.eB = None
      self.V = None;  self.eV = None
      self.u = None;  self.eu = None
      self.g = None;  self.eg = None
      self.r = None;  self.er = None
      self.i = None;  self.ei = None
      self.Y = None;  self.eY = None
      self.J = None;  self.eJ = None
      self.H = None;  self.eH = None
      self.K = None;  self.eK = None
      self.H_K = None;  self.eH_K = None
      self.J_K = None;  self.eJ_K = None

   def mktemplate(self, dm15, dm15_int=None, dm15_colors='int', generate=0):
      '''2nd and 3rd arguments ignored.'''
      self.dm15 = dm15

      if dm15_colors == 'int':
         self.normalize = 1
      else:
         self.normalize = 0
      if generate:
         # generate the model light-curve from -10 to 80 in 1 day increments.
         self.t = num.arange(-15,81, 1.0)
         for band in ['B','V','u','g','r', 'i','Y','J','H','K','J_K','H_K']:
            self.__dict__[band],self.__dict__['e'+band], mask = \
               self.eval(band, self.t)
            self.__dict__['e'+band] = num.where(mask, self.__dict__['e'+band], -1.0)

   def deltaTmax(self, band, generation=1):
      '''Given the current dm15, what is the time of maximum of [band]
      relative to B-band.'''
      if band == 'V':  return 1.07
      if band == 'u':  return -1.7
      if band == 'g':  return 0.48
      if band == 'r':  return poly(self.dm15, 1.1, [1.24, -2.55, 10.1])
      if band == 'i':  return breakpoly(self.dm15, 1.51, [-2.97, 0.41, 44.44], False)
      if band == 'Y':  return breakpoly(self.dm15, 1.75, [-4.14, 0.53, 210.4], False)
      if band == 'J':  return breakpoly(self.dm15, 1.63, [-3.69, -0.12, 87.1], False)
      if band == 'H':  return breakpoly(self.dm15, 1.72, [-5.24, -1.46, 417.9], False)
      return 0

   def teval(self, band, time, dm15, gen=1, extrap=False):
      f,ef,mask = finterp(band, time, dm15, 'dm15', gen, extrap=extrap)

      return(num.where(mask,-2.5*num.log10(f),-1))

   def domain(self, band, gen=1):
      return get_t_lim(band, 'dm15', gen)

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

      # This provides a template for JHK photometry based on Kevin Krisciunas' 
      #  polynomial
      s = dm152s(self.dm15)
      if band == 'J_K':
         return(0.080 + evt/s*0.05104699 + 0.007064257*(evt/s)**2 - 0.000257906*(evt/s)**3,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 
      elif band == 'H_K':
         return(0.050 + evt/s*0.0250923 + 0.001852107*(evt/s)**2 - 0.0003557824*(evt/s)**3,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 
      elif band == 'K':
         return(0.042 + evt/s*0.02728437+ 0.003194500*(evt/s)**2 - 0.0004139377*(evt/s)**3,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 

      # Check to see if we are extrapolating beyond the training boundaries
      # If so, and we go beyond min/max dm15, use root-finding to figure out
      # what stretch would be required to map the fastest/slowest templates
      # the the requested dm15
      dmmin,dmmax = get_p_lim(band,'dm15', gen)
      if sextrap and not (dmmin < self.dm15 < dmmax):
         tmin,tmax = get_t_lim('B', 'dm15', gen)
         test_t = num.linspace(0,tmax,10)
         if self.dm15 <= dmmin:
            start = self.teval('B', 15.0, dmmin)
            target = start + (self.dm15 - dmmin)
            dmlim = dmmin
         else:
            start = self.teval('B', 15.0, dmmax)
            target = start + (self.dm15 - dmmax)
            dmlim = dmmax
         test_vals = self.teval('B', test_t, dmlim) - target
         id = num.nonzero(num.greater(test_vals,0))[0][0]
         if debug:  print("start=",start,"dm15=",self.dm15,"target = ",target)
         t0 = test_t[id-1]
         t1 = test_t[id]
         if debug:  print("t0 = ",t0,"t1 = ",t1)
         root = scipy.optimize.brentq(\
               lambda x:  self.teval('B', x, dmlim) - target, t0, t1)
         if debug:  print("root = ",root)
         s = 15./root
      else:
         s = 1.0

      # Make sure we stay within the training limits
      dm15 = self.dm15
      if sextrap and dm15 <= dmmin:
         dm15 = dmmin + 0.001
      if sextrap and dm15 >= dmmax:
         dm15 = dmmax - 0.001

      tmin,tmax = get_t_lim(band, 'dm15', gen)
      tmask = num.greater_equal(evt, tmin)*num.less_equal(evt, tmax)

      # Apply any stretch needed for extrapolating
      evt = evt*s
      evd,eevd,mask = finterp(band, evt, dm15, 'dm15', gen, extrap=extrap)

      if not extrap:
         mask = mask*tmask

      if not mag:
         return(evd, eevd, mask)
      else:
         return(-2.5*num.log10(evd), eevd/evd*1.0857, mask)


class st_template:
   def __init__(self):
      self.st = None
      self.normalize = 1   # Do we force max of lightcurve = 0?
      self.t = None
      self.B = None;  self.eB = None
      self.V = None;  self.eV = None
      self.u = None;  self.eu = None
      self.g = None;  self.eg = None
      self.r = None;  self.er = None
      self.i = None;  self.ei = None
      self.Y = None;  self.eY = None
      self.J = None;  self.eJ = None
      self.H = None;  self.eH = None
      self.K = None;  self.eK = None
      self.H_K = None;  self.eH_K = None
      self.J_K = None;  self.eJ_K = None

   def mktemplate(self, st, dm15_int=None, dm15_colors='int', generate=0):
      '''2nd and 3rd arguments ignored.'''
      self.st = st

      if generate:
         # generate the model light-curve from -10 to 80 in 1 day increments.
         self.t = num.arange(-15,81, 1.0)
         for band in ['B','V','u','g','r', 'i','Y','J','H','J_K','H_K','K']:
            self.__dict__[band],self.__dict__['e'+band], mask = \
               self.eval(band, self.t)
            self.__dict__['e'+band] = num.where(mask, self.__dict__['e'+band], -1.0)

   def deltaTmax(self, band, gen=1):
      '''Given the current dm15, what is the time of maximum of [band]
      relative to B-band.'''
      if gen > 2:
         return Tmax_splines[band](self.st)
      if band == 'V':  return 1.33
      if band == 'u':  return -1.60
      if band == 'g':  return 0.18
      if band == 'r':  return poly(self.st, 1.0, [1.56, 4.24, 25.62])
      if band == 'i':  return breakpoly(self.st, 1.00, [-3.28, 2.02, 16.69], True)
      if band == 'Y':  return breakpoly(self.st, 0.95, [-4.69, -0.08, 25.43], True)
      if band == 'J':  return breakpoly(self.st, 0.92, [-4.03, -2.42, 14.40], True)
      if band == 'H':  return breakpoly(self.st, 1.02, [-4.30, 4.26, 20.40], True)
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
         return (evt*0, evt*0, num.zeros(evt.shape, dtype=bool))

      if band == 'J_K':
         s = self.st
         return(0.080 + evt/s*0.05104699 + 0.007064257*(evt/s)**2 - 0.000257906*(evt/s)**3,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 
      elif band == 'H_K':
         s = self.st
         return(0.050 + evt/s*0.0250923 + 0.001852107*(evt/s)**2 - 0.0003557824*(evt/s)**3,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 
      elif band == 'K':
         s = self.st
         return(0.042 + evt/s*0.02728437+ 0.003194500*(evt/s)**2 - 0.0004139377*(evt/s)**3,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 

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
