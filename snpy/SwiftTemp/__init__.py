#!/usr/bin/env python
'''A python module to generate lightcurve templates for the Swift
UV filters UVM2, UVW1, and UVW2
Author:  Chris Burns
Version:  0.1

This is a pure stretch template based on SN2011fe.
'''
import sys,os,string
import numpy as num
from scipy.interpolate import interp1d
import scipy.optimize
import pickle

debug=0

template_bands = ['UVW1','UVW2','UVM2']

base = os.path.dirname(globals()['__file__'])
if base == '':  base = '.'

def dm152s(dm15):
   '''Convert from dm15 parameter to stretch.'''
   return((3.06-dm15)/2.04)

class dm15_template:
   def __init__(self):
      self.dm15 = None
      self.normalize = 1   # Do we force max of lightcurve = 0?
      self.t = None
      self.interp = {}

      self.setup_interpolators()

      self.UVM2 = None;  self.eUVM2 = None
      self.UVW1 = None;  self.eUVW1 = None
      self.UVW2 = None;  self.eUVW2 = None

   def __getstate__(self):

      d = self.__dict__.copy()
      d['interp'] = {}
      return d

   def setup_interpolators(self):

      x,y = num.loadtxt(os.path.join(base, "SNIa_w1_template_fe.dat"),
            unpack=True)
      self.interp['UVW1'] = interp1d(x,y, bounds_error=False)
      x,y = num.loadtxt(os.path.join(base, "SNIa_w2_template_fe.dat"),
            unpack=True)
      self.interp['UVW2'] = interp1d(x,y, bounds_error=False)
      x,y = num.loadtxt(os.path.join(base, "SNIa_m2_template_fe.dat"),
            unpack=True)
      self.interp['UVM2'] = interp1d(x,y, bounds_error=False)

   def mktemplate(self, dm15, dm15_int=None, dm15_colors='int', generate=0):
      '''2nd and 3rd arguments ignored.'''
      self.dm15 = dm15

      if generate:
         # generate the model light-curve from -10 to 80 in 1 day increments.
         self.t = num.arange(-15,81, 1.0)
         for band in ['UVW1','UVW2','UVM2']:
            self.__dict__[band],self.__dict__['e'+band], mask = \
               self.eval(band, self.t)
            self.__dict__['e'+band] = num.where(mask, self.__dict__['e'+band], -1.0)

   def deltaTmax(self, band):
      '''Given the current dm15, what is the time of maximum of [band]
      relative to B-band.'''
      if band == 'UVW1':  return -2.3
      if band == 'UVW2':  return -1.7
      if band == 'UVM2':  return -0.2
      return 0

   def domain(self, band):
      '''returns the valid domain of the template'''
      if band not in self.interp:
         self.set_interpolators()
      s = dm152s(self.dm15)
      return (self.interp[band].x.min()*s, self.interp[band].x.max()*s)

   def eval(self, band, times, z=0, mag=1, sextrap=1, gen=1, toff=True):
      '''Evaluate the template in band [band] at epochs [times].  Optionally
      redshift by (1+[z]).  If [mag]=1, return in magnitudes, otherwise return
      in flux units.  If [sextrap]=1, extrapolate beyond the training sample
      by using a stretch.  Use [gen] to specifiy the generation of the template.
      If you want the Tmax - Tmax(B) offset applied, set [toff] to True,
      otherwise, Tmax will be at 0 for every filter.'''

      if band not in self.interp:
         self.setup_interpolators()

      if toff:
         evt = (times - self.deltaTmax(band))/(1+z)
      else:
         evt = times/(1+z)

      # Apply a stretch consistent with fits to dm15
      s = dm152s(self.dm15)
      evm = self.interp[band](evt/s)
      eevm = num.zeros(evm.shape)   # No uncertainties
      mask = ~num.isnan(evm)
      evm[~mask] = -1
      if mag:
         return evm,eevm,mask
      else:
         return num.power(10, -0.4*evm),eevm,mask

class st_template:
   def __init__(self):
      self.st = None
      self.normalize = 1   # Do we force max of lightcurve = 0?
      self.t = None
      self.interp = {}

      self.setup_interpolators()

      self.UVM2 = None;  self.eUVM2 = None
      self.UVW1 = None;  self.eUVW1 = None
      self.UVW2 = None;  self.eUVW2 = None

   def __getstate__(self):

      d = self.__dict__.copy()
      d['interp'] = {}
      return d

   def setup_interpolators(self):

      x,y = num.loadtxt(os.path.join(base, "SNIa_w1_template_fe.dat"),
            unpack=True)
      self.interp['UVW1'] = interp1d(x,y, bounds_error=False)
      x,y = num.loadtxt(os.path.join(base, "SNIa_w2_template_fe.dat"),
            unpack=True)
      self.interp['UVW2'] = interp1d(x,y, bounds_error=False)
      x,y = num.loadtxt(os.path.join(base, "SNIa_m2_template_fe.dat"),
            unpack=True)
      self.interp['UVM2'] = interp1d(x,y, bounds_error=False)

   def mktemplate(self, st, dm15_int=None, dm15_colors='int', generate=0):
      '''2nd and 3rd arguments ignored.'''
      self.st = st

      if generate:
         # generate the model light-curve from -10 to 80 in 1 day increments.
         self.t = num.arange(-15,81, 1.0)
         for band in ['UVW1','UVW2','UVM2']:
            self.__dict__[band],self.__dict__['e'+band], mask = \
               self.eval(band, self.t)
            self.__dict__['e'+band] = num.where(mask, self.__dict__['e'+band], -1.0)

   def deltaTmax(self, band):
      '''Given the current dm15, what is the time of maximum of [band]
      relative to B-band.'''
      if band == 'UVW1':  return -2.3
      if band == 'UVW2':  return -1.7
      if band == 'UVM2':  return -0.2
      return 0

   def domain(self, band):
      '''returns the valid domain of the template'''
      s = dm152s(self.dm15)
      if band not in self.interp:
         self.set_interpolators()
      return (self.interp[band].x.min()*s, self.interp[band].x.max()*s)

   def eval(self, band, times, z=0, mag=1, sextrap=1, gen=1, toff=True):
      '''Evaluate the template in band [band] at epochs [times].  Optionally
      redshift by (1+[z]).  If [mag]=1, return in magnitudes, otherwise return
      in flux units.  If [sextrap]=1, extrapolate beyond the training sample
      by using a stretch.  Use [gen] to specifiy the generation of the template.
      If you want the Tmax - Tmax(B) offset applied, set [toff] to True,
      otherwise, Tmax will be at 0 for every filter.'''

      if band not in self.interp:
         self.setup_interpolators()

      if toff:
         evt = (times - self.deltaTmax(band))/(1+z)
      else:
         evt = times/(1+z)

      # Apply a stretch consistent with fits to dm15
      evm = self.interp[band](evt/self.st)
      eevm = num.zeros(evm.shape)   # No uncertainties
      mask = ~num.isnan(evm)
      evm[~mask] = -1.0
      if mag:
         return evm,eevm,mask
      else:
         return num.power(10, -0.4*evm),eevm,mask
