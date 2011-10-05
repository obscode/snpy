## Automatically adapted for numpy.oldnumeric Feb 04, 2009 by ipython

'''Create a dm15 template from the Hsiao dm15 sequene I created.'''

import FITS
from numpy.oldnumeric import *
from snpy.filters import fset
import os
import pickle
from scipy.interpolate import bisplev

base = os.path.dirname(os.path.abspath(globals()['__file__']))
base = os.path.join(os.path.dirname(base), 'typeIa')
ff = FITS.FITS(os.path.join(base,'Hsiao_sequence_0.2.fits.gz'))
ff2 = FITS.FITS(os.path.join(base,'91bg_sequence_0.2.fits.gz'))

try:
   f = open(os.path.join(base, 'tck_destiny.pickle'), 'r')
   btck = pickle.load(f)
   f.close()
   use_fits = 0
except:
   use_fits = 1

spec_data = ff.data()
dm15s = ff['CRVAL3'] + ff['CDELT3']*arange(ff['NAXIS3'])
days = ff['CRVAL2'] + ff['CDELT2']*arange(ff['NAXIS2'])
waves = ff['CRVAL1'] + ff['CDELT1']*arange(ff['NAXIS1'])

spec_data2 = ff2.data()
dm15s2 = ff2['CRVAL3'] + ff2['CDELT3']*arange(ff2['NAXIS3'])
days2 = ff2['CRVAL2'] + ff2['CDELT2']*arange(ff2['NAXIS2'])
waves2 = ff2['CRVAL1'] + ff2['CDELT1']*arange(ff2['NAXIS1'])

def interp_spec2(day, dm15):
   i0 = (dm15 - ff2['CRVAL3'])/ff2['CDELT3']
   j0 = (day - ff2['CRVAL2'])/ff2['CDELT2']

   i = int(i0)
   j = int(j0)
   if i < 0:  i == 0
   if j < 0:  j == 0
   if i+1 >= len(dm15s2):  i = len(dm15s2) - 2
   if j+1 >= len(days2):  j = len(days2) - 2
   ff00 = spec_data2[i,j]
   ff10 = spec_data2[i+1,j]
   ff01 = spec_data2[i, j+1]
   ff11 = spec_data2[i+1, j+1]
   day0 = days2[j];  day1 = days2[j+1]
   dm0 = dm15s2[i];  dm1 = dm15s2[i+1]
   w11 = (dm15 - dm0)*(day - day0)
   w00 = (dm1 - dm15)*(day1 - day)
   w01 = (dm1 - dm15)*(day - day0)
   w10 = (dm15 - dm0)*(day1 - day)
   tw = (w00 + w01 + w10 + w11)

   data = (ff00*w00 + ff01*w01 + ff10*w10 + ff11*w11)/tw

   return(data)

def interp_spec(day, dm15):
   i0 = (dm15 - ff['CRVAL3'])/ff['CDELT3']
   j0 = (day - ff['CRVAL2'])/ff['CDELT2']

   i = int(i0)
   j = int(j0)
   if i+1 >= len(dm15s):  i = len(dm15s) - 2
   if j+1 >= len(days):  j = len(days) - 2
   ff00 = spec_data[i,j]
   ff10 = spec_data[i+1,j]
   ff01 = spec_data[i, j+1]
   ff11 = spec_data[i+1, j+1]
   day0 = days[j];  day1 = days[j+1]
   dm0 = dm15s[i];  dm1 = dm15s[i+1]
   w11 = (dm15 - dm0)*(day - day0)
   w00 = (dm1 - dm15)*(day1 - day)
   w01 = (dm1 - dm15)*(day - day0)
   w10 = (dm15 - dm0)*(day1 - day)
   tw = (w00 + w01 + w10 + w11)

   data = (ff00*w00 + ff01*w01 + ff10*w10 + ff11*w11)/tw

   return(data)

def make_lc(band, dm15, day_evals):
   if dm15 > 1.7:
      specs = [interp_spec2(day, dm15) for day in day_evals]
      mags = [fset[band].synth_mag(waves2, specs[i]) for i in range(len(specs))]
   else:
      specs = [interp_spec(day, dm15) for day in day_evals]
      mags = [fset[band].synth_mag(waves, specs[i]) for i in range(len(specs))]
   mags = array(mags)
   return(mags)


class template:

   def __init__(self):
      self.dm15 = 1.1
      self.dm15max = dm15s[-1]
      self.dm15min = dm15s[0]

      self.t = days
      for i in range(15):
         self.__dict__['d%d'%i] = None
         self.__dict__['ed%d'%i] = None
      for i in range(10):
         self.__dict__['dw%d'%i] = None
         self.__dict__['edw%d'%i] = None
   
   def mktemplate(self, dm15, dm15_int=None, dm15_colors='int', generate=0):
      self.dm15 = dm15
      if generate:
         for i in range(15):
            self.__dict__['d%d'%i],self.__dict__['ed%d'%i],mask= \
                  self.eval('d%d' % i, self.t)
            self.__dict__['ed%d'%i] = where(mask, self.__dict__['ed%d'%i], -1.0)
         for i in range(10):
            self.__dict__['dw%d'%i],self.__dict__['edw%d'%i],mask= \
                  self.eval('dw%d' % i, self.t)
            self.__dict__['edw%d'%i] = where(mask, self.__dict__['edw%d'%i], -1.0)


   def MMax(self, band):
      if self.dm15 < 1.7:
         return (-19.319 + 0.634*(self.dm15 - 1.1), 
                 sqrt(0.024**2 + (self.dm15-1.1)**2*0.082**2))
      else:
         return (-23.4506 + 7.52*(self.dm15 - 1.1),
                 sqrt(0.024**2 + (self.dm15-1.1)**2*0.082**2))

   def eval(self, band, times, z=0, mag=1):
      if len(shape(times)) == 0:
         evt = array([times/(1+z)])
         scalar = 1
      else:
         evt = times/(1+z)
         scalar = 0

      if band not in self.__dict__:
         raise AttributeError, "Sorry, band %s is not supported by this module" %\
               band

      if self.dm15 > 1.7:
         gids = greater_equal(evt, days2[0])*less_equal(evt, days2[-1])
      else:
         gids = greater_equal(evt, days[0])*less_equal(evt, days[-1])
      eevd = gids*0.01
      if use_fits:
         evd = make_lc(band, self.dm15, evt)
      else:
         gids = gids*greater_equal(evt, -10)*less_equal(evt, 70)
         if isinstance(evt, not ArrayType):
            evd = array([bisplev(evt, self.dm15, btck[band])])
         else:
            if len(evt) == 1:
               evd = array([bisplev(evt, self.dm15, btck[band])])
            else:
               evd = bisplev(evt, self.dm15, btck[band])[:,0]
      
      if scalar:
         return(evd[0], eevd[0], gids[0])
      else:
         return(evd,eevd,gids)   

