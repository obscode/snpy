#!/usr/bin/env python
'''A python module to generate lightcurve templates.
Author:  Chris Burns
Version:  0.1 (extreme-alpha)

This module uses Barry Madore's GLoEs algorithm to interpolate over a surface
with heterogeneous data.  In this case, we have a surface with x=time, y=dm15
and z=flux.

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

NEW:  Using GLOES is rather slow and expensive.  However, you can always 
generate the surface and fit it with a bi-variate spline, capturing the 
information, but making things far faster.  Now, we check to see if the tck
in a file tck.pickle is available.  If so, then use this instead of calling
gloes.
'''
import sys,os,string
import numpy as num
import dm15temp2c as dm15tempc
import pickle
from scipy.interpolate import bisplev

base = os.path.dirname(globals()['__file__'])
if base == '':  base = '.'
form = "%5.1lf %10.3lf %10.6lf %10.3lf %10.6lf %10.3lf %10.6lf %10.3lf %10.6lf"
filter_numbers = {'B':0, 'V':1, 'u_s':2, 'g_s':3, 'r_s':4, 'i_s':5,
                  'Bs':6, 'Vs':7, 'u':8, 'g':9, 'r':10, 'i':11,
                  'Y':12,'J':13,'H':14}

if not os.path.isfile(os.path.join(base, 'templates.dat')):
   raise IOError, "templates.dat file not found where expected:  %s" % \
         (os.path.join(base, 'templates.dat'))
dm15tempc.load_data(base)
use_gloes = 0
use_gp = 1
sigx0 = 3.0
sigy0 = 0.3
xscale = 0.1
maxsigmax = 10.0

NIR_range = [-12,10]
try:
   f = open(os.path.join(base, 'tck.pickle'), 'r')
   btck = pickle.load(f)
   f.close()
   if os.path.isfile(os.path.join(base, 'tck2.pickle')):
      f = open(os.path.join(base, 'tck2.pickle'), 'r')
      btck2 = pickle.load(f)
      f.close()
   else:
      btck2 = None
except:
   use_gloes = 1

def dm152s(dm15):
   '''Convert from dm15 parameter to stretch.'''
   return((3.06-dm15)/2.04)

class template:
   def __init__(self):
      self.dm15 = None
      self.Rv = 2.0
      self.normalize = 1   # Do we force max of lightcurve = 0?
      self.t = None
      self.B = None;  self.eB = None
      self.V = None;  self.eV = None
      self.u_s = None;  self.eu = None
      self.g_s= None;  self.eg = None
      self.r_s = None;  self.er = None
      self.i_s = None;  self.ei = None
      self.Bs = None;  self.eBs = None
      self.Vs = None;  self.eVs = None
      self.u = None;  self.eu = None
      self.g = None;  self.eg = None
      self.r = None;  self.er = None
      self.i = None;  self.ei = None
      self.Y = None;  self.eY = None
      self.J = None;  self.eJ = None
      self.H = None;  self.eH = None
      self.K = None;  self.eK = None

      # for when we have more than one
      # for when we have more than one
      self.calibration = 'CRB'


   def mktemplate(self, dm15, dm15_int=None, dm15_colors='int', generate=0):
      '''2nd and 3rd arguments ignored.'''
      self.dm15 = dm15
      #if not use_gloes:
      #   if self.dm15 < 0.7:  self.dm15 = 0.7
      #   if self.dm15 > 2.0:  self.dm15 = 2.0

      if dm15_colors == 'int':
         self.normalize = 1
      else:
         self.normalize = 0
      if generate:
         # generate the model light-curve from -10 to 80 in 1 day increments.
         self.t = num.arange(-15,81, 1.0)
         for band in ['B','V','u_s','g_s','r_s','i_s','Bs','Vs','u','g','r',
                      'i','Y','J','H','K']:
            self.__dict__[band],self.__dict__['e'+band], mask = \
               self.eval(band, self.t)
            self.__dict__['e'+band] = num.where(mask, self.__dict__['e'+band], -1.0)

   def eval(self, band, times, z=0, mag=1):
      if len(num.shape(times)) == 0:
         evt = num.array([times/(1+z)])
         scalar = 1
      else:
         evt = times/(1+z)
         scalar = 0
      if band not in filter_numbers and band != 'K':
         raise AttributeError, "Sorry, band %s is not supported by dm15temp2" % \
               band

      # This provides a template for JHK photometry based on Kevin Krisciunas' polynomial
      s = dm152s(self.dm15)
      #if band == 'J':
      #   return(0.080 + evt/s*0.05104699 + 0.007064257*(evt/s)**2 - 0.000257906*(evt/s)**3,
      #         num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
      #         num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 
      #elif band == 'H':
      #   return(0.050 + evt/s*0.0250923 + 0.001852107*(evt/s)**2 - 0.0003557824*(evt/s)**3,
      #         num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
      #         num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 
      if band == 'K':
         return(0.042 + evt/s*0.02728437+ 0.003194500*(evt/s)**2 - 0.0004139377*(evt/s)**3,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])*0.08,
               num.greater_equal(evt/s, NIR_range[0])*num.less_equal(evt/s, NIR_range[1])) 
      if use_gloes or not (0.6 < self.dm15 < 2.0) or band in ['Y','J','H']:
         evd = evt*0
         eevd = evt*0
         dm15tempc.dm15temp(filter_numbers[band], self.dm15, evt, evd, eevd, 
               self.normalize, sigx0, xscale, maxsigmax, sigy0)
      else:
         if not isinstance(evt, num.ndarray):
            evd = num.array([bisplev(evt, self.dm15, btck[band])])
            eevd = num.array([bisplev(evt, self.dm15, btck["e_"+band])])
            if btck2 is not None and use_gp:
               evd += num.array([bisplev(evt, self.dm15, btck2[band])])
               eevd += num.array([bisplev(evt, self.dm15, btck2["e_"+band])])

         else:
            if len(evt) == 1:
               evd = num.array([bisplev(evt, self.dm15, btck[band])])
               eevd = num.array([bisplev(evt, self.dm15, btck["e_"+band])])
               if btck2 is not None and use_gp:
                  evd += num.array([bisplev(evt, self.dm15, btck2[band])])
                  eevd += num.array([bisplev(evt, self.dm15, btck2["e_"+band])])
            else:
               evd = bisplev(evt, self.dm15, btck[band])[:,0]
               eevd = bisplev(evt, self.dm15, btck["e_"+band])[:,0]
               if btck2 is not None and use_gp:
                  evd += bisplev(evt, self.dm15, btck[band])[:,0]
                  eevd += bisplev(evt, self.dm15, btck["e_"+band])[:,0]
      # Now mask out [-10,70]
      gids = num.greater_equal(evt, -10)*num.less_equal(evt,70)

      # Here is some extra masking for high dm15s where the bisplev goes wonky
      if band[0] == 'i':
         if self.dm15 > 2:  
            gids = gids*num.greater(evt, -4.0)
         elif 1.85 < self.dm15 <= 2.0:  
            gids = gids*num.greater(evt, -5.0)
         elif 1.80 < self.dm15 <= 1.85:
            gids = gids*num.greater(evt, -6.0)
         elif 1.70 < self.dm15 <= 1.80:
            gids = gids*num.greater(evt, -7.0)
      if band[0] == 'r':
         if self.dm15 > 2:
            gids = gids*num.greater(evt, -5.0)
         elif 1.9 < self.dm15 <=2.0:
            gids = gids*num.greater(evt, -6.0)
         elif 1.83 < self.dm15 <=1.9:
            gids = gids*num.greater(evt, -7.0)
         elif 1.76 < self.dm15 <=1.83:
            gids = gids*num.greater(evt, -8.0)
      if band[0] == 'g':
         if self.dm15 > 2:
            gids = gids*num.greater(evt, -6.0)
         elif 1.88 < self.dm15 <=2.0:
            gids = gids*num.greater(evt, -7.0)
         elif 1.83 < self.dm15 <=1.88:
            gids = gids*num.greater(evt, -8.0)
      if band[0] in ['Y','J','H']:
         gids = gids*num.greater(evt, -8)*num.less(evt, 57)
         if self.dm15 > 1.7:  gids = gids*num.less(evt, 40)

      if not mag:
         return(evd, eevd, gids)
      zids = num.greater(evd, 0)
      flux = num.where(zids, evd, 1.0)
      evd = -2.5*num.log10(flux)
      eevd = num.where(zids, eevd/flux*1.0857, 9.99)
      gids = gids*zids
      if scalar:
         return(evd[0], eevd[0], gids[0])
      else:
         return(evd, eevd, gids)

   def MMax(self, band):

      # These have been solved as a function of Rv
      As = {'Bs':[-19.148, -0.065, 0.002],
            'Vs':[-19.146, -0.067, 0.002],
            'u': [-18.722, -0.068, 0.003],
            'g': [-19.187, -0.067, 0.002],
            'r': [-19.061, -0.065, 0.002],
            'i': [-18.486, -0.057, 0.001]}
      eAs = {'Bs':0.03, 'Vs':0.03, 'g':0.03, 'r':0.03, 'i':0.03, 'u':0.04}
      eBs = {'Bs':0.08, 'Vs':0.08, 'g':0.09, 'r':0.08, 'i':0.08, 'u':0.12}
      Bs = {'Bs':[0.561, 0.055, 0.001],
            'Vs':[0.421, 0.070, -0.002],
            'u': [0.891,  0.073, -0.001],
            'g': [0.487,  0.061, -0.000],
            'r': [0.246,  0.080, -0.004],
            'i': [0.044,  0.075, -0.004]}

      if band not in ['Bs','Vs','u','g','r','i']:
         return(-19.0, 0.02)
      else:
         A = As[band][0] + As[band][1]*self.Rv + As[band][2]*self.Rv**2
         B = Bs[band][0] + Bs[band][1]*self.Rv + Bs[band][2]*self.Rv**2
         return(A + B*(self.dm15 - 1.1), 
               num.sqrt(eAs[band]**2 + eBs[band]**2*(self.dm15-1.1)**2))



if __name__ == "__main__":
   use_gloes = 1
   if len(sys.argv) < 3:
      print "Usage:  dm15temp2 dm15 filter"
      print "        dm15 = decline rate"
      print "        filter = one of u,g,r,i,Bs,Vs, u_s,g_s,r_s,i_s,B,V"
      sys.exit(1)

   dm15 = float(sys.argv[1])
   filter = sys.argv[2]
   t = template()
   t.mktemplate(dm15, generate=0, dm15_colors='int')
   print "# This is dm15temp2 v. 0.1, using GLoEs with built-in weights"
   print "# Constructed light curve template for dm15=%.3f mag, filter=%s" % (dm15,filter)
   print "# Columns:  \n#  1    t(Bmax) \n#  2-3  M-M(max)   sigma[M-M(max)]"
   ts = num.arange(-10, 70, 1.0)
   m,em,mask = t.eval(filter, ts)
   for i in range(len(m)):
      if mask[i] == 1:
         print ts[i], m[i], em[i]
