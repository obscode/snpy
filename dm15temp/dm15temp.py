#!/usr/bin/env python
'''A python module (and script) for generating lightcurve templates based on
Prieto et al. 2006.  Uses Prieto's own C-code, which he made available, wrapped
up for use in Python.

August 14, 2007:  added a hack to produce NIR stretch templates valid on 
[-12,10].
Feb. 2014:  The method used up to now was really inefficient. Now we do what
             what we do with the other templates:  spline the surface.'''
import sys,os,string
import numpy as num
import scipy
#import dm15tempc
import pickle
from scipy.optimize import brentq

template_bands = ['Bs','Vs','Rs','Is']

form = "%5.1lf %10.3lf %10.6lf %10.3lf %10.6lf %10.3lf %10.6lf %10.3lf %10.6lf"
#dm15_path = '/Users/burns/CSP/k-corrections/dm15temp/'
dm15_path = os.path.dirname(globals()['__file__'])

def dm152s(dm15):
   '''Convert from dm15 parameter to stretch.'''
   return((3.06-dm15)/2.04)

class template:
   def __init__(self):
      self.dm15 = None
      self.dm15max = 1.93
      self.dm15min = 0.83
      #/* Reading data of the templates */

      self.t = None
      self.et = None
      self.B = None
      self.eB = None
      self.R = None
      self.eR = None
      self.I = None
      self.eI = None
      self.V = None
      self.eV = None
      self.J = None
      self.eJ = None
      self.H = None
      self.eH = None
      self.K = None
      self.eK = None

      f = open(os.path.join(dm15_path, "tck.pickle"))
      self.tck = pickle.load(f)
      f.close()

      self.shiftV = 0
      self.shiftR = 0
      self.shiftI = 0
   def dm152s(self, dm15):
      '''Given a value of dm15, find the stretch that will convert
      the current B-lc to that value of dm15.'''

      try:
         root = brentq(lambda x: self.tck['B'].ev(x, self.dm15)-dm15,0,30)
      except ValueError:
         return (1)

      return(15.0/root)

   def mktemplate(self, dm15, method=1, colors='none', generate=0):

      self.dm15 = dm15
      if dm15 < 0.83:  
         self.dm15=0.83
      if dm15 > 1.93:  
         self.dm15=1.93
      if dm15 < 0.83 or dm15 > 1.93:
         self.s = self.dm152s(dm15)
      else:
         self.s = 1.0
      if generate:
         self.t = arange(-10,71, 1.0)
         for band in ['B','V','R','I']:
            self.__dict__[band] = \
                  self.tck[band].ev(self.t/self.s, self.t*0+self.dm15)
            self.__dict__['e'+band] = \
                  self.tck['e_'+band].ev(self.t/self.s, self.t*0+self.dm15)
         self.eB = num.where(num.less(self.eB, 0.001), 0.001, self.eB)
         self.eV = num.where(num.less(self.eV, 0.001), 0.001, self.eV)
         self.eR = num.where(num.less(self.eR, 0.001), 0.001, self.eR)
         self.eI = num.where(num.less(self.eI, 0.001), 0.001, self.eI)
         # Here is a hack to produce a template in the same manner as Jose-Louis'
         # code.
         #gids = num.greater_equal(self.t, -12)*num.less_equal(self.t, 10)
         #t = num.compress(gids, self.t)
         self.J,self.eJ,maskJ = self.eval('J', self.t)
         self.H,self.eH,maskH = self.eval('H', self.t)
         self.K,self.eK,maskK = self.eval('K', self.t)
         self.eJ = num.where(maskJ, self.eJ, -1.0)
         self.eH = num.where(maskH, self.eH, -1.0)
         self.eK = num.where(maskK, self.eK, -1.0)

      if colors=='none':
         BmV = 0
         BmR = 0
         BmI = 0
      else:
         BmV = self.shiftV
         BmR = self.shiftR
         BmI = self.shiftI
         self.V -= BmV
         self.R -= BmR
         self.I -= BmI

   def eval(self, band, times, z=0, k=1):
      '''Evaluate, using a spline, the value of the template at specific
      times, optionally with a redshift (in the sense that the times should
      be blueshifted before interpolating.  Also returns a mask indicating
      the interpolated points (1) and the extrapolated points (0)'''
      if len(num.shape(times)) == 0:
         evt = num.array([times/(1+z)])
         scalar = 1
      else:
         evt = times/(1+z)
         scalar = 0
      if band not in self.__dict__ and band not in ['J','H','K']:
         raise AttributeError, "Sorry, band %s is not supported by dm15temp" % \
               band
      s = dm152s(self.dm15)
      if band == 'J':
         return(0.080 + evt/s*0.05104699 + 0.007064257*(evt/s)**2 - 0.000257906*(evt/s)**3,
               0.0*evt/s + 0.06, num.greater_equal(evt/s, -12)*num.less_equal(evt/s, 10)) 
      elif band == 'H':
         return(0.050 + evt/s*0.0250923 + 0.001852107*(evt/s)**2 - 0.0003557824*(evt/s)**3,
               0.0*evt/s + 0.08, num.greater_equal(evt/s, -12)*num.less_equal(evt/s, 10)) 
      elif band == 'K':
         return(0.042 + evt/s*0.02728437+ 0.003194500*(evt/s)**2 - 0.0004139377*(evt/s)**3,
               0.0*evt/s + 0.08, num.greater_equal(evt/s, -12)*num.less_equal(evt/s, 10)) 
      evd = self.tck[band].ev(evt/self.s, evt*0+self.dm15)
      eevd = self.tck['e_'+band].ev(evt/self.s, evt*0+self.dm15)
      mask = num.greater_equal(evt/self.s, -10)*num.less_equal(evt/self.s,70)

      if scalar:
         return(evd[0], eevd[0], mask[0])
      else:
         return(evd, eevd, mask)

   def MMax(self, band):
      '''Given a value of dm15, return the maximum magnitude in each filter, 
      based on Prieto et al. (2006).
      Added a bit of a hack for dm15 > 0.7:  do another linear fit, but it's going to
      be noisy, as it's based on fewer points (fix when we get more dm15's > 1.7).  Now
      the routine also returns the formal error from the fits.
      '''
      a1 = {'B':-19.319, 'V':-19.246, 'R':-19.248, 'I':-18.981, 'J':-18.57, 'H':-18.24, 'K':-18.42}
      ea = {'B':0.024, 'V':0.022, 'R':0.031, 'I':0.023, 'J':0.03, 'H':0.04, 'K':0.04}
      b1 = {'B':0.634,   'V':0.606,   'R':0.566,   'I':0.524, 'J':0, 'H':0,'K':0}
      eb = {'B':0.082, 'V':0.075, 'R':0.012, 'I':0.090,'J':0, 'H':0, 'K':0}
      a2 = {'B':-23.4506, 'V':-20.320, 'R':-20.350, 'I':-19.654}
      b2 = {'B':7.52,  'V':2.53,   'R':2.4,   'I':1.64}
      if self.dm15 < 1.7:
         a = a1; b=b1
      else:
         a = a2; b=b2
      # For those bands with not calibration yet:
      if band not in a:  return(-19.0, 0.0)
      return (a[band] + b[band]*(self.dm15 - 1.1),
              num.sqrt(ea[band]**2 + (self.dm15-1.1)**2*eb[band]**2))
  

if __name__ == "__main__":
   dm15 = float(sys.argv[1])
   t = template()
   t.mktemplate(dm15)
   print "# Constructed light curve template for dm15=%.3f mag" % (dm15)
   print "# Columns: \n#  1    t(Bmax) \n#  2-3  B-B(max)   sigma[B-B(max)]^2"
   print "#  4-5  V-V(max)   sigma[V-V(max)]^2 \n#  6-7  R-R(max)   sigma[R-R(max)]^2 \n"
   print "#  8-9  I-I(max)   sigma[I-I(max)]^2 \n"
   
   sys.stderr.write('%f %f %f\n' % (t.shiftV,t.shiftR,t.shiftI))
   for i in range(len(t.B)):
     print form % (t.t[i],t.B[i],t.eB[i],t.V[i],t.eV[i],t.R[i], t.eR[i],
           t.I[i],t.eI[i])
