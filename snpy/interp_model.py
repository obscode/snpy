'''inter_model.py:  a class that defines a SN model to be fit by SNOOPY.  Unlike the
other module (model.py), this model is based on a interpolate1d.py and so doesn't have
a least-squares fit routine (at least not at the user level).  This supplants the
template-fitting with splines.'''

from __future__ import print_function
import os
from snpy import kcorr
from snpy.filters import fset
from numpy import *
from scipy import stats
from scipy.optimize import brentq
from numpy import median, diag, atleast_1d, std
from numpy.linalg import inv

#from snpy.utils import fit1dcurve
#from snpy.utils import InteractiveFit
import fit1dcurve
import InteractiveFit

if 'gp' in fit1dcurve.functions:
   default_type = 'gp'
else:
   default_type = 'hyperspline'

debug = 0

class model:


   def __init__(self, parent):
      self.parameters = {}       # needed for compatibility.
      self.errors = {}
      self.parent = parent       # make a link to the SN object
      self.parent.model = self   # make a link to the model
      self.do_kcorr = 1       # Do we perform initial k-corrections?
      self.interps = {}       # Keep an Interpolator object here

      self.model_in_mags = {}

   def __getattr__(self, key):
      if key == '_fbands':
         return list(self.__dict__['interps'].keys())
      raise AttributeError("Model object has no attribute %s" % (key))

   def fit_lc(self, band, type=None, fit_flux=False, **args):
      '''Fit a single filter with an interpolator of type [type], passing
      [args] to its __init__ function.  If [fit_flux] is True, fit the 
      data in flux space (warning:  you could end up with negative
      fluxes in this case).'''

      # check type for consistency
      if type is None:
         type = default_type
      if type not in fit1dcurve.functions:
         print("Error:  specified type is not supported.  Try one of:")
         fit1dcurve.list_types()
         raise ValueError

      x = self.parent.data[band].MJD
      if fit_flux:
         y = self.parent.data[band].flux
         ey = self.parent.data[band].e_flux
         self.model_in_mags[band] = False
      else:
         y = self.parent.data[band].mag
         ey = self.parent.data[band].e_mag
         self.model_in_mags[band] = True
      self.interps[band] = fit1dcurve.Interpolator(type, x, y, ey, **args)


   def fit(self, bands, type=None, fit_flux=False, **args):
      '''Globally fit the [bands] data with an interpolator of type [type]
      using the same arguments [args] for each.  Simply calls self.fit_lc()
      repeatedly with the same arguments.
      '''
      for b in bands:
         self.fit_lc(b, type, fit_flux, **args)

   def __call__(self, band, t):
      '''Call the interpolator for band [band] at time [t].  Returns the tuple
      (y,ey,m) where y is the interpolated point, ey is the error, and m is
      a validity mask (data is good where True)'''
      scalar = (len(shape(t)) == 0)
      t = atleast_1d(t)
      f,mask = self.interps[band](t)
      ef = self.interps[band].error(t)

      if not self.model_in_mags[band]:
         bids = less_qual(ev, 0)
         ev[bids] = 1.0
         mask = mask*bids
         ef = 1.0857*ef/f
         f = -2.5*log10(ev) + self.parent.data[band].filter.zp

      return f, ef, mask

   def compute_lc_params(self, N=50):
      '''Compute dm15, Tmax, Mmax, and covariances for each band that has an 
      interpolator setup.'''

      # rest-frame of the SN
      day15 = 15*(1 + self.parent.z)
      for band in self.interps:
         zp = self.parent.data[band].filter.zp
         Tmaxs = []
         Mmaxs = []
         dm15s = []
         inter = self.interps[band]
         for i in range(N):
            # Find Tmax
            if i:  
               inter.draw()
               if inter.deriv(Tmaxs[0]-1)*inter.deriv(Tmaxs[0]+1) > 0:
                  Tmaxs.append(-1)
                  Mmaxs.append(-1)
                  dm15s.append(-1)
                  continue
               else:
                  # self.deriv() for GP is really slow, so try this faster 
                  #    approach
                  xs = array([brentq(inter.deriv, Tmaxs[0]-1, Tmaxs[0]+1)])
                  ys = array(inter(xs[0]))
                  if self.model_in_mags[band]:
                     curvs = array([1])
                  else:
                     curvs = array([-1])
                  #xs,ys,curvs = inter.find_extrema(xmin=Tmaxs[0]-10, xmax=Tmaxs[0]+10)
            else:
               xs,ys,curvs = inter.find_extrema()
            if self.model_in_mags[band]:
               xs = xs[greater(curvs,0)]
               ys = ys[greater(curvs,0)]
            else:
               xs = xs[less(curvs,0)]
               ys = ys[less(curvs,0)]
            if len(xs) == 0:
               if i == 0:
                  # If we can't find maximum on the data, we're done
                  for attrib in ['Tmax','Mmax','dm15','e_Tmax','e_Mmax',
                        'e_dm15','cov_Tmax_dm15','cov_Tmax_Mmax',
                        'cov_Mmax_dm15']:
                     self.parent.data[band].__dic__[attrib] = None
                  return
               Tmaxs.append(-1)
               Mmaxs.append(-1)
               dm15s.append(-1)
               continue
            Tmaxs.append(xs[0])

            # Find Mmax/dm15
            y15,m15 = inter(Tmaxs[-1] + day15)
            if self.model_in_mags[band]:
               Mmaxs.append(ys[0])
               if m15:
                  dm15s.append(y15 - Mmaxs[-1])
               else:
                  dm15s.append(-1)
            else:
               if ys[0] < 0:
                  Mmaxs.append(-1)
                  dm15s.append(-1)
               else:
                  Mmaxs.append(-2.5*log10(ys[0]) + zp)
                  if m15:
                     dm15s.append(-2.5*log10(y15) + zp - Mmaxs[-1])
                  else:
                     dm15s.append(-1)
         #Now compute stats
         inter.reset_mean()
         Tmaxs,Mmaxs,dm15s = list(map(array, [Tmaxs, Mmaxs, dm15s]))

         Tgids = greater(Tmaxs, 0)
         if sum(Tgids)*1.0/len(Tgids) < 0.8:
            print("Warning!  More than 20% of MC realizations had bad Tmax")
         Mgids = greater(Mmaxs, 0)
         if sum(Mgids)*1.0/len(Mgids) < 0.8:
            print("Warning!  More than 20% of MC realizations had bad Mmax")
         dgids = greater(dm15s, 0)
         if sum(dgids)*1.0/len(dgids) < 0.8:
            print("Warning!  More than 20% of MC realizations had bad dm15")
         self.parent.data[band].Tmax = Tmaxs[0]
         self.parent.data[band].e_Tmax = std(Tmaxs[Tgids]) 
         self.parent.data[band].Mmax = Mmaxs[0]
         self.parent.data[band].e_Mmax = std(Mmaxs[Mgids]) 
         self.parent.data[band].dm15 = dm15s[0]
         self.parent.data[band].e_dm15 = std(dm15s[dgids]) 
         self.parent.data[band].cov_Tmax_dm15 = \
               sum(compress(Tgids*dgids, (Tmaxs-Tmaxs[0])*(dm15s-dm15s[0])))/\
               (sum(Tgids*dgids) - 1)
         self.parent.data[band].cov_Tmax_Mmax = \
               sum(compress(Tgids*Mgids, (Tmaxs-Tmaxs[0])*(Mmaxs-Mmaxs[0])))/\
               (sum(Tgids*Mgids) - 1)
         self.parent.data[band].cov_Mmax_dm15 = \
               sum(compress(Mgids*dgids, (Mmaxs-Mmaxs[0])*(dm15s-dm15s[0])))/\
               (sum(Mgids*dgids) - 1)
         return

