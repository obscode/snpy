from numpy import *
from filters import fset
import numpy.random as RA
from utils import stats
from utils import fit1dcurve
import plotmod
from scipy.optimize import brentq

if 'gp' in fit1dcurve.functions:
   default_method = 'gp'
else:
   default_method = 'hyperspline'

class lc:
   '''This class contains the data representing a lightcurve.  This means 
   the time (in MJD), magnitudes, uncertainties, fluxes, etc.  It is meant 
   to be contained by the super class, to make for somewhat more convenient 
   access to the data.  The member data and functions are as follows:
      - self.parent         The parent super class
      - self.band           Which band does the data represent
      - self.MJD            The time data in Modified Julian Day (JD - 2400000.5)
      - self.t              (self.MJD - parent.Tmax)  Quick way to get epoch
      - self.mag            magnitudes
      - self.e_mag          uncertainties therein
      - self.covar          full error matrix, if known (otherwise, diagonal based on
                            self.e_mag
      - self.flux           flux-based brightnesses
      - self.e_flux         uncertainties therein
      - self.mask           boolean array, used to mask out bad/
                            inappropriate values
      - self.filter         reference to the filter instance'''

   def __init__(self, parent, band, MJD, mag, e_mag, restband=None, K=None):
      self.parent = parent     # pointer to the SN containing class
      self.band = band         # observed filter
      self.MJD = MJD           # time (in Modified Julian Day... or whatever)
      self.magnitude = mag     # magnitude
      self.e_mag = e_mag       # error in magnitude
      if sometrue(equal(self.e_mag,0)):
         print "Warning:  you have errors that are zero:  setting them to 0.001 mag."
         print "If you don't like that, fix your errors."
         self.e_mag = where(equal(self.e_mag, 0), 0.001, self.e_mag)
      self.covar = None        # full error matrix (if errors are correlated)
      self.restband = restband # the observed band in the frame of the SN
      self.K = K               # k-correction to convert mag to restband
      self.debug=0             # Of course, we'll never need this!

      # a way to mask out bad data from fitting
      self.mask = ones(self.mag.shape, dtype=bool) 
      self.filter = fset[self.band] # instance of filter object for this band

      #self.tck = None          # Spline solution
      #self.model = None        # a model
      #self.model_t = None      #  time for the model
      #self.model_ks = None     #  k-corrections for the model
      #self.model_type = None   #  How the model was constructed:  'template',
      #                         #  'template_colors', or 'spline'
      self.interp = None
      self.model_flux = 0       # is the model in flux-space?
      self.tmin = None         # range over which the model is valid
      self.tmax = None

      self.Tmax = None
      self.Mmax = None
      self.e_Tmax = None
      self.e_Mmax = None
      self.cov_Tmax_Mmax = None
      self.cov_Tmax_dm15 = None
      self.cov_Mmax_dm15 = None

      self.mp = None     # reference to multiplot or PanelPlot

   def get_t(self):
      try:
         t = self.MJD - self.parent.Tmax
      except:
         t = self.MJD
      return(t)

   def get_flux(self):
      return(power(10.0, -0.4*(self.mag - self.filter.zp)))

   def get_e_flux(self):
      return(self.get_flux()*self.e_mag/1.0857)

   def get_covar(self, flux=1):
      '''returns the error matrix in flux units (unless flux=0).  If this was
      not setup by the user, the variance will be returned as a 1-D array.'''
      if self.covar is not None:
         if flux:
            f = self.get_flux()
            return self.covar*f[NewAxis,:]*f[:,NewAxis]/1.17874
         else:
            return self.covar
      else:
         if flux:
            return power(self.get_e_flux(),2)
         else:
            return power(self.e_mag,2)

   def __getattr__(self, name):
      if name == "flux":
         return self.get_flux()
      elif name == "e_flux":
         return self.get_e_flux()
      elif name == "t":
         return self.get_t()
      elif name == 'mag':
         if self.K is None:
            return(self.magnitude)
         else:
            try:
               return(self.magnitude - self.K)
            except:
               raise AttributeError, "Error:  k-corrections and magnitudes incompatible"
      elif 'interp' in self.__dict__ and self.__dict__['interp'] is not None:
         if name in self.__dict__['interp'].pars:
            return getattr(self.__dict__['interp'],name)
      raise AttributeError, "Error:  attribute %s not defined" % (name)

   def __setattr__(self, name, value):
      if 'interp' in self.__dict__:
         if self.__dict__['interp'] is not None:
            if name in self.__dict__['interp'].pars:
               setattr(self.__dict__['interp'], name, value)
               self.replot()
               return
      if name == 'mask':
         if 'mask' in self.__dict__.keys():
            raise TypeError, "lc instance's mask attribute must be modifed in-place"
      self.__dict__[name] = value

   def __getstate__(self):
      # Need this because MPL objects are not pickleable.
      odict = self.__dict__.copy()
      if 'mp' in odict:  del odict['mp']
      if 'pts' in odict:  del odict['pts']
      return odict

   def time_sort(self):
      ids = argsort(self.MJD)
      self.MJD = take(self.MJD, ids)
      self.magnitude = take(self.magnitude, ids)
      self.e_mag = take(self.e_mag, ids)
      if self.K is not None:
         self.K = take(self.K, ids)

   def mask_epoch(self, tmin, tmax):
      '''Update the lc's mask to only include data between tmin and tmax.'''
      self.mask *= greater_equal(self.t, tmin)
      self.mask *= less_equal(self.t, tmax)

   def mask_emag(self, max):
      '''Update the lc's mask to only include data with e_mag < max.'''
      self.mask *= less_equal(self.e_mag, max)

   def eval(self, times, t_tol=-1, epoch=0):
      '''Interpolate (if required) the data to time 'times'.  If there is a data point 
      closer than t_tol away from a requested time, that value is used without 
      interpolation.  If epoch is nonzero, then times are interpreted relative
      to parent.Tmax'''
      if self.interp is None:
         raise AttributeError, "Error.  To interpolate, you need to fit a template first."

      times = atleast_1d(times)
      if epoch: times = times + self.parent.Tmax

      evm,mask = self.interp(times)
      if self.model_flux:
         mask = mask*greater(evm,0)
         evm = where(evm <=0, 1.0, evm)
         evm = -2.5*log10(evm) + self.filter.zp

      if t_tol > 0:
         # Now, we scan t and eval_t and find where they are less than tol.  
         # In these cases, we take the average of any matching times
         delta = absolute(self.MJD[NewAxis,:] - times[:,NewAxis])
         cond = less(delta, t_tol)
         values = array([self.mag]*len(times))*cond
         N = sum(cond, axis=1)
         w = where(N > 0, N, 1)
         s = sum(values, axis=1)
         evm = where(N > 0, s/w, evm)

      return(evm,mask)

   def spline_fit(self, fitflux=0, do_sigma=1, Nboot=100, keep_boot=1, method='spline2',
         **args):
      '''Make a spline template of the lightcurve.  If fitflux=1, then fit in 
      flux-space.  If do-sigma=1, do a monte-carlo simulation to get an idea of
      the uncertainties in the final parameters.  Nboot controls the number of
      bootstrap interations.  You can choose to use a simple
      spline (default) or "hyper-spline" using the method keyword.   The rest
      of the arguments are passed to self.mkspline (if method='spline') or
      self.mkspline2 (if method= 'spline2'.  Upon successful completion of the
      routine, the following member variables will be populated: 
         Tmax, e_Tmax: time of maximum (and error if do_sigma=1) 
         Mmax, e_Mmax:   peak magnitude
         dm15, e_dm15:   delta-m 15
         model_type:     "spline" or "spline2"
         tck:            spline info '''
      return self.template(fitflux, do_sigma, Nboot, keep_boot, method, **args)

   def list_types(self):
      print "the following methods are available for constructing templates:"
      fit1dcurve.list_types()

   def template(self, fitflux=False, do_sigma=True, Nboot=50, method=default_method,
         compute_params=True, **args):
      '''Make an interpolating template of the lightcurve.  If fitflux, then fit in 
      flux-space.  If do-sigma, do a monte-carlo simulation to get an idea of
      the uncertainties in the final parameters.  Nboot controls the number of
      bootstrap interations.  There are several interpolating methods that
      can be chosen, depending on your python distribution.  To find the 
      methods available, run self.list_types().  The rest of the arguments are 
      passed to the interpolating method.  Upon successful completion of the
      routine, the following member variables will be populated: 
         Tmax, e_Tmax: time of maximum (and error if do_sigma=1) 
         Mmax, e_Mmax:   peak magnitude
         dm15, e_dm15:   delta-m 15
         model_type:     'spline','hyperspline','polynomial','chebyshev','hermite',
                         'hermiteE','gp'  (run self.list_types to see available)'''

      if self.parent.Tmax > 0:
         evt = arange(int(self.t[0]), int(self.t[-1])+1)*1.0 + self.parent.Tmax
      else:
         evt = arange(self.MJD[0], self.MJD[-1])*1.0

      x = self.MJD
      if fitflux:
         y = self.flux
         ey = self.e_flux
      else:
         y = self.mag
         ey = self.e_mag

      self.interp = fit1dcurve.Interpolator(method, x, y, ey, self.mask, **args)
      if fitflux:
         self.model_flux = 1
      else:
         self.model_flux = 0

      if compute_params:
         self.compute_lc_params(N=Nboot)

   def compute_lc_params(self, N=50):
      '''Compute dm15, Tmax, Mmax, and covariances for the light-curve.'''

      # rest-frame of the SN
      day15 = 15*(1 + self.parent.z)
      zp = self.filter.zp
      Tmaxs = []
      Mmaxs = []
      dm15s = []
      inter = self.interp
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
               if self.model_flux:
                  curvs = array([-1])
               else:
                  curvs = array([1])
         else:
            xs,ys,curvs = inter.find_extrema()
         if self.model_flux:
            xs = xs[less(curvs,0)]
            ys = ys[less(curvs,0)]
         else:
            xs = xs[greater(curvs,0)]
            ys = ys[greater(curvs,0)]
         if len(xs) == 0:
            if i == 0:
               # If we can't find maximum on the data, we're done
               for attrib in ['Tmax','Mmax','dm15','e_Tmax','e_Mmax',
                     'e_dm15','cov_Tmax_dm15','cov_Tmax_Mmax',
                     'cov_Mmax_dm15']:
                  self.__dic__[attrib] = None
               return
            Tmaxs.append(-1)
            Mmaxs.append(-1)
            dm15s.append(-1)
            continue
         Tmaxs.append(xs[0])

         # Find Mmax/dm15
         y15,m15 = inter(Tmaxs[-1] + day15)
         if not self.model_flux:
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
      Tmaxs,Mmaxs,dm15s = map(array, [Tmaxs, Mmaxs, dm15s])

      Tgids = greater(Tmaxs, 0)
      if sum(Tgids)*1.0/len(Tgids) < 0.8:
         print "Warning!  More than 20% of MC realizations had bad Tmax"
      Mgids = greater(Mmaxs, 0)
      if sum(Mgids)*1.0/len(Mgids) < 0.8:
         print "Warning!  More than 20% of MC realizations had bad Mmax"
      dgids = greater(dm15s, 0)
      if sum(dgids)*1.0/len(dgids) < 0.8:
         print "Warning!  More than 20% of MC realizations had bad dm15"
      self.Tmax = Tmaxs[0]
      self.e_Tmax = std(Tmaxs[Tgids]) 
      self.Mmax = Mmaxs[0]
      self.e_Mmax = std(Mmaxs[Mgids]) 
      self.dm15 = dm15s[0]
      self.e_dm15 = std(dm15s[dgids]) 
      self.cov_Tmax_dm15 = \
            sum(compress(Tgids*dgids, (Tmaxs-Tmaxs[0])*(dm15s-dm15s[0])))/\
            (sum(Tgids*dgids) - 1)
      self.cov_Tmax_Mmax = \
            sum(compress(Tgids*Mgids, (Tmaxs-Tmaxs[0])*(Mmaxs-Mmaxs[0])))/\
            (sum(Tgids*Mgids) - 1)
      self.cov_Mmax_dm15 = \
            sum(compress(Mgids*dgids, (Mmaxs-Mmaxs[0])*(dm15s-dm15s[0])))/\
            (sum(Mgids*dgids) - 1)
      return

   def plot(self, device='/XSERVE', interactive=True, epoch=1, flux=0, gloes=True,
         symbol=4):
      '''Plot this light-curve.  You can specify a PGPLOT device (ignored if using
      matplotlib as plotter), the default is an X server.  If epoch=1, plot times
      relative to self.Tmax.  If flux=1, plot in flux units.  use GLoEs to smooth 
      the data and produce a model.  You can specify the symbol to plot with 'symbol'.'''
      # first, do some cleanup if we're using matplotlib:
      try:
         self.mp.bc.disconnect()
      except:
         pass
      return plotmod.plot_lc(self, device, epoch, flux, symbol)

   def replot(self):
      '''Replot a figure, if it exists and belongs to this instance.'''
      plotmod.replot_lc(self)
