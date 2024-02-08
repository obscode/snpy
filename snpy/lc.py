from __future__ import print_function
from numpy import *
from .filters import fset
import numpy.random as RA
from .utils import fit1dcurve
from .utils import InteractiveFit
from . import plotmod
from scipy.optimize import brentq

if 'gp' in fit1dcurve.functions:
   default_method = 'gp'
else:
   default_method = 'hyperspline'

class lc:
   '''This class contains the data representing a lightcurve.  This means 
   the time (in MJD), magnitudes, uncertainties, fluxes, etc.  It is meant 
   to be contained by the :class:`snpy.sn` class, to make for somewhat more 
   convenient access to the data.

   Args:
      parent (snpy.sn instance): the parent objects to which the data belongs
      band (str): the passband the data was observed through
      MJD (float array): date of observation. Even though it's called MJD,
                         no assumptions are made for epoch zero-point.
      mag (float array): observations in magnitudes
      e_mag (float array): error of observations in magnitudes
      restband (str): If different than the observed passs band.
      K (float array):  If specified, K-corrections will be stored with the
                        magnitudes. Do this to k-correct the data indpendent
                        of fitting.
      SNR (float array):  optional signal-to-noise of the data, if different
                          than 1.087/e_mag.
      sids (int array): optional index array of into sources of photometry.
   '''

   def __init__(self, parent, band, MJD, mag, e_mag, restband=None, K=None, 
         SNR=None, sids=None):
      self.parent = parent     # pointer to the SN containing class
      self.band = band         # observed filter
      self.MJD = MJD           # time (in Modified Julian Day... or whatever)
      self.magnitude = mag     # magnitude
      self._flux = None        # Caching
      self._eflux = None       #    "
      self.e_mag = e_mag       # error in magnitude
      self._SNR = SNR           # Signal-to-noise ratio
      if any(equal(self.e_mag,0)):
         print("Warning:  you have errors that are zero:  setting them to 0.001 mag.")
         print("If you don't like that, fix your errors.")
         self.e_mag = where(equal(self.e_mag, 0), 0.001, self.e_mag)
      self.covar = None        # full error matrix (if errors are correlated)
      self.restband = restband # the observed band in the frame of the SN
      self.K = K               # k-correction to convert mag to restband
      self.debug=0             # Of course, we'll never need this!

      # a way to mask out bad data from fitting
      self.mask = ones(self.mag.shape, dtype=bool) 
      self.filter = fset[self.band] # instance of filter object for this band

      # A way to keep track of sources of photometry
      if sids is None:
         self.sids = zeros(self.mag.shape[0], dtype=int)
      else:
         try:
            self.sids = asarray(sids, dtype=int)
         except:
            raise ValueError("sids must be integer array or list")
         if self.sids.shape != self.mag.shape:
            raise ValueError("sids must have the same shape as photometry")

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

   def __eq__(self, other):
      '''Check to see if another lc instance is equal to this one.'''
      if self.MJD.shape != other.MJD.shape:
         return False
      return all(less(absolute(self.MJD-other.MJD),1e-5)) \
             and all(less(absolute(self.magnitude-other.magnitude),1e-5)) \
             and all(less(absolute(self.e_mag-other.e_mag),1e-5))

   def __ne__(self, other):
      return not self.__eq__(other)

   def get_t(self):
      '''Get the epoch of observations relative to parent.Tmax. This can also
      be accessed through the :attr:`.t` attribute.'''
      try:
         t = self.MJD - self.parent.Tmax
      except:
         t = self.MJD
      return(t)

   def get_flux(self):
      '''Get the flux of the observations. This can also be access through
      the :attr:`.flux` attribute.'''
      if getattr(self, '_flux', None) is None:
         self._flux = power(10.0, -0.4*(self.mag - self.filter.zp))
      return self._flux

   def get_SNR(self):
      '''Get the signal-to-noise ratio of the data. This can also be 
      accessed through the :attr:`.SNR` attribute.'''
      if getattr(self, '_SNR', None) is None:
         return 1.087/self.e_mag
      return self._SNR

   def get_e_flux(self):
      '''Get the flux error of the observations. This can also be access through
      the :attr:`.e_flux` attribute.'''
      if getattr(self, '_eflux', None) is None:
         self._eflux = self.get_flux()*self.e_mag/1.0857
      return self._eflux

   def get_covar(self, flux=1):
      '''returns the error matrix in flux units (unless flux=0).  If this was
      not setup by the user, the variance will be returned as a 1-D array.i
      
      Args:
         flux (bool):  If true, return in flux units, otherwise in mags.
      
      Returns:
         float array: the square covariance error matrix.'''
      if self.covar is not None:
         if flux:
            f = self.get_flux()
            return self.covar*f[newaxis,:]*f[:,newaxis]/1.17874
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
      elif name == "SNR":
         return self.get_SNR()
      elif name == "t":
         return self.get_t()
      elif name == 'mag':
         if self.K is None:
            return(self.magnitude)
         else:
            try:
               return(self.magnitude - self.K)
            except:
               raise AttributeError("Error:  k-corrections and magnitudes incompatible")
      elif 'interp' in self.__dict__ and self.__dict__['interp'] is not None:
         if name in self.__dict__['interp'].pars:
            return getattr(self.__dict__['interp'],name)
      raise AttributeError("Error:  attribute %s not defined" % (name))

   def __setattr__(self, name, value):
      if 'interp' in self.__dict__:
         if self.__dict__['interp'] is not None:
            if name in self.__dict__['interp'].pars:
               setattr(self.__dict__['interp'], name, value)
               # Now check to see if there is an interactive fit to update
               mp = getattr(self, 'mp', None)
               if mp is not None and \
                     isinstance(mp, InteractiveFit.InteractiveFit):
                  mp.redraw()

               #self.replot()
               return
      if name == 'mask':
         if 'mask' in list(self.__dict__.keys()):
            raise TypeError("lc instance's mask attribute must be modifed in-place")
      self.__dict__[name] = value

   def __getstate__(self):
      # Need this because MPL objects are not pickleable.
      odict = self.__dict__.copy()
      if 'mp' in odict:  del odict['mp']
      if 'pts' in odict:  del odict['pts']
      if '_lc_labs' in odict: del odict['_lc_labs']
      if '_eflux' in odict: del odict['_eflux']
      if '_flux' in odict: del odict['_flux']
      return odict

   def __setstate__(self, state):
      # We need this because GP's need mean function, but pickling a 
      # member function is not supported. So if we unpickle a GP, set
      # the mean accordingly.
      self.__dict__.update(state)
      if getattr(self, 'interp',None) is not None:
         if fit1dcurve.GaussianProcess is not None and \
            isinstance(self.interp, fit1dcurve.GaussianProcess):
            if getattr(self.interp, 'mean', None) is None:
               self.interp.mean = self.mean
               self.interp.setup = False

   def time_sort(self):
      ids = argsort(self.MJD)
      self.MJD = take(self.MJD, ids)
      self.magnitude = take(self.magnitude, ids)
      self.e_mag = take(self.e_mag, ids)
      self.sids = take(self.sids, ids)
      if self.K is not None:
         self.K = take(self.K, ids)

   def mask_epoch(self, tmin, tmax):
      '''Update the lc's mask to only include data between tmin and tmax.i
      
      Args:
         tmin,tmax (float): range over which data is considered good.
         
      Returns:
         None
         
      Effects:
         The mask attribute is updated.'''
      self.mask[greater_equal(self.t, tmax)] = False
      self.mask[less_equal(self.t, tmin)] = False

   def mask_emag(self, emax):
      '''Update the lc's mask to only include data with e_mag < max.
      
      Args:
         emax (float): maximum error for good data
         
      Returns:
         None
         
      Effects:
         The mask attribute is updated.'''
      self.mask[greater(self.e_mag, emax)] = False

   def mask_SNR(self, minSNR):
      '''Update the lc's mask to only include data with SNR > min.
      
      Args:
         minSNR (float): minimum signal-to-noise for good data
         
      Returns:
         None
         
      Effects:
         The mask attribute is updated.'''
      self.mask[less(self.SNR, minSNR)] = False

   def eval(self, times, t_tol=-1, epoch=0):
      '''Interpolate (if required) the data to specific times.  If there is a
      data point "close enough", that value is
      used without interpolation.  If epoch is nonzero, then times are
      interpreted relative to parent.Tmax
      
      Args:
         times (float array): times at which to interpolate
         t_tol (float):  If a data point is less than t_tol away from requested
                         time, the data will be returned. Use -1 to always
                         interpolate.
         epoch (bool): If True, ``lc.parent.Tmax`` is added to times
         
      Returns:
         2-tuple:  (mag,mask)
         
         * mag (float array):  interpolated magnitudes
         * mask (bool array):  True if interpolation is valid
         
      Raises:
         AttributeError:  if an interpolator for the data has not been assigned
      '''

      if self.interp is None:
         # Try to use a model
         if self.band not in self.parent.model._fbands: raise AttributeError("Error.  To interpolate, you need to fit a template or model first.")
         
      times = atleast_1d(times) 
      if epoch: times = times + self.parent.Tmax

      if self.interp is not None:
         evm,mask = self.interp(times)
         if self.model_flux:
            mask = mask*greater(evm,0)
            evm = where(evm <=0, 1.0, evm)
            evm = -2.5*log10(evm) + self.filter.zp
      else:
         evm,eevm,mask = self.parent.model(self.band, times)
         if not self.parent.model.model_in_mags:
            mask = mask*greater(evm, 0)
            evm = where(evm <=0, 1.0, evm)
            evm = -2.5*log10(evm) + self.filter.zp

      if t_tol > 0:
         # Now, we scan t and eval_t and find where they are less than tol.  
         # In these cases, we take the average of any matching times
         delta = absolute(self.MJD[newaxis,:] - times[:,newaxis])
         cond = less(delta, t_tol)
         values = array([self.mag]*len(times))*cond
         N = sum(cond, axis=1)
         w = where(N > 0, N, 1)
         s = sum(values, axis=1)
         evm = where(N > 0, s/w, evm)

      return(evm,mask)

   def spline_fit(self, fitflux=0, do_sigma=1, Nboot=100, keep_boot=1, 
         method='spline2', **args):
      '''Make a spline template of the lightcurve.  The name is a bit misleading
      as many interpolators are available. Once this is done, the light-curve
      can be interpolated and properties (Tmax, Mmax, etc) can be computed.
      
      Args:
         fitflux (bool):  If True, fit the flux instead of magnitudes
         do_sigma (bool):  If True, perform MC iterations to estimate the
                           errors in the interpolations.
         Nboot (int):  Number of MC iterations to perform.
         keep_boot (bool): If True, keep the MC realizations for future use
         method (str):  The interpolation method. Use :meth:`.list_types`
                        to find out what is available.
         args (dict):  all other arguments are sent to the constructor of
                       the interpolator. See :mod:`snpy.utils.fit1dcurve`
                       for documentation on this.

      Returns:
         None

      Effects:
         Upon successful completion of the
         routine, the following member variables will be populated: 
         * Tmax, e_Tmax: time of maximum (and error if do_sigma=1) 
         * Mmax, e_Mmax:   peak magnitude
         * dm15, e_dm15:   delta-m 15
         * model_type:     "spline" or "spline2"
         * tck:            spline info '''
      return self.template(fitflux=fitflux, do_sigma=do_sigma, Nboot=Nboot, 
            method=method, **args)

   def curve_fit(self, fitflux=0, do_sigma=1, Nboot=100, keep_boot=1, 
         method='spline2', **args):
      '''Make a spline template of the lightcurve.  The name is a bit misleading
      as many interpolators are available. Once this is done, the light-curve
      can be interpolated and properties (Tmax, Mmax, etc) can be computed.
      
      Args:
         fitflux (bool):  If True, fit the flux instead of magnitudes
         do_sigma (bool):  If True, perform MC iterations to estimate the
                           errors in the interpolations.
         Nboot (int):  Number of MC iterations to perform.
         keep_boot (bool): If True, keep the MC realizations for future use
         method (str):  The interpolation method. Use :meth:`.list_types`
                        to find out what is available.
         args (dict):  all other arguments are sent to the constructor of
                       the interpolator. See :mod:`snpy.utils.fit1dcurve`
                       for documentation on this.

      Returns:
         None

      Effects:
         Upon successful completion of the
         routine, the following member variables will be populated: 
         * Tmax, e_Tmax: time of maximum (and error if do_sigma=1) 
         * Mmax, e_Mmax:   peak magnitude
         * dm15, e_dm15:   delta-m 15
         * model_type:     "spline" or "spline2"
         * tck:            spline info '''
      return self.template(fitflux=fitflux, do_sigma=do_sigma, Nboot=Nboot, 
            method=method, **args)


   def list_types(self):
      print("the following methods are available for constructing templates:")
      fit1dcurve.list_types()

   def mean(self, x, flux=False, verbose=False):
      '''A convenience function used with the GP interpolator. If a model
      exists for the LC, use it where it is valid, otherwise, return
      an extrapolation.'''
      scalar = False
      if len(shape(x)) == 0:
         scalar = True
      x = atleast_1d(x)
      if x.shape[0] == 0:
         '''empty array, just return some sensible value'''
         if flux:
            return x*0 + median(self.flux)
         else:
            return x*0 + median(self.mag)


      if len(x.shape) == 2 and x.shape[1] == 1:
         x = x[:,0]
      if self.band in self.parent.model._fbands:
         sids = argsort(x)
         m,em,f = self.parent.model(self.band, x[sids], extrap=False)
         if not all(f):
            # need to do linear interpolation between end of template and
            # remaining points
            m0,em0,f0 = self.parent.model(self.band, self.MJD)
            # First valid point
            mid = argmin(self.MJD[f0])
            xmin = self.MJD[f0][mid]
            if verbose:
               print("extrapolating early data up to ",xmin)
            lids = less(x[sids], xmin)
            if any(lids):
               m1,em1,f1 = self.parent.model(self.band, x[sids][lids],
                     extrap=True)
               m[lids] = m1
            mid = argmax(self.MJD[f0])
            # last valid point on the model
            xmax = self.MJD[f0][mid]
            ymax = m0[f0][mid]

            if verbose:
               print("extrapolating late data from (%f,%f)" % (xmax,ymax))

            lids = greater(x[sids], xmax)
            if any(lids):
               gids = greater_equal(self.MJD, xmax)
               if verbose:
                  print("   using points at", self.MJD[gids])
               a,b = polyfit(self.MJD[gids], self.mag[gids], 1)
               if verbose:
                  print("   median slope: ", a)
               m[lids] = ymax + a*(x[sids][lids]-xmax)
         y = x*0
         put(y, sids, m)
         if flux:
            return power(10.0, -0.4*(y - self.filter.zp))
         return y
      # No model, so we'll just give a reasonable constant mean
      if flux:
         return x*0 + median(self.flux)
      return x*0 + median(self.mag)

   def template(self, fitflux=False, do_sigma=True, Nboot=50, 
         method=default_method, compute_params=True, interactive=False, **args):
      '''A backwrd-compatibility alias for :meth:`.spline_fit`.'''

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

      xx,yy,ee = fit1dcurve.regularize(x, y, ey)
      if len(xx) < 2:
         raise ValueError("Cannot interpolate data with less than two distinct data points")
      if method == 'gp':
         args['mean'] = self.mean
      self.interp = fit1dcurve.Interpolator(method, x, y, ey, self.mask, **args)
      if interactive:
         self.mp = plotmod.launch_int_fit(self, fitflux=fitflux)
      if fitflux:
         self.model_flux = 1
      else:
         self.model_flux = 0

      if compute_params:
         self.compute_lc_params(N=Nboot)

   def compute_lc_params(self, N=50):
      '''Compute dm15, Tmax, Mmax, and covariances for the light-curve.
      
      Args:
         N (int):  the nubmer of Monte-Carlo iterations for computing
                   errors in parameters.
                   
      Returns:
         None

      Effects:
         The foloowing LC instance's member variables are updated:
          - Tmax:  time of primary maximum. None if such does not exist
          - e_Tmax:  error in Tmax
          - dm15:  The decline rate parameter, or None if not measurable
          - e_dm15:  error in dm15
          - Mmax:  the magnitude at maximum, or None if no maximum.
          - e_Mmax:  the error in Mmax
          - cov_Tmax_dm15:  covariance between Tmax and dm15
          - cov_Mmax_dm15:  covariance between Mmax and dm15
          - cov_Tmax_Mmax:  covariance between Tmax and Mmax
      '''

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
            xs,ys,curvs = inter.find_extrema(xmin=Tmaxs[0]-5, xmax=Tmaxs[0]+5)
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
                  self.__dict__[attrib] = None
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
      Tmaxs,Mmaxs,dm15s = list(map(array, [Tmaxs, Mmaxs, dm15s]))

      Tgids = greater(Tmaxs, 0)
      goodfrac = sum(Tgids)*1.0/len(Tgids)
      if sum(Tgids)*1.0/len(Tgids) < 0.8:
         print("Warning!  %f%% of MC realizations had good Tmax" % \
               (sum(Tgids)*100.0/len(Tgids),))
      Mgids = greater(Mmaxs, 0)
      if sum(Mgids)*1.0/len(Mgids) < 0.8:
         print("Warning!  %f%% of MC realizations had good Mmax" %\
               (sum(Mgids)*100.0/len(Mgids),))
      dgids = greater(dm15s, 0)
      if sum(dgids)*1.0/len(dgids) < 0.8:
         print("Warning!  %f%% of MC realizations had good dm15" %\
               (sum(dgids)*100.0/len(dgids),))
      if any(Tgids):
         self.Tmax = Tmaxs[0]
         self.e_Tmax = std(Tmaxs[Tgids]) 
      else:
         self.Tmax = self.e_Tmax = -1
      if any(Mgids):
         self.Mmax = Mmaxs[0]
         self.e_Mmax = std(Mmaxs[Mgids]) 
      else:
         self.Mmax = self.e_Mmax = -1
      if any(dgids):
         self.dm15 = dm15s[0]
         self.e_dm15 = std(dm15s[dgids]) 
      else:
         self.dm15 = self.e_dm15 = -1

      if sum(Tgids*dgids) > 3:
         self.cov_Tmax_dm15 = \
               sum(compress(Tgids*dgids, (Tmaxs-Tmaxs[0])*(dm15s-dm15s[0])))/\
               (sum(Tgids*dgids) - 1)
      
      if sum(Tgids*Mgids) > 3:
         self.cov_Tmax_Mmax = \
               sum(compress(Tgids*Mgids, (Tmaxs-Tmaxs[0])*(Mmaxs-Mmaxs[0])))/\
               (sum(Tgids*Mgids) - 1)
      if sum(Mgids*dgids) > 3:
         self.cov_Mmax_dm15 = \
               sum(compress(Mgids*dgids, (Mmaxs-Mmaxs[0])*(dm15s-dm15s[0])))/\
               (sum(Mgids*dgids) - 1)
      return

   def plot(self, epoch=1, flux=0, gloes=True, symbol=4, outfile=None,
         use_model=True):
      '''Plot this light-curve, including residuals if a model or interpolator
      is defined.  
      
      Args:
         epoch (Bool):  If True, plot times relative to self.Tmax
         flux (Bool):   If True, plot in flux units rather than magnitudes
         gloes (Bool):  If True, use GLoEs to smooth data and produce model
         symbol (misc): Specify matplotlib symbol to use for plotting points
         outfile (string or None):  If not None, output graph to specified file
         use_model (Bool): If True and both model and interpolator are defined,
                           use the model to plot residuals.

      Returns:
         PanelPlot or SimplePlot instance.
            A reference to the plot instance created to show the data.
      '''

      return plotmod.plot_lc(self, epoch, flux, symbol, outfile=outfile,
            use_model=use_model)

   def replot(self):
      '''Replot a figure, if it exists and belongs to this instance.'''
      plotmod.replot_lc(self)
