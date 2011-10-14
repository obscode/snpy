## Automatically adapted for numpy.oldnumeric Feb 04, 2009 by ipython

from numpy.oldnumeric import *
import scipy.interpolate
from utils import fit_spline
import spline2
from filters import fset
import numpy.oldnumeric.random_array as RA
from utils import stats
from utils import GLOEs
from numpy import bool
from numpy import diag
import plotmod

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
      - self.filter         reference to the filter instance
      - self.tck            Spline representation (knots, coefficients, order)
      - self.spline()       Return a spline interpolated representation of the 
                            data (experimental)'''
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

      self.mask = ones(self.mag.shape, typecode=bool) # a way to mask out bad data from fitting
      self.filter = fset[self.band] # instance of filter object for this band

      self.tck = None          # Spline solution
      #self.model = None        # a model
      #self.model_t = None      #  time for the model
      #self.model_ks = None     #  k-corrections for the model
      self.model_type = None   #  How the model was constructed:  'template',
      #                         #  'template_colors', or 'spline'
      self.model_flux = 0       # is the model in flux-space?
      self.tmin = None         # range over which the model is valid
      self.tmax = None

      self.Tmax = None
      self.Mmax = None
      self.e_Tmax = None
      self.e_Mmax = None

      self.line = GLOEs.line(self.MJD, self.flux, self.e_flux)

      try:
         self.line.N = self.line.find_N(3,10)
      except:
         self.line.sigma = abs(max(self.line.xdata)-min(self.line.xdata))

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
      #elif name == 'model_ks':
      #   if 'model_ks' in self.__dict__:
      #      return self.__dict__['model_ks']
      #   else:
      #      return None
      else:
         raise AttributeError, "Error:  attribute %s not defined" % (name)

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
      self.mask = greater_equal(self.t, tmin)*less_equal(self.t, tmax)

   def mask_emag(self, max):
      '''Update the lc's mask to only include data with e_mag < max.'''
      self.mask = self.mask*less_equal(self.e_mag, max)

   def mkspline(self, knots=None, k=3, s=None, fitflux=0, tmin=None, tmax=None, 
         task=-1, anchor_dist=[0,0], slopes=[None,None]):
      '''Generate a spline interpolation of the data.  You can specify the 
      knots, the order of the spline (k), the smoothing factor, and the task
      to be performed (see scipy documentation).  If you set fitflux nonzero,
      the fit is performed in flux space.  Use tmin and tmax to restrict 
      the portion of the lightcurve to be splined (default whole range).  If you 
      want to control the slopes at the endpoints, then specify them in slopes
      keyword (a 2-tuple).  This effectively adds artificial points at the
      beginning and end, anchor_dist away from the first(last) point.'''
      MJD = compress(self.mask, self.MJD)
      if len(MJD) < 2:
         raise RuntimeError, "Error:  You need at least two valid points to spline"
      mag = compress(self.mask, self.mag)
      if tmax is None: tmax = MJD[-1]
      if tmin is None: tmin = MJD[0]
      self.tmax = tmax
      self.tmin = tmin
      if s is not None:  self.s = s
      e_mag = compress(self.mask, self.e_mag)
      (self.tck,self.fp,ier,mesg) = fit_spline.make_spline(MJD, mag, 
            e_mag, knots=knots, k=k, s=s, fitflux=fitflux, 
            zpt=self.filter.zp, tmin=tmin, tmax=tmax, task=task,
            anchor_dist=anchor_dist, slopes=slopes)
      self.rchisq = self.fp/len(self.mag)

   def mkspline2(self, fitflux=0, k=3, tmin=None, tmax=None, acfsearch=0,
         acffunc='exp', ksi=None, n=None, allownonopt=1, lopt=None,
         rejlev=0.05, xlog=0, adaptive=0, max_curv_fac=10, verbose=0):
      '''Generate a spline interpolation of the data using the Thijsse et al.'s
      "hyper-spline" routine.  You can specify the order of the spline (k), and
      whether to fit in fluxspace (fitflux=1).  The rest of the arguments are more
      advanced.  See docstring of spline2.spline2 for more info.'''
      self.s = -1
      MJD = compress(self.mask, self.MJD)
      if len(MJD) < 2:
         raise RuntimeError, "Error:  You need at least two valid points to spline"
      mag = compress(self.mask, self.mag)
      e_mag = compress(self.mask, self.e_mag)
      if tmax is None: tmax = MJD[-1]
      if tmin is None: tmin = MJD[0]
      self.tmax = tmax
      self.tmin = tmin
      (self.tck,self.fp,ier,mesg) = fit_spline.make_spline2(MJD, mag, 
            e_mag, k=k, fitflux=fitflux, zpt=self.filter.zp, tmin=tmin, 
            tmax=tmax, acfsearch=acfsearch, acffunc=acffunc,
            ksi=ksi, n=n, allownonopt=allownonopt, lopt=lopt,rejlev=rejlev,
            xlog=xlog, adaptive=adaptive, max_curv_fac=max_curv_fac, 
            verbose=verbose)
      self.rchisq = self.fp/len(self.mag)
      if adaptive:  self.lopt = ier

   def eval(self, times, t_tol=0.1, epoch=0):
      '''Interpolate (if required) the data to time 'times'.  If there is a data point 
      closer than t_tol away from a requested time, that value is used without 
      interpolation.  If epoch is nonzero, then times are interpreted relative
      to parent.Tmax'''
      if not len(shape(times)):
         scalar = 1
         times = [times]
      else:
         scalar = 0
      times = asarray(times)
      if epoch: times = times + self.parent.Tmax

      if self.model_type == 'spline':
         if self.tck is None:
            raise AttributeError, "Error:  model_type is spline, but no tck defined."
         evm = scipy.interpolate.splev(times, self.tck)         
         mask = greater_equal(times, self.tck[0][0])*less_equal(times, self.tck[0][-1])
         # TUrns out splev will take a length 1 array and silently return a scalar!  
         #  Check for this and convert if needed         
         if len(shape(evm)) == 0: evm = array([evm])
         if self.model_flux:
            # handle zero or negative flux case
            mask = mask*greater(evm, 0)
            evm = where(evm <= 0, 1.0, evm)
            evm = -2.5*log10(evm) + self.filter.zp
      elif self.model_type == 'spline2':
         if self.tck is None:
            raise AttributeError, "Error:  model_type is spline, but no tck defined."
         evm = spline2.evalsp(times, self.tck)
         mask = greater_equal(times, self.tck[0][0])*less_equal(times, self.tck[0][-1])
         if self.model_flux:
            # handle zero of negaive flux case
            mask = mask*greater(evm, 0)
            evm = where(evm <= 0, 1.0, evm)
            evm = -2.5*log10(evm) + self.filter.zp
      else:
         # do it via the GLOEs interpolation
         if sometrue(self.line.ydata - self.flux):
            self.line.ydata = self.flux
            self.line.dydata = self.e_flux
         dum = self.line.eval(times)
         # Look for negative flux and strictly accept only interpolation
         mask = greater(self.line.y, 0)*less_equal(times, max(self.line.xdata))*\
                greater_equal(times, min(self.line.xdata))
         evm = where(self.line.y <= 0, 1.0, self.line.y)
         evm = -2.5*log10(evm) + self.filter.zp

      if t_tol is not None:
         # Now, we scan t and eval_t and find where they are less than tol.  
         # In these cases, we take the average of any matching times
         delta = absolute(self.MJD[NewAxis,:] - times[:,NewAxis])
         cond = less(delta, t_tol)
         values = array([self.mag]*len(times))*cond
         N = sum(cond, axis=1)
         w = where(N > 0, N, 1)
         s = sum(values, axis=1)
         evm = where(N > 0, s/w, evm)

      if scalar:
         return(evm[0],mask[0])
      else:
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

   def template(self, fitflux=0, do_sigma=1, Nboot=100, keep_boot=1, method='spline', 
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

      if self.parent.Tmax > 0:
         evt = arange(int(self.t[0]), int(self.t[-1])+1)*1.0 + self.parent.Tmax
      else:
         evt = arange(self.MJD[0], self.MJD[-1])*1.0
      magnitude = self.magnitude

      if do_sigma:
         # randomize in flux space:
         flux = power(10, -0.4*(self.mag - self.filter.zp))
         cov_mat = self.get_covar(flux=1)
         if len(shape(cov_mat)) == 1:
            # Make one out of the errors:
            cov_mat = diag(cov_mat)
         realizations = RA.multivariate_normal(flux, cov_mat, Nboot)
         realizations = where(realizations <= 0, 1.0, realizations)
         realizations = -2.5*log10(realizations)
         temp = realizations*0 + magnitude[NewAxis,:]
         realizations = where(realizations == 0, temp, realizations+self.filter.zp)
         realizations = concatenate([[magnitude],realizations])
         #[realizations.append(take(self.magnitude, ids[i])) for i in range(Nboot)]
         #[realizations_t.append(take(self.MJD, ids[i])) for i in range(Nboot)]
      else:
         realizations = [magnitude]

      # Store the results
      models = []
      Tmaxes = []
      dm15s = []
      Mmaxes = []
      tcks = []
      self.model_type = method
      self.model_args = args.copy()

      for i in range(len(realizations)):
         self.magnitude = realizations[i]
         #self.MJD = realizations_t[i]
         if i == 0:
            if method == 'spline':
               self.mkspline(fitflux=fitflux,**args)
            else:
               self.mkspline2(fitflux=fitflux,**args)
            otck = self.tck
         else:
            if method == 'spline':
               args['task'] = -1
               args['knots'] = otck[0][4:-4]
               self.mkspline(fitflux=fitflux, **args)
            else:
               if 'adaptive' in args:
                  args['adaptive'] = 0
                  args['lopt'] = self.lopt
               self.mkspline2(fitflux=fitflux, **args)
            tcks.append(self.tck)
         if method == 'spline':
            evm = scipy.interpolate.splev(evt, self.tck)
         else:
            evm = spline2.evalsp(evt, self.tck)

         if fitflux:
            self.model_flux = 1
         else:
            self.model_flux = 0
 
         models.append(evm)
         if method == 'spline':
            derivs = scipy.interpolate.splev(evt, self.tck, 1)
            tck = scipy.interpolate.splrep(evt, derivs, k=3, s=0)
            roots = scipy.interpolate.sproot(tck)
            if len(roots) == 0:
               # on the decline
               self.Tmax = 0
               self.Mmax = -1
               self.dm15 = -1
            else:
               self.Tmax = None
               for root in roots:
                  curve = scipy.interpolate.splev(root, self.tck, 2)
                  if curve > 0 and not fitflux:
                     self.Tmax = root
                     break
                  elif curve < 0 and fitflux:
                     self.Tmax = root
                     break
               if self.Tmax is not None:
                  T15 = self.Tmax + 15*(1.0 + self.parent.z)
                  res,mask = self.eval(array([self.Tmax, T15]))
                  if mask[0] == 1:
                     self.Mmax = res[0]
                     if mask[1] == 1:
                        self.dm15 = res[1] - res[0]
                     else:
                        self.dm15 = -1
                  else:
                     self.Mmax = -1
                     self.dm15 = -1
               else:
                  self.Tmax = 0
                  self.dm15 = -1
                  self.Mmax = -1
         else:
            xs,ys,signs = spline2.eval_extrema(self.tck)
            if fitflux:
               gids = less(signs, 0)
            else:
               gids = greater(signs, 0)
            if len(nonzero(gids)) == 0:
               self.Tmax = 0
               self.dm15 = -1
               self.Mmax = -1
            else:
               id = nonzero(gids)[0]
               self.Tmax = xs[id]
               self.Mmax = ys[id]
               T15 = self.Tmax + 15*(1.0 + self.parent.z)
               res,mask = self.eval(array([self.Tmax, T15]))
               if mask[0] == 1:
                  self.Mmax = res[0]
                  if mask[1] == 1:
                     self.dm15 = res[1] - res[0]
                  else:
                     self.dm15 = -1
               else:
                  self.Mmax = -1
                  self.dm15 = -1
         if i == 0 or (self.Tmax > 0 and self.Mmax > 0 and self.dm15 > 0):
            Tmaxes.append(self.Tmax)
            Mmaxes.append(self.Mmax)
            dm15s.append(self.dm15)

      models,Tmaxes,Mmaxes,dm15s = map(array,[models,Tmaxes,Mmaxes,dm15s])
      # Now we update 
      gids = greater_equal(evt, self.tmin)*less_equal(evt, self.tmax)
      self.model = compress(gids, models[0])
      self.model_t = compress(gids, evt)
      self.model_mask = ones(self.model.shape, typecode=bool)
      if fitflux:
         try:
            self.model = -2.5*log10(self.model) + self.filter.zp
         except:
            print "Looks like your spline has gone crazy and is predicting"
            print " negative fluxes.  Try increasing the smoothing s"
            return

      self.magnitude = realizations[0]
      self.tck = otck
      self.Tmax = Tmaxes[0]
      self.dm15 = dm15s[0]
      self.Mmax = Mmaxes[0]

      #if we're doing errors:
      if do_sigma:
         # One-sigma index when sorted
         #sig1 = int(len(models)*0.68)
         #medians = stats.median(models)
         #self.model_sigmas = sort(absolute(models-medians), axis=0)[sig1,:]
         self.model_sigmas = stats.sigma(models)
         if fitflux:
            self.model_sigmas = 1.0857*self.model_sigmas/models[0]
         self.model_sigmas = compress(gids, self.model_sigmas)
         gids = greater(Tmaxes, 0)
         self.e_Tmax = 1.49*stats.median(absolute(compress(gids, Tmaxes)-self.Tmax))
         gids = greater(dm15s, 0)
         self.e_dm15 = 1.49*stats.median(absolute(compress(gids, dm15s)-self.dm15))
         gids = greater(Mmaxes, 0)
         self.e_Mmax = 1.49*stats.median(absolute(compress(gids,Mmaxes)-self.Mmax))

         # Now report on covariances:
         gids = greater(Tmaxes,0)*greater(dm15s,0)
         self.cov_Tmax_dm15 = sum(compress(gids, 
            (Tmaxes-self.Tmax)*(dm15s-self.dm15)))/(sum(gids) - 1)
         gids = greater(Tmaxes,0)*greater(Mmaxes,0)
         self.cov_Tmax_Mmax = sum(compress(gids, 
            (Tmaxes-self.Tmax)*(Mmaxes-self.Mmax)))/(sum(gids) - 1)
         gids = greater(Mmaxes,0)*greater(dm15s,0)
         self.cov_Mmax_dm15 = sum(compress(gids, 
            (dm15s-self.dm15)*(Mmaxes-self.Mmax)))/(sum(gids) - 1)

         # Lastly, if we are to save the spline realizations:
         if keep_boot:
            self.tcks = tcks
      return

   def plot(self, device='/XSERVE', interactive=0, epoch=1, flux=0, gloes=0,
         symbol=4):
      '''Plot this light-curve.  You can specify a PGPLOT device (ignored if using
      matplotlib as plotter), the default is an X server.  I epoch=1, plot times
      relative to self.Tmax.  If flux=1, plot in flux units.  If gloes=1 and there
      is no spline or model defined, plot GLOEs model.  If interactive=1, the
      following bindings are set:
      - 'click':  print coordinates
      - 'a':      add spline knot at location (only for method='spline')
      - 'd':      delete closest knot point (only for method='spline')
      - 'm':      move knot point (only for method='spline')
      - 's/S':    make smoothing smaller (s) or larger (S) (only for spline)
      - 'n/N':    make GLOEs N smaller (n) or larger
      - 'x':      select x range
      - 'b':      select box range
      - 'y':      select y range
      - 'q':      quit interactive mode.
      plot in flux space.  If epoch=1, plot time relative to Tmax.  If gloes=1, 
      use GLoEs to smooth the data and produce a model.  You can specify the
      symbol to plot with 'symbol'.'''
      # first, do some cleanup if we're using matplotlib:
      try:
         self.mp.bc.disconnect()
      except:
         pass
      return plotmod.plot_lc(self, device, interactive, epoch, flux, gloes, symbol)
