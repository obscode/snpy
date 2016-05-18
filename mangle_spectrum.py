'''
Color-matching:  modifying a spectrum by multiplication of a smooth function
in order to reproduce the observed photometry. A mangler class is used to
do all the heavy lifting of fitting a smooth function. A base class called 
function is provided and should be sub-classed to implement new smoothing 
functions.  So far, there's tension splines and CCM.'''

import types
import numpy as num
from snpy.filters import fset
from snpy.filters import filter
from snpy.filters import vegaB
import scipy.interpolate
from scipy.optimize import minimize
from scipy.optimize import leastsq
try:
   import snpy.tspack as tspline
except:
   tspline = None
from snpy.utils import deredden
import copy

filts = fset
debug=False

class function:
   '''A function object that represents a way of making a smooth multiplication
   on the spectrum, thereby "mangling" it.  The function has an evaluate
   function as well as a wrapper that is suitable for using with mpfit.'''

   def __init__(self, parent):
      '''
      Args:
         parent (Mangler instance): the mangler that will implement this
                                    function.
      '''
      self.parent = parent

   def setup(self):
      '''Do any initial setup.'''
      pass

   def init_pars(self):
      '''Reset the parameters of the function to initial values.  This should
      cause the function to return the identity ( f(x) = 1.0).'''
      pass

   def set_pars(self, pars):
      pass

   def __eval__(self, x):
      pass

class f_tspline(function):
   '''Tension spline mangler. Tension splines are controlled by a tension
   parameter. The higher the tension, the less curvature the spline can
   have. This is good for keeping the spline from "wigglig" too much.'''

   def __init__(self, parent, tension=None, slopes=None, 
         verbose=False, log=False):
      '''
         Args:
            parent (instance): Mangler that will fit the function.
            tension (float): A tension parameter
            slopes (2-tuple): The end slopes can be kept fixed using this
                              argument (e.g., set to (0,0) to flatten out
                              the spline at both ends
            verbose (bool): debug info
            log (bool):  Not implemented yet.
      '''
      function.__init__(self,parent)
      self.tension = tension       # specify a global tension parameter
      self.slopes = slopes
      if self.slopes is not None:
         try:
            self.slopes = list(self.slopes)         # End slopes if not using gradient
         except:
            raise ValueError, "slopes must be a list type of length 2"
      self.pars = None
      self.verbose = verbose
      self.log = log

   def init_pars(self, nid=0, init=None):
      # one parameter for each knot points
      np = len(self.parent.bands) - 1
      pi = [1.0 for i in range(np)]
      if init is not None:
         if len(init) == len(pi):
            pi[:] = init[:]
         elif len(init) == len(pi)+1:
            pi[:] = init[:-1]
         else:
            print "Warning: initial mangler parameters inconsistent"
      return(pi)

   def set_pars(self, pars):
      self.pars = num.asarray(pars)

   def __call__(self, x):
      if self.parent.ave_waves is None or self.pars is None:
         return x*0+1.0
      p = num.zeros((len(self.pars)+1,))
      p[0] = 1.0
      p[1:] = self.pars
      if self.slopes is None:
         xyds = tspline.tspsi(self.parent.ave_waves, p,
               tension=self.tension, curvs=[0,0])
      else:
         xyds = tspline.tspsi(self.parent.ave_waves, p,
               tension=self.tension, slopes=self.slopes)

      res = num.array([tspline.tsval1(x[i], xyds) for i in range(x.shape[0])])

      y0 = tspline.tsval1(xyds[0][0], xyds)
      y1 = tspline.tsval1(xyds[0][-1], xyds)
      if self.slopes is not None:
         #dy0 = tspline.tsval1(self.parent.ave_waves[0], xyds, 1)[0]
         #dy1 = tspline.tsval1(self.parent.ave_waves[-1], xyds, 1)[0]
         res = num.where(x < xyds[0][0], y0 + self.slopes[0]*(x-xyds[0][0]), res)
         res = num.where(x > xyds[0][-1], y1 + self.slopes[1]*(x-xyds[0][-1]), res)
      else:
         # end slopes
         s0 = tspline.tsval1(xyds[0][0], xyds, 1)
         s1 = tspline.tsval1(xyds[0][-1], xyds, 1)
         res = num.where(x < xyds[0][0], y0 + s0*(x-xyds[0][0]), res)
         res = num.where(x > xyds[0][-1],y1 + s1*(x-xyds[0][-1]), res)
      return(res)

class f_spline(function):
   '''Spline mangler. Splines are controlled by a smoothing
   parameter. The higher the smoothing, the less curvature the spline can
   have. This is good for keeping the spline from "wigglig" too much.'''

   def __init__(self, parent, s=0, slopes=None, verbose=False, log=False):
      '''
         Args:
            parent (instance): Mangler that will fit the function.
            s (float): A smoothing parameter
            extrapolate (string): 'constant':  extrapolate using constants
                                  'gradient':  extrapolate using end slopes
            verbose (bool): debug info
            log (bool):  Not implemented yet.
      '''
      function.__init__(self,parent)
      self.s = s                   # specify a global smoothness parameter
      #self.gradient = gradient    # constrain the anchors by slope?
      self.slopes = slopes         # End slopes ((0,0) to flatten at ends)
      self.pars = None
      self.verbose = verbose
      self.log = log

   def init_pars(self, nid=0, init=None):
      # one parameter for each knot points, but keep last fixed at 1.0
      np = len(self.parent.bands) - 1
      pi = num.array([1.0 for i in range(np)])
      if init is not None:
         if len(init) == len(pi):
            pi[:] = init[:]
         elif len(init) == len(pi)+1:
            pi[:] = init[:-1]
         else:
            print "Warning: initial mangler parameters inconsistent"
      return(pi)

   def set_pars(self, pars):
      self.pars = num.asarray(pars)

   def __call__(self, x):
      if self.parent.ave_waves is None or self.pars is None:
         return x*0+1.0
      p = num.zeros((len(self.pars)+3,))
      p[1] = 1.0    # normalized to first band (since only interested in colors)
      p[2:-1] = self.pars  # these control the shape. p[0] and p[-1] control slopes
      ws = num.zeros((len(self.parent.ave_waves)+2,))
      ws[1:-1] = self.parent.ave_waves
      ws[0] = ws[1]-(ws[2]-ws[1])
      ws[-1] = ws[-2] + (ws[-2]-ws[-3])
      print ws

      dx1 = ws[2] - ws[1]
      dx2 = ws[-3] - ws[-2]
      if self.slopes is None:
         # naturally extend the slope, make curvature 0
         p[0] = p[1] + (p[2]-p[1])/dx1*(ws[0] - ws[1])
         p[-1] = p[-2] + (p[-3]-p[-2])/dx2*(ws[-1] - ws[-2])
      else:
         # Artificially set the slopes
         p[0] = p[1] + self.slopes[0]*(ws[1] - ws[0])
         p[-1] = p[-2] + self.slopes[1]*(ws[-1] - ws[-2])
         
      xyds = scipy.interpolate.splrep(ws, p, w=p*0+1, s=self.s)

      res = num.array([scipy.interpolate.splev(x[i], xyds) \
            for i in range(x.shape[0])])

      return(res)

class f_ccm(function):
   '''The Cardelli, Clayton, and Mathis (1989) reddening law as a smooth
   fucntion. This would be what dust does to a spectrum, so is a good
   guess, but has much less freedom than tension splines.'''

   def __init__(self, parent, tension=None, gradient=False, slopes=None, verbose=False):
      '''
         Args:
            parent (instance): Mangler instance that will handle the function
            tension (bool):  not used
            gradient (bool):  not used
            slopes (bool):  not used
            verbose (bool):  extra info
      '''
      self.pars = [0.0, 3.1]
      self.verbose = verbose

   def init_pars(self, nid=0, init=None):
      # one parameter for each of ebv, rv
      if init is not None:
         if len(init) == 2:
            return init
      return([0.0, 3.1])

   def set_pars(self, pars):
      self.pars = pars

   def __call__(self, x):
      m = num.array([deredden.unred(x[i], x[i]*0+1, -self.pars[0], 
         R_V=self.pars[1])[0] for i in range(x.shape[0])])
      #m *= self.pars[0]
      return m

mangle_functions = {'spline':f_spline,
                    'tspline':f_tspline,
                    'ccm':f_ccm}


class mangler:
   '''Given a spectrum and spectrum object, find the best set of a function's 
   paramters, such that the function multiplied by the spectrum produces
   the colors specified.'''

   def __init__(self, wave, flux, method, normfilter=None, z=0, **margs):
      '''
         Args:
            wave (float array): input wavelength in Angstroms. Can be a
                                list for multiple spectra
            flux (float array): associated fluxes in arbitrary units
            method (str): The method used (should be a key of
                             mangle_functions:  'tspline' or 'ccm'
            normfilter (str): the mangled function will be normalized
                              to the observed flux through this filter.
                              Defaut is to take the filter with the
                              most valid observations
            z (float): The redshift of the object. The spectrum is redshifted
                       by this amount before doing the mangling.
            **margs (dict): All other argumements are sent to the mangling
                      function's __init__().'''

      self.wave = num.asarray(wave)
      self.flux = num.asarray(flux)
      self.normfilter = normfilter
      if len(num.shape(self.wave)) == 1:
         self.wave = self.wave.reshape((1,self.wave.shape[0]))
      if len(num.shape(self.flux)) == 1:
         self.flux = self.flux.reshape((1,self.flux.shape[0]))
      self.ave_waves = None
      if method not in mangle_functions:
         methods = ",".join(mangle_funcitons.keys())
         raise ValueError, "method must be one of the following:\n%s" % methods
      self.function = mangle_functions[method](self, **margs)
      self.z = z
      self.verbose=margs.get('verbose', False)

   def get_colors(self, bands):
      '''Given a set of filters, determine the colors of the mangled
      spectrum for the current set of function parameters.  You'll get
      N-1 colors for N bands and each color will be bands[i+1]-bands[i]
      
      Args:
         bands (list of str): The list of filters to construct colors from
         
      Returns
         2d array: colors'''
      mflux = self.function(self.wave)*self.flux
      cs = self.resp_rats*0
      for i in range(cs.shape[0]):
         for j in range(cs.shape[1]):
            m1 = fset[bands[j]].synth_mag(self.wave[i], mflux[i], z=self.z)
            m2 = fset[bands[j+1]].synth_mag(self.wave[i], mflux[i], z=self.z)
            cs[i,j] = m1-m2
      return cs

   def get_mflux(self):
      '''Given the current paramters, return the mangled flux.'''
      #if self.verbose:  print 'calling mangling function'
      mflux = self.function(self.wave)
      #if self.verbose:  print 'done'
      return(self.flux*mflux)

   def create_anchor_filters(self, bands, anchorwidth):
      '''Given a set of filters in bands, create two fake filters at
      either end to serve as anchor filters. These can be useful for
      controling splines that don't naturally allow their end slopes
      to be constrained.
      
      Args:
         bands (list of str): The list of filters to construct colors from
         anchorwidth (float): distance from reddest and bluest filter
                              in Angstroms.
      
      Returns:
         None

      Effects:
         fset is updated to include two new filters: 'red_anchor' and
         'blue_anchor'
      '''
      mean_wave = num.array([fset[b].ave_wave for b in bands])
      red = bands[num.argmax(mean_wave)]
      blue = bands[num.argmin(mean_wave)]

      wave0 = fset[blue].wave[0]
      wave1 = fset[red].wave[-1]

      # Create two new fake filters
      fset['blue'] = filter('blue_anchor')
      fset['red'] = filter('red_anchor')
      resp = num.array([0.,0.,1.,1.,0.,0.])
      dwave = num.array([anchorwidth+2., anchorwidth+1., anchorwidth, 3., 2., 1.])
   
      filts['blue'].wave = wave0 - dwave
      filts['blue'].resp = resp*1.0
      filts['red'].wave = wave1 + dwave[::-1]
      filts['red'].resp = resp
   
      # Setup zero points to reasonable values
      filts['red'].zp = filts['red'].compute_zpt(vegaB, 0.0)
      filts['blue'].zp = filts['blue'].compute_zpt(vegaB, 0.0)
   
      filts['red'].ave_wave = num.sum(filts['red'].wave)/len(filts['red'].wave)
      filts['blue'].ave_wave = num.sum(filts['blue'].wave)/len(filts['blue'].wave)
   
      wave0 = filts['blue'].wave[0]
      wave1 = filts['red'].wave[-1]
   
      if wave0 < self.wave.min()*(1+self.z) or wave1 > self.wave.max()*(1+self.z): 
         print 'Problem in mangle_spectrum: SED does not cover anchor filter '+\
               'definitions'
         if(self.verbose): 
            print 'Anchor filter definitions cover  ',wave0, 'A to ',wave1,'A'
            print 'SED covers ',self.wave.min()*(1+self.z),'A to ',\
                  self.wave.max()*(1+self.z),'A (observed)'
      if(self.verbose): print 'Anchor filter definitions cover  ', \
                               wave0,'A to ',wave1,'A'
      

   def solve(self, bands, colors, fixed_filters=None, 
         anchorwidth=100, xtol=1e-10, ftol=1e-4, gtol=1e-10,
         init=None):
      '''Solve for the mangling function that will produce the observed colors
      in the filters defined by bands.
      
      Args:
         bands (list of str): The list of filters to construct colors from
         colors (float array): the colors to match. If 1d, the N-1 colors
                               should correspond to bands[i-1]-bands[i],
                               if 2d, the first index should reflect spectrum
                               number.
         fixed_filters (2-tuple or None): The filters whose control points
                              should remain fixed. If None, two fake filters
                              will be created.
         anchorwidth (float): If making fake filters, they will be spaced
                              this many Angstroms from the bluest and reddest
                              filters (default 100).
         xtol, ftol, gtol (float): see scipy.optimize.leastsq.'''

      # going to try to do multiple-epoch colors simultaneously with one
      #  mangle function.  colors[i,j] = color for epoch i, color j
      colors = num.asarray(colors)
      if len(num.shape(colors)) == 1:
         colors = colors.reshape((1,colors.shape[0]))

      if colors.shape[0] != self.wave.shape[0]:
         raise ValueError, "colors.shape[0] and number of spectra must match"
      if colors.shape[1] != len(bands) - 1:
         raise ValueError, "colors.shape[1] must be len(bands)-1"
      self.gids = num.less(colors, 90)

      self.bands = bands
      self.mfluxes = num.zeros((colors.shape[0],len(bands)), dtype=num.float32)
      self.mfactors = num.array([num.power(10,-0.4*fset[b].zp) for \
            b in self.bands])

      if len(bands) != colors.shape[1] + 1:
         raise ValueError, "length of bands must be one more than colors"

      if self.normfilter is None:
         num_good = num.sum(self.gids*1, axis=0)
         self.normfilter = bands[num.argmax(num_good)+1]
      else:
         if self.normfilter not in bands:
            raise ValueError, "normfilter must be one of bands"
      self.nid = bands.index(self.normfilter)

      # The average wavelengths of the filters being used
      self.ave_waves = num.array([fset[b].ave_wave for b in self.bands])
      minb = self.bands[num.argmin(self.ave_waves)]
      maxb = self.bands[num.argmax(self.ave_waves)]
      
      ## Anchor wavelengths if method requires them
      #w0 = fset[minb].wave.min() - 100
      #w1 = fset[maxb].wave.max() + 100
      #self.anchorwaves = num.concatenate([[w0],self.ave_waves,[w1]])

      # Construct the flux levels we want from the colors
      flux_rats = num.power(10, 0.4*colors)    # M X N-1 of these
      dzps = num.array([fset[bands[i]].zp - fset[bands[i+1]].zp \
            for i in range(0,len(bands)-1)])
      self.resp_rats = num.power(10, -0.4*(colors - dzps[num.newaxis,:]))
      self.resp_rats = num.where(self.gids, self.resp_rats, 1)
      id = self.bands.index(self.normfilter)
      if self.verbose:
         print "You input the following colors:"
         for i in range(colors.shape[1]):
            print "%s - %s" % (bands[i],bands[i+1])
         print colors
         print "This translates to the following response ratios:"
         print self.resp_rats[self.gids]
         print "Initial colors for this spectrum are:"
         print self.get_colors(bands)
         self.init_colors = self.get_colors(bands)


      # Do some initial setup
      pi = self.function.init_pars(nid=id, init=init)


      if self.verbose:
         quiet = 0
      else:
         quiet = 1
      result = leastsq(self.leastsq, pi, args=(bands,), xtol=1e-4, ftol=1e-4,
            gtol=1e-4, full_output=True)
      self.function.set_pars(result[0])
      if self.verbose:
         print "The final colors of the SED are:"
         print self.get_colors(bands)
         print "Compared to initial colors:"
         print self.init_colors

      return(result)

   def leastsq(self, p, bands):
      self.function.set_pars(p)
      mflux = self.get_mflux()
      mresp = self.resp_rats*0
      filts = num.arange(mresp.shape[1])
      filts = num.repeat(filts, mresp.shape[0]).reshape(mresp.shape[::-1]).T
      for i in range(mresp.shape[0]):
         for j in range(mresp.shape[1]):
            if not self.gids[i,j]:  continue
            mresp[i,j] = fset[bands[j]].response(self.wave[i], mflux[i], z=self.z)/\
                         fset[bands[j+1]].response(self.wave[i], mflux[i], z=self.z)

      delt = self.resp_rats[self.gids] - mresp[self.gids]
      #return sum(num.absolute(delt))
      #return (0, num.ravel(delt))
      return num.ravel(delt)

messages = ['Bad input parameters','chi-square less than ftol',
      'paramters changed less than xtol',
      'chi-square less than ftol *and* paramters changed less than xtol',
      'cosine between fvec and jacobian less than gtol',
      'Exceeded maximum number of iterations',
      'ftol is too small','xtol is too small','gtol is too small']

def mangle_spectrum2(wave, flux, bands, mags, fixed_filters=None, 
      normfilter=None, z=0, verbose=0, anchorwidth=100,
      method='tspline', xtol=1e-6, ftol=1e-6, gtol=1e-6, init=None, **margs):
   '''Given an input spectrum, multiply by a smooth function (aka mangle)
   such that the synthetic colors match observed colors.

   Args:
      wave (float array):  Input wavelengths in Angstroms
      flux (float array):  Input fluxes in arbitrary units
      bands (list of str): list of observed filters
      mags (float array): Observed magnitudes
      m_mask (bool array): mask array indicating valid magnitudes.
      fixed_filters (str or None):  If not None, append fake filters on the
                                    red and/or blue end of the spectrum, keeping
                                    their mangled flux fixed. Specify 
                                    'red', 'blue', or 'both'.
      normfilter (str): If specified, the resulting mangled spectrum is 
                        normalized such that the flux through normfilter is
                        the same as the input spectrum.

      z (float):  redshift
      verbose (bool): output extra info.
      anchorwidth (float): width of the "fixed filters" at either end of
                           the SED in Angstroms.
      method (str): specify the method (function form) with which to mangle
                    the spectrum. 'tspline' for tension spline or 'ccm' for
                    the CCM reddening law.
      xtol,ftol,gtol (float):  used to define the tolorance for the 
                    goodness of fit. See ``scipy.optimize.leastsq`` for
                    meaning of these parameters.
      init (list/array): initial values for the mangle parameters.
      margs (dict): All additional arguments to function are sent to 
                    the :class:`mangle_spectrum.Mangler` class.
      '''
   m = mangler(wave, flux, method, z=z, verbose=verbose, 
         normfilter=normfilter, **margs)
   if len(num.shape(mags)) == 1:
      oned = True
      gids = num.less(mags[:-1],90)*num.less(mags[1:],90)
      colors = num.where(gids, mags[:-1]-mags[1:], 99.9)
   else:
      oned = False
      gids = num.less(mags[:-1,:],90)*num.less(mags[1:,:],90)
      colors = num.where(gids, mags[:-1,:]-mags[1:,:], 99.9)
   result = m.solve(bands, colors, fixed_filters=fixed_filters,
         anchorwidth=anchorwidth, xtol=xtol, ftol=ftol, gtol=gtol,
         init=init)
   if result[4] > 4 or verbose:
      print "Warning: %s" % result[5]

   # finally, normalize the flux
   mflux = m.get_mflux()
   id = bands.index(m.normfilter)
   if oned:
      omag = num.array([mags[id]])
   else:
      omag = mags[id,:]
   for i in range(len(mflux)):
      mmag = fset[m.normfilter].synth_mag(wave,mflux[i],z=z)
      mflux[i] = mflux[i]*num.power(10,-0.4*(omag[i]-mmag))
   if len(m.ave_waves) > len(m.function.pars):
      pars = num.concatenate([m.function.pars,[1.0]])
   else:
      pars = m.function.pars

   return (mflux, m.ave_waves, pars)

def apply_mangle(wave,flux,sw,sf,method='tspline', **margs):

   m = mangler(wave, flux, method, **margs)
   m.ave_waves = sw
   m.function.set_pars(sf)
   return(m.get_mflux())
