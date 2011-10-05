'''
Color-matching:  modifying a spectrum by multiplication of a smooth function
in order to reproduce the observed photometry.'''

import types
import numpy.oldnumeric as num
from snpy.filters import fset
from snpy.filters import filter
from snpy.filters import vegaB
import scipy.interpolate
try:
   import snpy.tspack as tspline
except:
   tspline = None
from snpy.utils import deredden
from snpy.utils import mpfit
import copy

filts = fset
debug=False

class function:
   '''A function object that represents a way of making a smooth multiplication
   on the spectrum, thereby "mangling" it.  The function has an evaluate
   function as well as a wrapper that is suitable for using with mpfit.'''

   def __init__(self, parent):
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

   def __init__(self, parent, tension=None, gradient=False, slopes=None, verbose=False):
      function.__init__(self,parent)
      self.knots = None
      self.factors = None
      self.tension = tension       # specify a global tension parameter
      self.gradient = gradient     # constrain the anchors by slope?
      self.slopes = None           # End slopes if not using gradient
      self.pars = None
      self.verbose = verbose

   #def setup(self):
   #   # First, we need to figure out the wavelengths for the control points
   #   self.knots = self.parent.ave_waves
   #   self.factors = self.knots*0 + 1.0

   def init_pars(self):
      # one parameter for each knot points
      p = []
      pi = [{'fixed':0, 'limited':[1,0], 'limits':[0.0,0.0], 'step':0.001, 'value':1.0} \
            for i in range(len(self.parent.allbands))]
      bs = self.parent.allbands
      nf = len(bs)
      knots = self.parent.ave_waves
      if bs[0] == 'blue':
         for i in range(1, len(self.parent.fluxes) + 1):
            pi[i]['value'] = self.parent.fluxes[i-1]
         if self.gradient:
            dw = (knots[0] - knots[1])/(knots[1] - knots[2])
            pi[0]['tied']='p[1] + %8.4f*(p[1]-p[2])' % (dw)
            pi[0]['value'] = self.parent.fluxes[0]
         else:
            pi[0]['tied'] = 'p[1]'
      else:
         for i in range(0, len(self.parent.fluxes)):
            pi[i]['value'] = self.parent.fluxes[i-1]

      if self.parent.allbands[-1] == 'red':
         pi[-1]['value'] = pi[-2]['value']
         if self.gradient:
            dw = (knots[-1] - knots[-2]) / (knots[-2] - knots[-3] )
            pi[-1]['tied']='p[%d] + %8.4f*(p[%d]-p[%d])' % (nf-2, dw, nf-2, nf-3)
         else:
            pi[-1]['tied'] = 'p[%d]' % (nf - 2)
      return(pi)

   def set_pars(self, pars):
      self.pars = num.asarray(pars)

   def __call__(self, x):
      if self.parent.ave_waves is None or self.pars is None:
         return x*0+1.0
      if self.gradient:
         xyds = tspline.tspsi(self.parent.ave_waves, self.pars, 
               tension=self.tension, curvs=[0,0])
      else:
         xyds = tspline.tspsi(self.parent.ave_waves, self.pars, 
               tension=self.tension, slopes=self.slopes)

      res = tspline.tsval1(x, xyds)

      y0 = tspline.tsval1(self.parent.ave_waves[0], xyds)[0]
      y1 = tspline.tsval1(self.parent.ave_waves[-1], xyds)[0]
      if self.gradient:
         dy0 = tspline.tsval1(self.parent.ave_waves[0], xyds, 1)[0]
         dy1 = tspline.tsval1(self.parent.ave_waves[-1], xyds, 1)[0]
         res = num.where(x < self.parent.ave_waves[0], 
               y0 + (x-self.parent.ave_waves[0])*dy0, res)
         res = num.where(x > self.parent.ave_waves[-1], 
               y1 + (x - self.parent.ave_waves[-1])*dy1, res)
      else:
         res = num.where(x < self.parent.ave_waves[0], y0, res)
         res = num.where(x > self.parent.ave_waves[-1], y1, res)
      return(res)

class f_ccm(function):

   def __init__(self, parent, tension=None, gradient=False, slopes=None, verbose=False):
      self.pars = [1.0, 0.0, 3.1]
      self.verbose = verbose

   def init_pars(self):
      # one parameter for each of scale, ebv, rv
      pi = [{'fixed':0, 'limited':[1,0], 'limits':[0.0, 0.0], 'value':1.0, 'step':0.0001},
            {'fixed':0, 'limited':[0,0], 'value':0.0, 'step':0.0001},
            {'fixed':0, 'limited':[1,0], 'limits':[0.0, 0.0], 'value':3.1, 'step':0.0001}]
      return(pi)

   def set_pars(self, pars):
      self.pars = pars

   def __call__(self, x):
      m,a,b = deredden.unred(x, x*0+1,-self.pars[1], R_V=self.pars[2])
      m *= self.pars[0]
      return m

mangle_functions = {'tspline':f_tspline,
                    'ccm':f_ccm}


class mangler:
   '''Given a spectrum and spectrum object, find the best set of a function's 
   paramters, such that the function multiplied by the spectrum produces
   the colors specified.'''

   def __init__(self, wave, flux, method, z=0, **margs):

      self.wave = wave
      self.flux = flux
      self.ave_waves = None
      if method not in mangle_functions:
         methods = ",".joing(mangle_funcitons.keys())
         raise ValueError, "method must be one of the following:\n%s" % methods
      self.function = mangle_functions[method](self, **margs)
      self.z = z
      self.verbose=margs.get('verbose', False)

   def get_colors(self, bands):
      '''Given a set of filters, determine the colors of the mangled
      spectrum for the current set of function parameters.  You'll get
      N-1 colors for N bands and each color will be bands[i+1]-bands[i]'''
      return self.get_mags(bands[:-1]) - self.get_mags(bands[1:])

   def get_mags(self, bands):
      '''Given a set of filters, determine the mags in bands of the mangled
      spectrum for the current set of function parameters.'''
      mflux = self.function(self.wave)*self.flux
      mags = []
      for b in bands:
         mags.append(fset[b].synth_mag(self.wave, mflux, z=self.z))
      return(num.array(mags))

   def get_mflux(self):
      '''Given the current paramters, return the mangled flux.'''
      if self.verbose:  print 'calling mangling function'
      mflux = self.function(self.wave)
      if self.verbose:  print 'done'
      return(self.flux*mflux)

   def create_anchor_filters(self, bands, anchorwidth):
      '''Given a set of filters in bands, create two fake filters at
      either end to serve as anchor filters.'''
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
   
      if wave0 < self.wave[0]*(1+self.z) or wave1 > self.wave[-1]*(1+self.z): 
         print 'Problem in mangle_spectrum: SED does not cover anchor filter '+\
               'definitions'
         if(self.verbose): 
            print 'Anchor filter definitions cover  ',wave0, 'A to ',wave1,'A'
            print 'SED covers ',self.wave[0]*(1+self.z),'A to ',\
                  self.wave[-1]*(1+self.z),'A (observed)'
            print "Try reducing ANCHORWIDTH or a mangling function " +\
                  "that doesn't need anchors."
      if(self.verbose): print 'Anchor filter definitions cover  ', \
                               wave0,'A to ',wave1,'A'
      

   def solve(self, bands, colors, norm_filter=None, fixed_filters=None, 
         anchorwidth=100, xtol=1e-10, ftol=1e-10, gtol=1e-10):
      '''Solve for the mangling function that will produce the observed colors
      in the filters defined by bands.'''

      self.bands = bands
      if len(bands) != len(colors) + 1:
         raise ValueError, "length of bands must be one more than colors"

      if norm_filter is None:
         norm_filter = bands[-1]
      else:
         if norm_filter not in bands:
            raise ValueError, "norm_filter must be one of bands"
      colors = num.asarray(colors)
      if fixed_filters is not None:
         if fixed_filters == "blue":
            fixed_filters = [bands[0]]
         elif fixed_filters == 'red':
            fixed_filters = [bands[-1]]
         elif fixed_fitlers == 'both':
            fixed_filters = [bands[0], bands[-1]]
         else:
            raise ValueError, "fixed_filters must be 'blue','red', or 'both'"
         self.allbands = bands
      else:
         self.create_anchor_filters(bands, anchorwidth)
         self.allbands = ['blue'] + bands + ['red']

      # The average wavelengths of the filters being used
      self.ave_waves = num.array([fset[b].ave_wave for b in self.allbands])

      # Construct the flux levels we want from the colors
      flux_rats = num.power(10, 0.4*colors)    # N-1 of these
      fluxes = num.zeros((len(bands),), dtype=num.Float64)
      fluxes[0] = fset[bands[0]].response(self.wave, self.flux, z=self.z)/\
            num.power(10,0.4*fset[bands[0]].zp)
      for i in range(1, len(bands)):
         fluxes[i] = fluxes[i-1]*flux_rats[i-1]
      id = bands.index(norm_filter)
      factor = fset[norm_filter].response(self.wave, self.flux, z=self.z)/\
            num.power(10, 0.4*fset[norm_filter].zp)/fluxes[id]

      self.fluxes = fluxes*factor
      if self.verbose:
         print "You input the following colors:"
         print colors
         print "This translates to the following fluxes:"
         print self.fluxes
         print "Which, for self-consistency, should be:"
         print [-2.5*num.log10(self.fluxes[i]/self.fluxes[i+1]) for i in range(len(bands)-1)]
         print "Initial colors for this spectrum are:"
         print self.get_colors(bands)

      # Do some initial setup
      pi = self.function.init_pars()

      if self.verbose:
         quiet = 0
      else:
         quiet = 1
      result = mpfit.mpfit(self.leastsq, parinfo=pi, quiet=quiet, maxiter=200,
            ftol=ftol, gtol=gtol, xtol=xtol, functkw={'bands':bands})
      if (result.status == 5) : print \
        'Maximum number of iterations exceeded in mangle_spectrum'
      self.function.set_pars(result.params)
      if self.verbose:
         print "The final colors of the SED are:"
         print self.get_colors(bands)
      return(result)

   def leastsq(self, p, fjac, bands):
      self.function.set_pars(p)
      mflux = self.get_mflux()
      mfluxes = num.array([fset[b].response(self.wave, mflux, z=self.z)/\
            num.power(10, 0.4*fset[b].zp) for b in self.bands])
      return (0, self.fluxes - mfluxes)

messages = ['Bad input parameters','chi-square less than ftol',
      'paramters changed less than xtol',
      'chi-square less than ftol *and* paramters changed less than xtol',
      'cosine between fvec and jacobian less than gtol',
      'Exceeded maximum number of iterations',
      'ftol is too small','xtol is too small','gtol is too small']

def mangle_spectrum2(wave,flux,bands, colors, fixed_filters=None, 
      normfilter=None, z=0, verbose=0, anchorwidth=100,
      method='tspline', xtol=1e-10, ftol=1e-10, gtol=1e-10, **margs):
   m = mangler(wave, flux, method, z=z, verbose=verbose, **margs)
   res = m.solve(bands, colors, norm_filter=normfilter, fixed_filters=fixed_filters,
         anchorwidth=anchorwidth, xtol=xtol, ftol=ftol, gtol=gtol)
   if res.status > 4:
      print "Warning:  %s" % messages[res.status]
   elif res.status < 0:
      print "Warning:  some unknown error occurred"
   elif verbose:
      print "mpfit finised with:  %s" % messages[res.status]

   return (m.get_mflux(), m.ave_waves, m.function.pars)

def apply_mangle(wave,flux,sw,sf,method='tspline', **margs):

   m = mangler(wave, flux, method, **margs)
   m.ave_waves = sw
   m.function.set_pars(sf)
   return(m.get_mflux())
