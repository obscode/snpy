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

   def __init__(self, parent, tension=None, gradient=False, slopes=None):
      function.__init__(self,parent)
      self.knots = None
      self.factors = None
      self.tension = tension       # specify a global tension parameter
      self.gradient = gradient     # constrain the anchors by slope?
      self.slopes = None           # End slopes if not using gradient
      self.pars = None

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

   def __init__(self, parent, tension=None, gradient=False, slopes=None):
      self.pars = [1.0, 0.0, 3.1]

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

   def __init__(self, wave, flux, method, z=0, verbose=0, **margs):

      self.wave = wave
      self.flux = flux
      self.ave_waves = None
      if method not in mangle_functions:
         methods = ",".joing(mangle_funcitons.keys())
         raise ValueError, "method must be one of the following:\n%s" % methods
      self.function = mangle_functions[method](self, **margs)
      self.z = z
      self.verbose=verbose

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
      mflux = self.function(self.wave)
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


##NAME
## colour_min2
##PURPOSE
## Thing to be minimized to get colours right.  Returns the fluxes of
## the spectrum scaled by scale_factors.  
##COMMON BLOCKS
## filter_info            Not modified.  Zero points must be set. 
##                         See calc_zps
## templates_colour_fit   Internal common block from mangle_spectrum.
##                         Not modified
##PROCEDURES CALLED
## apply_mangle
## filter_integ
##AUTHOR
## Mark Sullivan, August 2004
##-
#def colour_min2(scale_factors, fjac, scale_factors_wave, tobe_mangled_wave,
#      tobe_mangled_flux, method, rv, allfilters, z, wanted, useccm, sigma):
#
#   # This differs from colour_min in 3 respects: First, it allows the
#   # AB photometric system to be used.  Second, startnorm/finalnorm are not
#   # normalized by any zeropoint in this version.  Finally, it allows
#   # for an arbitrary redshift -- aconley 2005/06/16
#
#   global filts
#
#   #generate mangled spectrum
#   spectrum_mangle_flux = apply_mangle(tobe_mangled_wave,tobe_mangled_flux,
#                                    scale_factors_wave,scale_factors,
#                                    method, Rv=rv)
#
#   # calculate the new fluxes in each filter
#   if(useccm):
#      temp_flux = num.zeros((len(allfilters),), num.Float64)
#      for i in range(len(allfilters)):
#         temp_flux[i]=filts[allfilters[i]].response(tobe_mangled_wave,
#            spectrum_mangle_flux, z)
#         temp_flux[i] /= num.power(10.0, 0.4*filts[allfilters[i]].zp)
#   else:
#      temp_flux = num.zeros((len(allfilters)-2,), num.Float64)
#      for i in range(len(allfilters)-2):
#         temp_flux[i]=filts[allfilters[i+1]].response(tobe_mangled_wave,
#            spectrum_mangle_flux, z)
#         temp_flux[i] /= num.power(10.0, 0.4*filts[allfilters[i+1]].zp)
#      
#   return (0, (temp_flux - wanted)/sigma)
#   
##NAME
## apply_mangle
##PURPOSE
## Apply the mangling to a spectrum given pretabulated coefficients,
##  which are determined by mangle_spectrum2
##INPUTS
## tobe_mangled_wave      Input wavelength vector
## tobe_mangled_flux      Input flux vector
##OPTIONAL INPUTS
## sw                     Wavelengthts of knots (not used in CCM)
## sf                     Scale factors of spline (E(B-V) for CCM,
##                        first element only)
## Rv                     Rv to use (CCM only; def: 3.1)
## ebmv                   alternative way to pass ebmv for CCM
##KEYWORDS
## spline              Use spline interpolation (DEFAULT)
## linear              Use linear interpolation (uses linterp)
## poly                Use polynomial interpolation (Not yet implemented)
## ccm                 Use the CCM dust law
##RETURNS
## Mangled spectrum
##COMMON BLOCKS
##  None
##AUTHOR
## Mark Sullivan
##-
#def apply_mangle(tobe_mangled_wave,tobe_mangled_flux,sw,sf,method='spline',
#      verbose=0,Rv=3.1):
#
#   ccm = linear = spline = poly = tsp = 0
#   if method=='spline':  spline = 1
#   if method=='tspline':  tsp = 1
#   if method=='linear':  linear = 1
#   if method=='poly':  poly = 1
#   if method=='ccm':  ccm = 1
#   
#   if(verbose):
#      if(spline): print 'Using spline interpolation'
#      if(tsp): print 'Using tension spline interpolation'
#      if(linear): print 'Using linear interpolation'
#      if(poly): print 'Using polynomial interpolation'
#      if(ccm): print 'Using CCM dust law'
#   
#   if(not ccm):
#      nscale_factors=len(sf)
#      nwavelengths=len(tobe_mangled_wave)
#      
#      # extrapolate the scalefactors in the ends of the spectrum using
#      # linear extrapolation
#      start_scale=sf[0]
#      end_scale=sf[-1]
#      
#      manglefactor=num.zeros((nwavelengths,), num.Float64)
#      ispline_start = num.nonzero(num.greater(tobe_mangled_wave-sw[0], 0))
#      if len(ispline_start) == 0:
#         ispline_start = 0
#      else:
#         ispline_start = ispline_start[0]
#         
#      ispline_end = num.nonzero(num.greater(tobe_mangled_wave-sw[-1], 0))
#      if len(ispline_end) == 0:
#         ispline_end = nwavelengths - 1
#      else:
#         ispline_end = ispline_end[0]
#
#      if (ispline_start > 0) : manglefactor[0:ispline_start]=start_scale
#      if (ispline_end < (nwavelengths-1)): 
#         manglefactor[ispline_end:]=end_scale
#      
#      # do the mangling to the passed spectrum
#      
#      if(spline):
#         tck = scipy.interpolate.splrep(sw, sf, k=3, s=0)
#         res = scipy.interpolate.splev(tobe_mangled_wave[ispline_start:ispline_end], tck)
#         manglefactor[ispline_start:ispline_end] = res
#         # Now linearly extrapolate the rest:
#         s1 = scipy.interpolate.splev(tobe_mangled_wave[ispline_start+1], tck, 1)
#         s2 = scipy.interpolate.splev(tobe_mangled_wave[ispline_end-1], tck, 1)
#         manglefactor[ispline_end:] = manglefactor[ispline_end] + \
#               (tobe_mangled_wave[ispline_end:] - tobe_mangled_wave[ispline_end])*s1
#         manglefactor[0:ispline_start] = manglefactor[ispline_start] + \
#               (tobe_mangled_wave[0:ispline_start] - \
#               tobe_mangled_wave[ispline_start])*s2
#      elif (linear):
#         tck = scipy.interpolate.splrep(sw, sf, k=1, s=0)
#         res = scipy.interpolate.splev(tobe_mangled_wave[ispline_start:ispline_end], tck)
#         manglefactor[ispline_start:ispline_end] =  res
#      elif (poly) :
#         raise RuntimeError,"Polynomial mangling not supported"
#
#      elif (tsp):
#         xyds = tspline.tspsi(sw, sf)
#         res = tspline.tsval1(tobe_mangled_wave[ispline_start:ispline_end], xyds)
#         manglefactor[ispline_start:ispline_end] = res
#         s1 = tspline.tsval1(tobe_mangled_wave[ispline_start+1], xyds, 1)
#         s2 = tspline.tsval1(tobe_mangled_wave[ispline_end-1], xyds, 1)
#         manglefactor[ispline_end:] = manglefactor[ispline_end-1] + \
#               (tobe_mangled_wave[ispline_end:] - \
#               tobe_mangled_wave[ispline_end-1])*s1
#         manglefactor[0:ispline_start] = manglefactor[ispline_start] + \
#               (tobe_mangled_wave[0:ispline_start] - \
#               tobe_mangled_wave[ispline_start])*s2
#
#      spectrum_mangle_flux=tobe_mangled_flux*manglefactor
#   else:
#      #;Using CCM law
#      ebmv=sf[0]
#      scalefac = sf[1]
#      spectrum_mangle_flux,a,b = deredden.unred(tobe_mangled_wave,
#            tobe_mangled_flux,-ebmv, R_V=Rv)
#      spectrum_mangle_flux *= scalefac
#   
#   return(spectrum_mangle_flux)
#
##NAME
## mangle_spectrum2
##PURPOSE
## Mangles an input spectrum so that it has the desired colors
## This version avoids some of the common blocks
##USAGE
## mang = mangle_spectrum( spectrum_wave, spectrum_flux, passedfilters,
##                         wanted_colours, vega_struct )
##INPUTS
## spectrum_wave       Wavelength vector for input spectrum
## spectrum_flux       Flux vector for input spectrum
## passedfilters       Filter numbers.
## wanted_colours      Desired colors.  Corresponding to [0]-[1],
##                     [1]-[2], [2]-[3] where [] refers to the entries
##                     of passedfilters
## vega_struct         Structure holding info about vega.  See read_vega
##OPTIONAL INPUTS
## redshift            What redshift to  this at.  Default: 0.0.
##                      Note that the returned SED is still in the rest
##                      frame, as the filters are shifted instead of
##                      the SED.
## normfilternum       Index into passedfilters for the filter used to
##                      keep the normalization of the spectrum.  After
##                      mangling, the flux through this filter should
##                      be the same as before mangling.  Default: 0 --
##                      that is, the first filter in passedfilter is
##                      used to normalize
## fixedfilters        Numbers of anchor filters.  if not provided, the
##                      code will choose its own anchor filters
## photsys             Photometric system colors are defined in.  1 for
##                      Vega, 2 for AB.  Default: 1
## truncate_filters     This adjusts the filter response to neglect all
##                     parts of the filter where the throughput is
##                     <truncate_pc%. Some filters have very wide definitions
##                     even where there is essentially no throughput. Not
##                     recommed for light-curve fitting, but great
##                     for mangling real-life spectra! Default: NO
## truncate_pc         Default 1: percentage throughput at whihc to
##                     truncate the filter responses
## anchorwidth         Width of the automatic anchorfilters. Default: 100A
## Rv                  Rv to use (CCM only; def: 3.1)
## noublesample      Array of bytes specifying whether to not use
##                      ublesampling in filter_integ.  The default is
##                      all 0b's -- i.e., to use ublesampling.
##OPTIONAL OUTPUTS
## factors_w           The wavelengths of the spline knots.  Note that
##                      these will be in the rest frame of the SED you
##                      pass in, so if you set the redshift parameter
##                      these may be divided by a factor of 1+z
##                      you n't expect. Not used for CCM.
## factors_s           The scale factors of the applied function. for
##                     CCM the first element is the best-fit E(B-V)
## ebmv                Alternative way of getting best E(B-V); CCM only
##KEYWORDS
## spline              Use spline interpolation (DEFAULT)
## linear              Use linear interpolation (uses linterp)
## poly                Use polynomial interpolation (Not yet implemented)
## ccm                 Use the CCM dust law
## quiet               Run MPFIT quietly. Default:Yes. You must specify
##                     QUIET=0 to override this. (extreme debug).
## verbose             Verbose output
## usemeanwave         Use the mean wavelengths of the filters for the
##                      knot points rather than the effective
##                      wavelength
## gradient            Use the gradient to set the anchor boundary
##                      conditions instead of fixed levels
##RETURNS
## The output spectrum, normalized so that the flux in the first filter
## is the same before and after mangling
##NOTES
## if factors_w and factors_s are known ahead of time, it is far faster
##  to use apply_mangle.
##COMMON BLOCKS
## filter_info         Note this is modified internally by adding the
##                      anchor filters to the .  These are removed
##                      at the  of this code
## templates_color_fit Created by this routine for passing to
##                      colour_min2, the minimizer function
##PROCEDURES CALLED
## filter_integ
## colour_min2
## effective_wavelength
## calc_ab_zp
## calc_vega_zp
## linterp
##AUTHOR
## Mark Sullivan, August 2004
##-
#
##;Little function to calculate colour of spectrum through two filters
##;Kind of like kcorr, but both filters are integrated at the object
##;redshift.  Returns the color in magnitudes
#
#def mangle_spectrum2_getcolour(filt1,filt2,wave,flux,redshift):
#   global filts
#
#   f1 = filts[filt1]
#   f2 = filts[filt2]
#
#   flux1 = f1.response( wave, flux, redshift)
#   flux2 = f2.response( wave, flux, redshift)
#   
#   return(-2.5*num.log10(flux1/flux2) + f1.zp - f2.zp)
#
#
#def mangle_spectrum2(s_wave,s_flux,passedfilters, wanted_colours,
#      fixedfilters=None, normfilter=None,z=0,verbose=0,
#      factors_s=None, factors_w=None, gradient=0, usemeanwave=1, 
#      truncate_filters=0, truncate_pc=1.0, quiet=1, anchorwidth=100,
#      method='spline', Rv=3.1):
#   global filts
#
#   if truncate_filters:
#      filts = deepcopy(fset)
#
#   usespline=usetspline=usepoly=uselinear=useccm = 0
#   if method == 'spline':  usespline = 1
#   if method == 'tspline':  
#      if tspline is not None:
#         usetspline = 1
#      else:
#         print "Warning:  tension spline module not installed, using regular" +\
#               " spline."
#         usespline = 1
#         method = 'spline'
#   if method == 'linear':  uselinear = 1
#   if method == 'poly':  usepoly = 1
#   if method == 'ccm':  useccm = 1
#   if verbose:
#      print ''
#      print 'MANGLE_SPECTRUM2 starting.'
#      if(usespline): print 'Using spline interpolation'
#      if(usetspline): print 'Using tension spline interpolation'
#      if(uselinear): print 'Using linear interpolation'
#      if(usepoly): print 'Using polynomial interpolation'
#      if(useccm): print 'Using CCM dust law'
#   
#   # get into percentage
#   truncate_pc=truncate_pc/100.
#   
#   obj_redshift = z
#   opz = 1.0 + obj_redshift
#   
#   #;Set which filter will be used to normalize
#   if normfilter is None or normfilter not in passedfilters:
#      normfilter = passedfilters[-1]
#   
#   # check spectrum is OK
#   if len(s_wave) != len(s_flux):
#      raise RuntimeError, 'spectrum wave and flux must have same number' + \
#                            ' of elements in mangle_spectrum !'
#   if num.maximum.reduce(s_flux) <= 0.0:
#      raise RuntimeError, 'ERROR in mangle_spectrum: Your spectrum must ' + \
#                          ' contain some positive (non-zero) fluxes.'
#   
#   if len(passedfilters) != len(wanted_colours)+1:
#      raise RuntimeError, 'Number of colours passed must be equal to number '+ \
#         'filters - 1 in mangle_spectrum !'
#   
#   # truncate_filters assumes filter response are normalised to a peak
#   # throughput of 1 !
#   if(truncate_filters): 
#      if(verbose): print ''
#      if(verbose): print 'Truncated passed filters:'
#      for i in range(len(passedfilters)):
#   
#                                   # find peak throughput
#         #filts[passedfilters[i]] = filters[passedfilters[i]].copy()
#         f = filts[passedfilters[i]]
#         f.resp = num.where(num.greater(f.resp/num.maximum.reduce(f.resp), 
#            truncate_pc), f.resp, 0.0)
#         ids = num.equal(f.resp, 0.0)
#         f.resp = f.resp[ids[0]-1:ids[-1]+1]
#         f.wave = f.wave[ids[0]-1:ids[-1]+1]
#         #f.min = num.minimum.reduce(f.wave)
#         #f.max = num.maximum.reduce(f.wave)
#         
#   # check that the (possibly truncated) passedfilters and spectrum completely 
#   #  overlap
#   low_waves = []
#   high_waves = []
#   for i in range(len(passedfilters)):
#       f = filts[passedfilters[i]]
#       low_waves.append(f.wavemin)
#       high_waves.append(f.wavemax)
#   
#   lowest_filter_wavelength = min(low_waves)
#   highest_filter_wavelength = max(high_waves)
#   
#   bluest_spectrum_wavelength=s_wave[0]
#   reddest_spectrum_wavelength=s_wave[-1]
#   
#   if(verbose): 
#      print ''
#      print 'Filter definitions cover  ',lowest_filter_wavelength,'A to ', \
#            highest_filter_wavelength,'A'
#      print 'SED covers ',bluest_spectrum_wavelength,'A to ',\
#            reddest_spectrum_wavelength,'A'
#   
#   if(lowest_filter_wavelength < bluest_spectrum_wavelength or 
#      highest_filter_wavelength > reddest_spectrum_wavelength): 
#      print 'Problem in mangle_spectrum: SED does not cover filter definitions'
#      if(verbose): 
#         print 'Filter definitions cover  ',lowest_filter_wavelength,\
#               'A to ',highest_filter_wavelength,'A'
#         print 'SED covers ',bluest_spectrum_wavelength,'A to ', \
#               reddest_spectrum_wavelength,'A'
#   
#   # copy the passed arrays to the common block
#   #tobe_mangled_flux=spectrum_flux
#   #tobe_mangled_wave=spectrum_wave
#   
#   nwavelengths=len(s_wave)
#   
#   # fixed filters contain extra spline knots - these are the "anchor" filters
#   # if these are not specified, we must work out what they should be
#   
#   # make an array containg all the filters
#   internalfilters=0
#   allfilters=passedfilters[:]
#   
#   #;Choose the fixedfilters if not set by the caller
#   if(not fixedfilters and  not useccm): 
#       internalfilters=1
#       if(verbose) : print \
#         'No fixed filters passed - defining internally..'
#       nfilters=len(passedfilters)
#       mean_wave=num.zeros((nfilters,), num.Float64)
#       for i in range(len(passedfilters)):
#           pf = filts[passedfilters[i]]
#           if ( usemeanwave ) :  
#              mean_wave[i]=pf.ave_wave / opz
#           else: 
#              mean_wave[i]=pf.eff_wave(s_wave,s_flux,z) / opz
#       mean_wave = num.array(mean_wave)
#       #dummy=MIN(mean_wave,bluest_filter,SUBSCRIPT_MAX=reddest_filter)
#       reddest_filter = passedfilters[num.argmax(mean_wave)]
#       bluest_filter = passedfilters[num.argmin(mean_wave)]
#   
#       #extract_filter,passedfilters[bluest_filter],filter_wave,dummy
#       bluest_filter_wave=filts[bluest_filter].wave[0]
#       #extract_filter,passedfilters[reddest_filter],filter_wave,dummy
#       reddest_filter_wave=filts[reddest_filter].wave[-1]
#   
#       # Create two new fake filters
#       filts['blue'] = filter('blue_anchor')
#       filts['red'] = filter('red_anchor')
#   
#       filts['blue'].wave = num.array([bluest_filter_wave-anchorwidth-2.,
#             bluest_filter_wave-anchorwidth-1., bluest_filter_wave-anchorwidth, 
#             bluest_filter_wave-3., bluest_filter_wave-2.,bluest_filter_wave-1.])
#       filts['blue'].resp = num.array([0.,0.,1.,1.,0.,0.])
#   
#       #;red anchor filter
#       filts['red'].wave = num.array([reddest_filter_wave+1.,
#            reddest_filter_wave+2., reddest_filter_wave+3.,
#            reddest_filter_wave+anchorwidth, reddest_filter_wave+anchorwidth+1.,
#            reddest_filter_wave+anchorwidth+2.])
#       filts['red'].resp = num.array([0.,0.0,1.,1.,0.0,0.])
#   
#       filts['red'].zp = filts['red'].compute_zpt(vegaB, 0.0)
#       filts['blue'].zp = filts['blue'].compute_zpt(vegaB, 0.0)
#   
#       filts['red'].ave_wave = num.sum(filts['red'].wave)/len(filts['red'].wave)
#       filts['blue'].ave_wave = num.sum(filts['blue'].wave)/len(filts['blue'].wave)
#   
#       lowest_filter_wavelength=bluest_filter_wave-anchorwidth-2.
#       highest_filter_wavelength=reddest_filter_wave+anchorwidth+2.
#   
#       if(lowest_filter_wavelength < s_wave[0] or \
#          highest_filter_wavelength > s_wave[-1]): 
#          print 'Problem in mangle_spectrum: SED does not cover anchor filter '+\
#                'definitions'
#          if(verbose): 
#             print 'Anchor filter definitions cover  ',lowest_filter_wavelength,\
#                   'A to ',highest_filter_wavelength,'A'
#             print 'SED covers ',bluest_spectrum_wavelength,'A to ',\
#                   reddest_spectrum_wavelength,'A'
#             print "Try reducing ANCHORWIDTH. But you're probably screwed."
#       if(verbose): print 'Anchor filter definitions cover  ',\
#             lowest_filter_wavelength,'A to ',highest_filter_wavelength,'A'
#   
#   if(not useccm): 
#       #; add the fixedfilters to the allfilters array if we're not using CCM
#       allfilters=['blue'] + allfilters + ['red']
#   
#   nfilters=len(allfilters)
#   
#   # calculate effective wavelengths of spectrum for the passed filters
#   # these are the positions of the spline knots
#   mean_wave=num.zeros((nfilters), num.Float64)
#   for i in range(nfilters):
#       if ( usemeanwave ) : 
#           mean_wave[i] = filts[allfilters[i]].ave_wave / opz
#       else: 
#           mean_wave[i] = filts[allfilters[i]].eff_wave( s_wave,s_flux,z ) / opz
#   if(not useccm): 
#      passed_mean_wave=mean_wave[1:-1]
#   
#      if (mean_wave[0] >= min(passed_mean_wave)): 
#         raise RuntimeError,"Your lower fixed filter must be bluer than the bluest filter passsed!"
#      if (mean_wave[nfilters-1] <= max(passed_mean_wave)): 
#         raise RuntimeError,"Your upper fixed filter must be redder than the reddest filter passsed!"
#   else: 
#      passed_mean_wave=mean_wave
#   
#   
#   # B-band flux of tobe_mangled_wave/tobe_mangled_flux
#   # Note that we don't bother with the zeropoint here
#   startnorm = filts[normfilter].response( s_wave, s_flux, z)
#   if(verbose): 
#      #; print wanted colours
#      print 'The colours you want are:'
#      for i in range(len(wanted_colours)):
#         print "%s (%10.3f) minus %s (%10.3f) = %8.4f" % (passedfilters[i],\
#               passed_mean_wave[i],passedfilters[i+1],\
#               passed_mean_wave[i+1],wanted_colours[i])
#   
#   if(verbose): 
#   # print colours of spectrum passed
#      spectrum_colour=num.zeros((len(wanted_colours)), num.Float64)
#      print 'The colours of the spectrum you passed are:'
#      for i in range(len(wanted_colours)):
#         #;Kcorr no longer works here because we are generalizing
#         #; to other redshifts, and we have to integrate both through
#         #; the same filter
#         spectrum_colour[i] = \
#            mangle_spectrum2_getcolour( passedfilters[i],passedfilters[i+1],\
#                                        s_wave, s_flux, obj_redshift)
#   
#         print "%s (%10.3f) minus %s (%10.3f) = %8.4f" % (passedfilters[i],
#               passed_mean_wave[i], passedfilters[i+1],
#               passed_mean_wave[i+1],spectrum_colour[i])
#   
#   if(not useccm): 
#      # setup initial scale factors array same size as number of passed filters
#      # set to 1.0 at start of fit
#      scale_factors_or=num.ones((nfilters,), num.Float64)
#      nscale_factors=nfilters
#   
#      # setup an array of "wanted" fluxes
#      wantedflux_at_mean_wave=num.zeros((nfilters-2),num.Float64)
#   else: 
#      scale_factors_or=num.array([ 0.0, 1.0 ]) #;ebmv, fluxscale
#      nscale_factors=2
#      wantedflux_at_mean_wave=num.zeros((nfilters),num.Float64)
#   
#   
#   wantedflux_at_mean_wave[0]=filts[passedfilters[0]].response(s_wave,
#                                            s_flux,z)
#   wantedflux_at_mean_wave[0] /= num.power(10, 0.4*filts[passedfilters[0]].zp)
#   
#   # calculate relative flux in other filters to first passed filter
#   # using passed colours
#   for i in range(len(wanted_colours)):
#      wantedflux_at_mean_wave[i+1]=  num.power(10.0,0.4*wanted_colours[i]) * \
#                                     wantedflux_at_mean_wave[i]
#   
#   #;Get normfilter to have the same flux before and after
#   desired_flux_in_normfilter = startnorm / \
#         num.power(10.0, 0.4*filts[normfilter].zp)
#   
#   factor = desired_flux_in_normfilter / \
#            wantedflux_at_mean_wave[ passedfilters.index(normfilter) ]
#   wantedflux_at_mean_wave = wantedflux_at_mean_wave * factor
#   
#   # setup our fitting weights including the two  values used
#   # for constraining the spline
#   # sigma is not used but I think must be specified
#   #sigma = num.ones((len(mean_wave)),num.Float32)*wantedflux_at_mean_wave[0]/1000.0
#   sigma = wantedflux_at_mean_wave/1000.0
#   # weights = num.ones((len(mean_wave)),num.Float32)*1000.0
#   
#   # sort the wavelengths of the passed filters
#   # these are set as the wavelengths of the scale factors
#   sort_index=num.argsort(mean_wave)
#   scale_factors_wave=num.take(mean_wave, sort_index)
#   #sigma = num.take(sigma, sort_index)
#   #weights = num.take(weights, sort_index)
#   
#   if(useccm): 
#      pi = [{'parname':'E(B-V)', 'fixed':0, 'limited':[0,0], 
#                  'limits':[0.0,0.0]},
#                 {'parname':'scale', 'fixed':0, 'limited':[0,0], 
#                  'limits':[0.0,0.0]}]
#      pi[0]['value'] = 0.0
#      pi[1]['value'] = 1.0
#   else: 
#      #; Set the boundaries etc. of the fitted scales
#      #; nothing is limited except fluxscale cannot be negative
#      pi = [{'fixed':0, 'limited':[1,0], 'limits':[0.0,0.0], 'step':0.001} \
#            for i in range(nscale_factors)]
#   
#      if(gradient): 
#         #; use the gradient of the two filters nearest the edge to set the
#         #; fluxscale in the anchor filters
#         #; gradient = dy/dx i.e. d(fluxscale)/d(wavelength)
#         #; red anchor
#         wavediff=( mean_wave[-1] - mean_wave[-2] ) / \
#                  ( mean_wave[-2] - mean_wave[-3] )
#         pi[nfilters-1]['tied']='p[%d] + %8.4f*(p[%d]-p[%d])' % (nfilters-2,
#               wavediff, nfilters-2, nfilters-3)
#         #blue anchor
#         wavediff=(mean_wave[0]-mean_wave[1])/(mean_wave[1]-mean_wave[2])
#         pi[0]['tied']='p[1]+%8.4f*(p[1]-p[2])' % (wavediff)
#      else: 
#         # tie the scales at the fixed filters to the scale in the filter next or
#         pi[0]['tied']='p[1]'
#         pi[nfilters-1]['tied']='p[%d]' % (nfilters-2)
#   
#      for i in range(nscale_factors):
#         pi[i]['value'] = 1.0
#   
#   # FIT! Change XTOL to get better fits
#   status=0
#   if(verbose): print 'Calling MPFITFUN...'
#   result = mpfit.mpfit(colour_min2, parinfo=pi, quiet=quiet, maxiter=200,
#         ftol=1.e-10, gtol=1.0e-10, xtol=1.0e-10, functkw=\
#         {'scale_factors_wave':scale_factors_wave,
#          'tobe_mangled_wave':s_wave,
#          'tobe_mangled_flux':s_flux,
#          'method':method, 'rv':Rv, 'allfilters':allfilters, 'z':z,
#          'wanted':wantedflux_at_mean_wave, 'useccm':useccm, 'sigma':sigma})
#   if (result.status == 5) : print \
#     'Maximum number of iterations exceeded in mangle_spectrum'
#   scale_factors = num.array(result.params)
#   
#   # generate final spectrum
#   final_mangled_flux = apply_mangle(s_wave,s_flux,
#                                       scale_factors_wave,scale_factors,
#                                       method, Rv=Rv)
#   
#   # scale spectrum so it has the same B-flux as that passed
#   # No zeropoint here because it wasn't used when startnorm was set
#   final_norm=filts[normfilter].response(s_wave, final_mangled_flux,z)
#   if final_norm > 0.0: 
#      factor=startnorm/final_norm 
#   else: 
#      factor=1.0
#   final_mangled_flux = final_mangled_flux * factor
#
#   
#   if verbose: 
#      if not useccm: 
#          #;We have two extra filters -- the 'guard' ones
#          temp_flux = num.zeros((nfilters-2,), typecode=num.Float64 )
#          init_flux = num.zeros((nfilters-2,), typecode=num.Float64 )
#      else: 
#          temp_flux = num.zeros((nfilters,), typecode=num.Float64 )
#          init_flux = num.zeros((nfilters,), typecode=num.Float64 )
#      print "%-8s   %12s   %12s   %12s" % ("Filter:","InitFlux","WantedFlux:",
#            "FinalFlux:")
#      for i in range(len(passedfilters)):
#          cfiltnum = passedfilters[i]
#          init_flux[i] = fset[cfiltnum].response(s_wave,s_flux, obj_redshift)
#          temp_flux[i] = fset[cfiltnum].response(s_wave,final_mangled_flux, obj_redshift)
#          temp_flux[i] /= num.power(10, 0.4*fset[cfiltnum].zp)
#          init_flux[i] /= num.power(10, 0.4*fset[cfiltnum].zp)
#
#          print "%-8s   %12.5g   %12.5g   %12.5g" % (cfiltnum,init_flux[i],
#                wantedflux_at_mean_wave[i],temp_flux[i])
#   if verbose: 
#   # and check we did it right
#      print 'The final colours of the mangled spectrum are:'
#      for i in range(len(wanted_colours)):
#         spectrum_colour[i] = \
#            mangle_spectrum2_getcolour( passedfilters[i],passedfilters[i+1],
#                                        s_wave, final_mangled_flux,
#                                        obj_redshift)
#   
#         if not useccm: 
#            print "%s (%10.3f) minus %s (%10.3f) = %8.4f" % (passedfilters[i],
#                  mean_wave[i+1], passedfilters[i+1], mean_wave[i+2],
#                  spectrum_colour[i])
#         else: 
#            print "%s (%10.3f) minus %s (%10.3f) = %8.4f" % (passedfilters[i],
#                  mean_wave[i],passedfilters[i+1], mean_wave[i],
#                  spectrum_colour[i])
#
#   if(useccm): ebmv=scale_factors[0]
#   if truncate_filters:
#      filts = fset
#   
#   return(final_mangled_flux, scale_factors_wave, scale_factors)
#   
