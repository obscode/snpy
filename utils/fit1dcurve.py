'''This module provices an abstract class whose purpose is fitting a 1D curve
using several different methods, but ultimately provide estimates of the 
interpolated points, derivatives, and locations of maxima/minima.'''

import numpy as num
from snpy.utils import fit_spline
from scipy.interpolate import splrep,splev,sproot
from scipy.optimize import brentq, newton
from scipy.misc import derivative as deriv
try:
   from snpy.spline2 import spline2, evalsp,eval_extrema,eval_x
except:
   spline2 = None
try:
   from numpy import polynomial
except:
   polynomial = None
try:
   import pymc
   from pymc import gp as GP
   import os
   if 'OMP_NUM_THREADS' not in os.environ:
      os.environ['OMP_NUM_THREADS'] = '1'
except:
   pymc = None

try:
   import InteractiveFit
except:
   InteractiveFit = None

functions = {}

class oneDcurve:

   num_real_keep = 100

   def __init__(self, x, y, ey):

      x = num.atleast_1d(x)
      y = num.atleast_1d(y)
      ey = num.atleast_1d(ey)

      if not len(x.shape) == 1:
         raise ValueError, "x, y, ey must be 1D arrays"
      if not x.shape == y.shape:
         raise ValueError, "x, y must have same shape"
      if not y.shape == ey.shape:
         raise ValueError, "y, ey must have same shape"

      self.xdata = x
      self.ydata = y
      self.eydata = ey
      self.vardata = num.power(ey,2)
      self.mask = num.ones(x.shape, dtype=bool)  # mask for the data

      self.realization = None
      self.realizations = []

      self.pars = {}
      self.setup = False
      self.ifit = None


   def __getattr__(self, key):
      if 'pars' in self.__dict__:
         if key in self.pars:
            return self.pars[key]
      if key == 'x':
         return self.xdata[self.mask]
      if key == 'y':
         return self.ydata[self.mask]
      if key == 'ey':
         return self.eydata[self.mask]
      if key == 'var':
         return self.vardata[self.mask]
      raise AttributeError, "Instance has not attribute %s" % (key)

   def __setattr__(self, key, value):
      if 'pars' in self.__dict__:
         if key in self.__dict__['pars']:
            self.__dict__['pars'][key] = value
            self.__dict__['setup'] = False
            if self.ifit is not None:
               self.ifit.redraw()
            return
         else:
            self.__dict__[key] = value
            return
      self.__dict__[key] = value

   def maskpoint(self, x, y):
      '''Mask the point closest to (x,y).'''
      id = num.argmin(num.power(self.x-x,2) + num.power(self.y-y,2))
      self.mask[id] = False
      self.setup = False

   def maskid(self, i):
      self.mask[id] = False
      self.setup = False

   def maskresids(self, absclip=None, sigclip=None):
      '''Mask data outside a range of residuals.  If [absclip]
      is specified, any residual > [absclip] is masked.  If
      [sigclip] is specified, any residual more than [sigclip]
      times the estimated standard deviation is masked.'''
      absdev = num.aboslute(self.residuals())
      if absclip is not None:
         self.mask = num.greater(absdev, absclip)
         self.setup = False
      elif sigclip is not None:
         sigma = 1.49*num.median(absdev)
         self.mask = num.greater(absdev, sigclip*sigma)
         self.setup = False

   def __call__(self, x):
      raise NotImplementedError('Derived class must overide')

   def error(self, x, N=50):
      '''Estimate the error at the point [x].'''
      scalar = (len(num.shape(x)) == 0)
      x = num.atleast_1d(x)
      
      if len(self.realizations) < N:
         for i in range(N-len(self.realizations)):
            self.draw()
         self.reset_mean()
      earr = []
      for i in range(N):
         self.realization = self.realizations[i]
         earr.append(self.__call__(x)[0])
      self.realization = None
      earr = num.array(earr)
      print earr.shape
      err = num.std(earr, axis=0)
      if scalar:
         return err[0]
      else:
         return err

   def draw(self):
      raise NotImplementedError('Derived class must overide')
   
   def reset_mean(self):
      raise NotImplementedError('Derived class must overide')

   def residuals(self, mask=True):
      '''returns an array of residuals.'''
      if mask:
         return self.y - self.__call__(self.x)[0]
      else:
         return self.ydata - self.__call__(self.xdata)[0]

   def rms(self):
      '''Returns RMS of residuals.'''
      return num.sqrt(num.mean(num.power(self.residuals(),2)))

   def chisquare(self):
      '''Returns the chi-square statistic for current parameters.'''
      return num.sum(num.power(self.residuals(),2)*num.power(self.var,-1))

   def rchisquare(self):
      raise NotImplementedError('Derived class must overide')

   def DW(self):
      '''Returns the Durvin-Watson statistic.'''
      r = self.residuals()
      return num.sum(num.power(r[1:] - r[:-1],2))/num.sum(num.power(r,2))

   def deriv(self, x, n=1):
      '''Returns the nth derivative of the function at x.'''
      raise NotImplementedError('Derived class must overide')

   def find_extrema(self):
      '''Find the position and values of the maxima/minima.'''
      raise NotImplementedError('Derived class must overide')

   def intercept(self, y):
      '''Find the value of x for which the interpolator goes through [y]'''
      raise NotImplementedError('Derived class must overide')

   def domain(self):
      '''Return the valid domain for this model'''

   def interact(self):
      '''If we have the InteractiveFit module, spawn an interactive fitter.'''
      if InteractiveFit is not None:
         return InteractiveFit.InteractiveFit(self)
      else:
         print "Sorry, you need to have matplotlib installed to use this feature"
         return None

   def help(self):
      '''Provide a help string.'''
      raise NotImplementedError('Derived class must overide')

polytypes = {'polynomial':polynomial.Polynomial,
             'chebyshev':polynomial.Chebyshev,
             'laguerre':polynomial.Laguerre,
             'hermite':polynomial.Hermite,
             'hermiteE':polynomial.HermiteE}

if polynomial is not None:
   class Polynomial(oneDcurve):
   
      def __init__(self, x, y, dy, **args):
         '''Fit an Nth order polynomial to the data.  The only arguments are
         [n], the order, [x0] the zero-point, xmin and xmax the lower and
         upper limits of the fit, respectively.'''
   
         oneDcurve.__init__(self, x, y, dy)
   
         self.pars = {
               'n':3,
               'type':'poly',
               'xmin':None,
               'xmax':None}
         for key in args:
            if key not in self.pars:
               raise TypeError, \
                     "%s is an invalid keyword argument for this method" % key
            self.pars[key] = args[key]
   
         if self.xmin is None:
            self.xmin = self.xdata.min()
         if self.xmax is None:
            self.xmax = self.xdata.max()
   
         # no need to regularize
         self._setup()
         self.realization = None

      def help(self):
         print 'The following parameters affect the fit:'
         print 'n:        order of the polynomial'
         print 'type:     type of polynomial:  polynomial, chebyshev, laguerre,'
         print '          hermite, or hermiteE'
         print 'xmin:     lower bound on data to interpolate'
         print 'xmax:     upper bound on data to interpolate'
   
      def __str__(self):
         return self.type + " polynomial"

      def _setup(self):
         '''Given the current set of params, setup the interpolator.'''
         self.mask = self.mask*num.greater_equal(self.xdata, self.xmin)*\
               num.less_equal(self.xdata, self.xmax)
   
         if self.type not in polytypes:
            raise ValueError, "Error:  the polynomial type must be one of " +\
                  ",".join(polytypes.keys())
         self.poly = polytypes[self.type].fit(self.x, self.y, deg=self.n, 
               w=num.power(self.ey,-1))
         self.setup = True
         self.realization = None
   
      def __call__(self, x):
         '''Interpolate at point [x].  Returns a 3-tuple: (y, mask) where [y]
         is the interpolated point, and [mask] is a boolean array with the same
         shape as [x] and is True where interpolated and False where 
         extrapolated'''
         if not self.setup:  self._setup()
         if self.realization is not None:
            res = self.realization(x)
         else:
            res = self.poly(x)
         return res, num.greater_equal(x, self.poly.domain[0])*\
               num.less_equal(x, self.poly.domain[1])
         
      def draw(self):
         '''Generate a random realization of the spline, based on the data.'''
         y_draw = num.random.normal(self.y, self.ey)
         self.realizations.append(\
               polytypes[self.type].fit(self.x, y_draw, deg=self.n,
               w=num.power(self.ey,-1)))
         if len(self.realizations) > self.num_real_keep:
            self.realizations = self.realizations[1:]
         self.realization = self.realizations[-1]


      def reset_mean(self):
         self.realization = None
   
      def rchisquare(self):
         chisq = self.chisquare()
         return chisq/(len(self.x) - 1 - len(self.poly.coef))
   
      def deriv(self, x, n=1):
         '''Returns the nth derivative of the function at x.'''
         dpoly = self.poly.deriv(m=n)
         return dpoly(x)
   
      def domain(self):
         '''Returns the valid domain of the polynomial.'''
         dom = self.poly.domain
         return (dom[0],dom[1])

      def find_extrema(self, xmin=None, xmax=None):
         '''Find the position and values of the maxima/minima.  Returns a tuple:
            (roots,vals,ypps) where roots are the x-values where the extrema
            occur, vals are the y-values at these points, and ypps are the
            2nd derivatives.  optionally, restrict roots to between xmin,
            and xmax'''
         if self.realization is not None:
            poly = self.realization
         else:
            poly = self.poly

         if xmin is None:  xmin = self.poly.domain[0]
         if xmax is None:  xmax = self.poly.domain[1]
         d1 = poly.deriv(m=1)
         d2 = poly.deriv(m=2)
         roots = d1.roots()
         # Roots can be complex.  Want only the real ones
         gids = num.iscomplex(roots)
         roots = num.real(roots[gids])
         gids = num.greater_equal(roots, xmin)*num.less_equal(roots, xmax)
         roots = roots[gids]
         if len(roots) == 0:
            return None,None,None
         vals = self.__call__(roots)
         curvs = d2(roots)
         curvs = num.where(curvs < 0, -1, curvs)
         curvs = num.where(curvs > 0, 1, curvs)
   
         return roots,vals[0],curvs
   
      def intercept(self, y):
         '''Find the value of x for which the interpolator goes through [y]'''
         if self.realization is not None:
            poly = self.realization - y
         else:
            poly = self.poly - y
         # Roots can be complex.  Want only real ones
         roots = poly.roots()
         gids = num.isreal(roots)
         roots = num.real(roots[gids])
         gids = num.greater_equal(roots, self.poly.domain[0])*\
               num.less_equal(roots, self.poly.domain[1])
         roots = roots[gids]
         if len(roots) == 0:
            return None
         else:
            return roots

   for t in polytypes:
      if t == 'polynomial':
         functions[t] = (Polynomial, "Nth order simple polynomial (numpy.Polynomial)")
      else:
         functions[t] = (Polynomial, "Nth order %s polynomial (numpy.Polynomial)" % t)


if spline2 is not None:
   class HyperSpline(oneDcurve):
   
      def __init__(self, x, y, dy, **args):
         '''Fit a spline2 (Thijsse) to the data.  [args] can be any argument
         recognized by spline2.spline2()'''
   
         oneDcurve.__init__(self, x, y, dy)
   
         self.pars = {
               'xrange':None,
               'degree':3,
               'acfsearch':0,
               'acffunc':'exp',
               'ksi':None,
               'n':None,
               'allownonopt':1,
               'lopt':None,
               'rejlev':0.05,
               'xlog':0}
         for key in args:
            if key not in self.pars:
               raise TypeError, \
                     "%s is an invalid keyword argument for this method" % key
            self.pars[key] = args[key]
   
         # Make sure the data conform to the spine requirements
         self._regularize()
         self._setup()
         self.realization = None
   
      def __str__(self):
         return "Hyperspline"

      def help(self):
         print "The following parameters affect the fit.  See help(spline2) for"
         print "more information"
         print "xrange:      tuple of (xmin,xmax) over which to fit"
         print "lopt:        Force knot optimization to start with lopt knots"
         print "degree:      order of the spline (default 3)"
         print "xlog:        0/1, apply log10() to the x values before fitting?"
         print "rejlev:      Use rejection level on statistical tests of rejlev"
         print "allownonopt: 0/1 Allow splines with non-optimized breakpoints?"
         print "acfsearch:   0/1, weather to search for auto-correlation"
         print "acffunc:     functional form of autocorrelation (default exp)"
         print "ksi:         specify auto-correlation length"
         print "n:           only search for autocorrelation on index scale n"

      def _regularize(self):
         # x-values need to be strictly ascending.
         sids = num.argsort(self.x)
         self.x = self.x[sids]
         self.y = self.y[sids]
         self.ey = self.ey[sids]
   
         # here's some Numeric magic.  first, find where we have repeating x-values
         Nmatrix = num.equal(self.x[:,num.newaxis], self.x[num.newaxis,:])
         val_matrix = self.y[:,num.newaxis]*Nmatrix
         e_matrix = self.ey[:,num.newaxis]*Nmatrix
    
         average = num.sum(val_matrix, axis=0)/sum(Nmatrix)
         e_average = num.sum(e_matrix, axis=0)/sum(Nmatrix)
         
         # at this point, average is the original data, but with repeating data points
         # replaced with their average.  Now, we just pick out the unique x's and
         # the first of any repeating points:
         gids = num.concatenate([[True], num.greater(self.x[1:] - self.x[:-1], 0.)])
         self.x = self.x[gids]
         self.y = average[gids]
         self.ey = e_average[gids]
   
      def _setup(self):
         '''Given the current set of params, setup the interpolator.'''
         self.tck = spline2(self.x, self.y, w=1.0/self.ey, **self.pars)
         self.setup = True
         self.realization = None
   
      def __call__(self, x):
         '''Interpolate at point [x].  Returns a 3-tuple: (y, mask) where [y]
         is the interpolated point, and [mask] is a boolean array with the same
         shape as [x] and is True where interpolated and False where extrapolated'''
         if not self.setup:
            self._setup()
   
         if len(num.shape(x)) < 1:
            scalar = True
         else:
            scalar = False
   
         x = num.atleast_1d(x)
         if self.realization:
            evm = num.atleast_1d(evalsp(x, self.realization))
            mask = num.greater_equal(x, self.realization[0][0])*\
                   num.less_equal(x,self.realization[0][-1])
         else:
            evm = num.atleast_1d(evalsp(x, self.tck))
            mask = num.greater_equal(x, self.tck[0][0])*num.less_equal(x,self.tck[0][-1])
   
         if scalar:
            return evm[0],mask[0]
         else:
            return evm,mask
         
      def domain(self):
         return (self.tck[0][0], self.tck[0][-1])

      def draw(self):
         '''Generate a random realization of the spline, based on the data.'''
         y_draw = num.random.normal(self.y, self.ey)
         self.realizations.append(\
               spline2(self.x, y_draw, w=1.0/self.ey, **self.pars))
         if len(self.realizations) > self.num_real_keep:
            self.realizations = self.realizations[1:]
         self.realization = self.realizations[-1]
    
      def reset_mean(self):
         self.realization = None
   
      def rchisquare(self):
         chisq = self.chisquare()
         return chisq/(len(self.x) - len(self.tck[0]) - 1)
   
      def deriv(self, x, n=1):
         '''Returns the nth derivative of the function at x.'''
         if self.realization:
            tck = self.realization
         else:
            tck = self.tck
         if len(num.shape(x)) < 1:
            scalar = True
         else:
            scalar = False
         x = num.atleast_1d(x)
         if self.realization:
            evm = num.atleast_1d(evalsp(x, self.realization, deriv=n))
         else:
            evm = num.atleast_1d(evalsp(x, self.tck, deriv=n))
   
         if scalar:
            return evm[0]
         else:
            return evm
   
      def find_extrema(self, xmin=None, xmax=None):
         '''Find the position and values of the maxima/minima.  Returns a tuple:
            (roots,vals,ypps) where roots are the x-values where the extrema
            occur, vals are the y-values at these points, and ypps are the
            2nd derivatives.  Optionally specify the range over which maxima
            are valid.'''
         if self.realization:
            vals = eval_extrema(self.realization)
         else:
            vals = eval_extrema(self.tck)
         gids = num.ones(vals.shape, dtype=num.bool)
         if xmin is not None:
            gids = gids*num.greater_equal(vals[0],xmin)
         if xmax is not None:
            gids = gids*num.less_equal(vals[0],xmax)
         return (vals[0][gids], vals[1][gids], vals[2][gids])
   
      def intercept(self, y):
         '''Find the value of x for which the interpolator goes through [y]'''
   
         if self.realization:
            return eval_x(y, self.realization)
         else:
            return eval_x(y, self.tck)
   functions['hyperspline'] = (HyperSpline, "B. Thijsse style spline (snpy.spline2)")


class Spline(oneDcurve):

   def __init__(self, x, y, dy, **args):
      '''Fit a scipy (Dierkx) spline to the data.  [args] can be any argument
      recognized by scipy.interpolate.splrep.'''

      oneDcurve.__init__(self, x, y, dy)

      self.pars = {
            't':None,
            'k':3,
            's':None,
            'xb':None,
            'xe':None,
            'task':0}
      for key in args:
         if key not in self.pars:
            raise TypeError, \
                  "%s is an invalid keyword argument for this method" % key
         self.pars[key] = args[key]

      # Make sure the data conform to the spine requirements
      self._regularize()
      self._setup()
      self.realization = None

   def __str__(self):
      return "Spline"

   def help(self):
      print "Parameters that affect the fit are:"
      print "k:      order of the spline (default 3)"
      print "task:   0,1,-1  (see scipy.interpolate.splrep)"
      print "s:      Smoothing length"
      print "t:      specify array of knots (task=-1)"
      print "xb      lower bound on x for fitting"
      print "xe      upper bound on x for fitting"

   def _regularize(self):
      # x-values need to be strictly ascending.
      sids = num.argsort(self.x)
      self.x = self.x[sids]
      self.y = self.y[sids]
      self.ey = self.ey[sids]

      # here's some Numeric magic.  first, find where we have repeating x-values
      Nmatrix = num.equal(self.x[:,num.newaxis], self.x[num.newaxis,:])
      val_matrix = self.y[:,num.newaxis]*Nmatrix
      e_matrix = self.ey[:,num.newaxis]*Nmatrix
 
      average = num.sum(val_matrix, axis=0)/sum(Nmatrix)
      e_average = num.sum(e_matrix, axis=0)/sum(Nmatrix)
      
      # at this point, average is the original data, but with repeating data points
      # replaced with their average.  Now, we just pick out the unique x's and
      # the first of any repeating points:
      gids = num.concatenate([[True], num.greater(self.x[1:] - self.x[:-1], 0.)])
      
      self.x = self.x[gids]
      self.y = average[gids]
      self.ey = e_average[gids]

   def _setup(self):
      '''Given the current set of params, setup the interpolator.'''
      self.tck = splrep(self.x, self.y, 1.0/self.ey, **self.pars)
      self.setup = True
      self.realization = None

   def __call__(self, x):
      '''Interpolate at point [x].  Returns a 3-tuple: (y, mask) where [y]
      is the interpolated point, and [mask] is a boolean array with the same
      shape as [x] and is True where interpolated and False where extrapolated'''
      if not self.setup:
         self._setup()

      if len(num.shape(x)) < 1:
         scalar = True
      else:
         scalar = False

      x = num.atleast_1d(x)
      if self.realization:
         evm = num.atleast_1d(splev(x, self.realization))
         mask = num.greater_equal(x, self.realization[0][0])*\
                num.less_equal(x,self.realization[0][-1])
      else:
         evm = num.atleast_1d(splev(x, self.tck))
         mask = num.greater_equal(x, self.tck[0][0])*num.less_equal(x,self.tck[0][-1])

      if scalar:
         return evm[0],mask[0]
      else:
         return evm,mask
      
   def domain(self):
      return (self.tck[0][0], self.tck[0][-1])

   def draw(self):
      '''Generate a random realization of the spline, based on the data.'''
      k = self.tck[2]
      y_draw = num.random.normal(self.y, self.ey)
      args = self.pars.copy()
      args['task'] = -1
      args['t'] = self.tck[0][k+1:-(k+1)]
      self.realizations.append(splrep(self.x, y_draw, self.ey, **args))
      if len(self.realizations) > self.num_real_keep:
         self.realizations = self.realizations[1:]
      self.realization = self.realizations[-1]
 
   def reset_mean(self):
      self.realization = None

   def rchisquare(self):
      chisq = self.chisquare()
      return chisq/(len(self.x) - len(self.tck[0]) - 1)

   def deriv(self, x, n=1):
      '''Returns the nth derivative of the function at x.'''
      if self.realization:
         tck = self.realization
      else:
         tck = self.tck
      if len(num.shape(x)) < 1:
         scalar = True
      else:
         scalar = False
      x = num.atleast_1d(x)
      if self.realization:
         evm = num.atleast_1d(splev(x, self.realization, der=n))
      else:
         evm = num.atleast_1d(splev(x, self.tck, der=n))

      if scalar:
         return evm[0]
      else:
         return evm

   def find_extrema(self, xmin=None, xmax=None):
      '''Find the position and values of the maxima/minima.  Returns a tuple:
         (roots,vals,ypps) where roots are the x-values where the extrema
         occur, vals are the y-values at these points, and ypps are the
         2nd derivatives.  Optionall, search only betwwen xmin and xmax.'''
      #evaluate the 1st derivative at k+1 intervals between the knots

      if self.realization:
         t,c,k = self.realization
      else:
         t,c,k = self.tck
      if xmax is None:  xmax = t[-1]
      if xmin is None:  xmin = t[0]
      x0s = t[k:-k]
      xs = []
      for i in range(len(x0s)-1):
         xs.append(num.arange(x0s[i],x0s[i+1],(x0s[i+1]-x0s[i])/(k+1)))
      xs = num.concatenate(xs)
      yps = self.deriv(xs, n=1)
      # now find the roots of the 1st derivative
      tck2 = splrep(xs, yps, k=3, s=0)
      roots = sproot(tck2)
      curvs = []
      vals = []
      for root in roots:
         vals.append(self.__call__(root)[0])
         curvs.append(self.deriv(root, n=2))
      gids = num.greater_equal(roots,xmin)*num.less_equal(roots,xmax)
      curvs = num.where(num.equal(curvs,0), 0, curvs/num.absolute(curvs))
      return roots[gids],num.array(vals)[gids],num.array(curvs)[gids]

   def intercept(self, y):
      '''Find the value of x for which the interpolator goes through [y]'''

      # use a fun little trick:
      if self.realization:
         tck = self.realization[0][::],self.realization[1]-y,\
               self.realization[2]
      else:
         tck = self.tck[0][::],self.tck[1]-y,self.tck[2]
      roots = sproot(tck)
      gids = num.greater_equal(roots,tck[0][0])*num.less_equal(roots,tck[0][-1])

      return sproot(tck)

functions['spline'] = (Spline, "Dierckx style splines (FITPACK)")

if pymc is not None:
   class GaussianProcess(oneDcurve):
   
      def __init__(self, x, y, dy, **args):
         '''Fit a GP (Gaussian Process) spline to the data.  [args] can be any argument
         recognized by scipy.interpolate.splrep.'''
   
         oneDcurve.__init__(self, x, y, dy)
   
         self.pars = {
               'diff_degree':None,
               'scale':None,
               'amp':None}
         for key in args:
            if key not in self.pars:
               raise TypeError, \
                     "%s is an invalid keyword argument for this method" % key
            self.pars[key] = args[key]
         #if 'func' in args:
         #   self.func = args['func']
         #else:
         #   self.func = lambda x:  x*0 + num.median(self.y)
   
         # Make sure the data conform to the spine requirements
         self.median = num.median(self.y)
         self._regularize()
         self._setup()
         self.realization = None
   
      def __str__(sef):
         return "Gaussian Process"

      def func(self, x):
         return x*0 + self.median

      def __getstate__(self):
         # we need to define this because Mean and Cov are not pickleable
         dict = self.__dict__.copy()
         if 'M' in dict:  del dict['M']
         if 'C' in dict:  del dict['C']
         # Setting setup to None will force re-generation of M and C
         #  when we are unpickled
         dict['setup'] = False
         return dict

      def help(self):
         print "Parameters that affect the fit:"
         print "scale:       Scale over which the function varies"
         print "amp:         Amplitude of typical function variations"
         print "diff_degree: Roughly, the degree of differentiability"

      def _regularize(self):
         # x-values need to be strictly ascending.  This isn't strictly required for
         #  GP's, but we want to elminiate duplicates, so might as well.
         sids = num.argsort(self.x)
         self.xdata = self.xdata[sids]
         self.ydata = self.ydata[sids]
         self.eydata = self.eydata[sids]
   
         # here's some Numeric magic.  first, find where we have repeating x-values
         Nmatrix = num.equal(self.xdata[:,num.newaxis], self.xdata[num.newaxis,:])
         val_matrix = self.ydata[:,num.newaxis]*Nmatrix
         e_matrix = self.eydata[:,num.newaxis]*Nmatrix
    
         average = num.sum(val_matrix, axis=0)/sum(Nmatrix)
         e_average = num.sum(e_matrix, axis=0)/sum(Nmatrix)
         
         # at this point, average is the original data, but with repeating data points
         # replaced with their average.  Now, we just pick out the unique x's and
         # the first of any repeating points:
         gids = num.concatenate([[True], num.greater(self.xdata[1:] - self.xdata[:-1], 0.)])
         
         self.xdata = self.xdata[gids]
         self.ydata = average[gids]
         self.eydata = e_average[gids]
   
         self.vardata = num.power(self.ey,2)
   
      def _setup(self):
         '''Given the current set of params, setup the interpolator.'''
   
         if self.diff_degree is None:
            self.diff_degree = 2
   
         if self.amp is None:
            self.amp = num.std(self.y)
   
         if self.scale is None:
            self.scale = (self.x.max() - self.x.min())/2
   
         self.M = GP.Mean(self.func)
         self.C = GP.Covariance(GP.matern.euclidean, diff_degree=self.diff_degree,
                           amp=self.amp, scale=self.scale)
   
         GP.observe(self.M, self.C, obs_mesh=self.x, obs_vals=self.y, obs_V=self.var)
         self.setup = True
         self.realization = None
   
      def __call__(self, x):
         '''Interpolate at point [x].  Returns a 3-tuple: (y, mask) where [y]
         is the interpolated point, and [mask] is a boolean array with the same
         shape as [x] and is True where interpolated and False where extrapolated'''
         if not self.setup:
            self._setup()
   
         if len(num.shape(x)) < 1:
            scalar = True
         else:
            scalar = False
   
         x = num.atleast_1d(x)
         if self.realization is not None:
            res = self.realization(x)
         else:
            res = self.M(x)
   
         if scalar:
            return res[0],self.x.min() <= x[0] <= self.x.max()
         else:
            return res,num.greater_equal(x, self.x.min())*\
                  num.less_equal(x, self.x.max())
         
      def domain(self):
         return (self.x.min(),self.x.max())

      def error(self, x):
         '''Returns the error in the interpolator at points [x].'''
         if not self.setup:
            self._setup()
   
         if len(num.shape(x)) < 1:
            scalar = True
         else:
            scalar = False
   
         x = num.atleast_1d(x)
         res = num.sqrt(self.C(x))
   
         if scalar:
            return res[0]
         else:
            return res
   
      def draw(self):
         '''Generate a random realization of the spline, based on the data.'''
         self.realization = GP.Realization(self.M, self.C)
    
      def reset_mean(self):
         self.realization = None
   
      def rchisquare(self):
         chisq = self.chisquare()
         return chisq/(len(self.x) - 4)
   
      def deriv(self, x, n=1):
         '''Returns the nth derivative of the function at x.'''
         if len(num.shape(x)) < 1:
            scalar = True
         else:
            scalar = False
         xs = num.atleast_1d(x)
         f = lambda x:  self.__call__(x)[0]
         res = deriv(f, xs, dx=self.scale/100., n=n)
   
         if scalar:
            return res[0]
         else:
            return res
   
      def find_extrema(self, xmin=None, xmax=None):
         '''Find the position and values of the maxima/minima.  Returns a tuple:
            (roots,vals,ypps) where roots are the x-values where the extrema
            occur, vals are the y-values at these points, and ypps are the
            2nd derivatives.  Optionally, only search for roots between
            xmin and xmax'''
         #evaluate the 1st derivative at sacle/10 intervals (that should be 
         #   enough)
         if xmin is None:  xmin = self.x.min()
         if xmax is None:  xmax = self.x.max()
         dx = min(self.scale/20, (xmax-xmin)/5.0)
         xs = num.arange(xmin, xmax+dx, dx)
         #f = lambda x: self.__call__(x)[0]
         dys = self.deriv(xs, n=1)
         #dys = num.diff(self.__call__(xs)[0])
         #adys = (dys[1:] + dys[:-1])/2
         #dys = num.concatenate([[dys[0]],adys,[dys[-1]]])
         pids = num.greater(dys, 0)
         inds = num.nonzero(pids[1:] - pids[:-1])[0]
   
         if len(inds) == 0:
            return (None,None,None)
         ret = []
         for i in range(len(inds)):
            #try:
            res = brentq(self.deriv, xs[inds[i]], xs[inds[i]+1])
            #res = newton(self.deriv, xs[inds[i]])
            ret.append(res)
            #except:
            #   continue
         ret = num.array(ret)
         vals = self.__call__(ret)[0]
         curvs = self.deriv(ret, n=2)
         curvs = num.where(curvs > 0, 1, curvs)
         curvs = num.where(curvs < 0, -1, curvs)
   
         return ret,vals,curvs
   
      def intercept(self, y):
         '''Find the value of x for which the interpolator goes through [y]'''
   
         xs = num.arange(self.x.min(), self.x.max(), self.scale/10)
         f = lambda x:  self.__call__(x)[0] - y
         ys = f(xs)
   
         pids = num.greater(ys, 0)
         if num.alltrue(pids) or num.alltrue(-pids):
            return None
   
         ret = []
         inds = num.nonzero(pids[1:] - pids[:-1])[0]
         for i in range(len(inds)):
            ret.append(brentq(f, xs[inds[i]], xs[inds[i]+1]))
         ret = num.array(ret)
         return ret
   functions['gp'] = (GaussianProcess, "Gaussian Process (pymc.GP)")
else:
   GaussianProcess = None


def Interpolator(type, x, y, dy, **args):
   '''Convenience function that returns a 1D interpolator of the given [type]
   if possible.'''
   if type not in functions.keys():
      raise ValueError,"Error:  the type %s is not defined.  Try list_types()"
   else:
      interp = functions[type][0]
      if type in polytypes:
         args['type'] = type
      return interp(x, y, dy, **args)

def list_types():
   '''Returns a list of 1D interpolators that are defined at load-time.'''
   l = functions.keys()
   l.sort()
   for t in l:
      print "%-15s%s" % (t+":",functions[t][1])

