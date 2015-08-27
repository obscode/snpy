''' Spline2.py: wrapper for B. Thijsse et al.'s hyper-spline routines.

Yet another spline interpolation routine.  The problem:  given a set of
experimental data with noise, find the spline with the optimal number of
knots.

Solution: 
    They use the usual kind of routines to determine least-squares
    splines from a given set of knot points.  The problem REALLY
    boils down to:  how many knots do you use?  There are two 
    extremes:  put a knot point on each data point to get an
    interpolating spline (which sucks for experimental data with
    noise).  The other extreme is to have the minimal set of knots
    to define a polynomial of order k (e.g., a cubic).  This also
    sucks.  Somewhere between the two extremes is a number of
    knots that optimally recovers the information in the data and
    smooths out the noise.

    spline2 starts with a large number of knots (interpolating
    spline) and iteratively removes knots until a figure of merit
    reaches some prescribed value.  In this case, this figure of
    merit is the Durbin-Watson statistic, which measures the auto-
    correlation between the residuals of the spline fit.

For more details, see:
 *  Barend J. Thijsse et al., "A Practical Algorithm for Least-Squares 
    spline Approximation of Data Containing Noise", Computers in Physics, 
    vol 12 no. 4 July 1998
 *  http://structureandchange.3me.tudelft.nl/
'''
import spline2c
import numpy as num

acffuncs = ['exp','gauss','linear','sinc']

def spline2(x, y, w=None, sigma=None, rsigma=None, xrange=None, degree=3, 
      acfsearch=0, acffunc='exp', ksi=None, n=None, 
      allownonopt=1, lopt=None, rejlev=0.05, xlog=0, 
      interactive=0, verbose=0, full_output=0):
   '''Solve for the optimal spline given a set of (x,y) data.  

      Args:
         w (float array): specify weights (1/sigma) for each data point
         sigma (float): specify an absolute sigma for all data points:  
                   w[i] = 1/sigma
                   this will then force a regular chi-square fit instead of DW
         rsigma (float): specify a relative sigma for all data ponts:  
                   w[i] = 1/(y*rsigma)
         xrange (2-tuple): Only use data in interval (xrange[0], xrange[1])
         degree (int): degree of the spline (default:  3 for cubic spline)
         acfsearch (bool): perform an automated search for autocorrelation.
         acffunc (str):  Use acffunc as the autocorrelation function.  
                   can be one of 'exp','gauss','linear', or 'sinc'.  Default:
                   'exp'
         ksi (float): Use a specific autocorrelation length equal to ksi
         n (int):  Only search for autocorrelation on index interval n
         allownonopt (bool):  Allow splines with non-optimized breakpoints 
                   (default True)
         lopt (int): Force knot optimization to start with lpot knots.
         rejlev (float): Use rejection level on statistical tests of rejlev 
                   (default 0.05)
         xlog (bool): Take the log10(x) before spline fitting. Default: False
         verbose (bool):  Lots of output to stdout. Default: False
         interactive (bool):  Allows the user to choose the optimal spline 
                      manually.
         full_output (bool):  along with tck, return the following statistics:
                       rms, dws (Durbin-Watson statistic), lfin (final number
                       of knot points), ksi (computed auto-correlation length),
                       acffit (indicator of how well the assumed auto-
                       correlation function represents the data),
 
      Returns:
        (tuple):  (t, c, k)   if full_output=0 ((t,c,k), rms, dws, lfin, ksi, 
        acffit)  if full_output=1

        - t: array of lfin+1 knot points
        - c: array of lfin+k-1 spline coefficients
        - k: order of the spline (note:  order = degree+1, so this is 4 
             for a cubic spline!)

        The tuple (t,c,k) can be input to routines such as evalsp().

        - rms: dispersion of the fitted spline
        - dws: Durbin-Watson statistic
        - lfin: final number of knot points
        - ksi: computed auto-correlation length
        - acffit: how well the correlations agree with assumed funcion.
   '''

   x = num.asarray(x).astype(num.float64)
   y = num.asarray(y).astype(num.float64)
   xin = x.astype(num.float64)
   yin = y.astype(num.float64)
   if w is not None:
      w = num.asarray(w).astype(num.float64)
      win = w.astype(num.float64)
   else:
      win = x*0.0 + 1.0    # assume no weight info.
   if not (len(xin) == len(yin) == len(win)):
      raise IndexError, "Arrays x,y, and w must have same length"

   rel=0
   fixed_sigma=0
   fixval=0
   if sigma is not None or rsigma is not None:
      fixed_sigma = 1
      if sigma is not None:
         fixval = sigma
         rel = 0
      else:
         fixval = rsigma
         rel = 1


   acf_ind = acffuncs.index(acffunc) + 1
   if xrange is not None:
      xbegin = xrange[0]
      xend = xrange[1]
      xflag = 1
   else:
      xbegin = 0
      xend = 0
      xflag = 0

   if n is not None:
      nset = 1
      n_max = n_min = n
   else:
      nset = 0
      n_max = n_min = 1

   if ksi is not None:
      ksiset = 1
      ksibegin = ksiend = ksi
   else:
      ksiset = 0
      ksibegin = ksiend = 0.0

   if lopt is not None:
      lind = 1
   else:
      lind = 0
      lopt = 0
   


   result = spline2c.spline2(xin, yin, win, degree, acfsearch, acf_ind, xflag,
         xbegin,xend, 0, 0, 0, xlog, nset, ksiset, n_min, n_max, ksibegin,
         ksiend, rejlev, rel, fixed_sigma, fixval, lind, lopt, allownonopt,
         interactive, verbose)

   if full_output:
      return((result[0:3]), result[3:])
   else:
      return(result[0:3])

def evalsp(x, tck, deriv=0):
   '''Evaluate a spline computed by spline2.  
   
   Args:
      x (float array or scalar): The spline is evaluated at the points x.  
      tck (3-tuple):  a tuple of (knots, coefficents, k) that are returned 
                      as the first output of spline2.  
      deriv (int): if  > 0, compute the deriv-th derivative of the spline 
      
   Returns: 
      float array: evaluated interpolant or derivative thereof.
   '''
   if len(num.shape(x)) == 0:
      x = num.array([x]).astype(num.float64)
   else:
      x = num.asarray(x).astype(num.float64)

   xin = x.astype(num.float64)
   t = tck[0]
   c = tck[1]
   l = len(t) - 1
   k = tck[2]

   result = spline2c.evalsp(xin, k, l, t, c, deriv)
   return(result)

def eval_extrema(tck):
   '''Attempts to find the extrema of the spline (where S'(x)==0). 
   
   Args:
      tck (3-tuple): tuple of (knots, coefficients, k) that are returned as 
                     the first output of spline2.  
                     
   Returns:
      3-tuple:  (xextr, yextr, signs)  where 
      
      - xextr are the x-positions of the extrema, 
      - yextr are S(extr), and 
      - signs provice the sign of the 2nd derivative S''(extr):  
      
         - signs[i] < 0 --> maximum, 
         - signs[i] > 0 --> minimum, 
         - signs[i] close to 0 --> inflection point
   '''
   t = tck[0]
   c = tck[1]
   l = len(t) - 1
   k = tck[2]
   
   if k <= 2:
      raise ValueError, "Spline order must be at least 2 for finding extrema"
   
   result = spline2c.eval_extrema(k, l, t, c)
   return(result)

def eval_inflect(tck):
   '''Attempts to find the inflection points of the spline (where S''(x)==0).

   Args:
      tck (3-tuple):  tuple of (knots, coefficients, k) that are returned as 
                       the first output of spline2.  
                       
   Returns:
      3-tuple:  (xinflect, yinflect, dyinflect) where 
      
      - xinflect are the x-positions of the inflection points, 
      - yminflect are S(xinflect), and 
      - dyinflect are S'(xinflect).
   '''

   t = tck[0]
   c = tck[1]
   l = len(t) - 1
   k = tck[2]
   if k <= 3:
      raise ValueError, "Spline order must be at least 3 for finding inflections"

   result = spline2c.eval_inflect(k,l,t,c)
   return(result)

def eval_integ(x0, x1, tck):
   '''Evaluates the integral from x0 to x1 of the spline defined by tck.
   
   Args:
      x0 (float):  lower limit of integration
      x1 (float):  upper limit of integration
      tck (3-tuple): tuple of (knots, coefficients, k) that are returned as 
                     the first output of spline2.  

   Returns:
      float:  the integration.
   '''
   t = tck[0]
   c = tck[1]
   l = len(t) - 1
   k = tck[2]
   if x0 < t[0] or x1 > t[1]:
      raise ValueError, "integration limits beyond spline definition"

   result = spline2c.eval_integ(x0, x1, k, l, t, c)
   return(result)

def eval_x(value, tck):
   '''Attempts to find the roots of the equation S(x) = value for the spline
   defined by tck.  
   
   Args:
      value (float):  value to solve the root for
      tck (3-tuple): tuple of (knots, coefficients, k) that are returned as 
                     the first output of spline2.  

   Returns:
      float array: roots of the equation.
   '''
   t = tck[0]
   c = tck[1]
   l = len(t) - 1
   k = tck[2]
   result = spline2c.eval_x(value, k, l, t, c)
   return(result)
