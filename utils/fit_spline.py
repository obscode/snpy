'''This module provides some helper functions for fitting a spline to Supernova
data.'''

from scipy.interpolate import splrep
from scipy.interpolate import splev
from scipy.integrate import quad
import numpy as num
import sys,os,string
from snpy import spline2
#import fit_lightcurves as lc
import scipy


def make_spline(t, m, e_m, knots=None, k=1, s=None, fitflux=0, zpt=0, 
      tmin=None, tmax=None, task=-1, anchor_dist=[5.0,5.0], slopes=[None,None]):
   '''A wrapper around splrep that makes sure the independent variable is
   monotonic and non-repeating.  Required arguments:  time (t), magnitudes
   (m) and errors (e_m).  If knots are specified, use them (if task==-1), otherwise,
   they are computed from -10 days to 100 days in 10-day increments.  k is the 
   spline order (default 1) and s is the smoothing factor, as per splrep.  If fitflux
   is nonzero, convert magnitudes to flux (using provided zpt).  tmin and tmax should
   be set to the limits of the spline.'''
   # first, we make sure that t is monotonically increasing with no repeated
   # elements
   sids = num.argsort(t)
   tt = t[sids]      #num.take(t, sids)
   mm = m[sids]      # num.take(m, sids)
   ee_m = e_m[sids]  #num.take(e_m, sids)

   if tmin is None:
      tmin = t.min()
   if tmax is None:
      tmax = t.max()

   # here's some Numeric magic.  first, find where we have repeating x-values
   Nmatrix = num.equal(tt[:,num.newaxis], tt[num.newaxis,:])
   #val_matrix = mm[:,num.newaxis]*num.ones((len(mm),len(mm)))
   #e_matrix = ee_m[:,num.newaxis]*num.ones((len(mm),len(mm)))
   val_matrix = mm[:,num.newaxis]*Nmatrix
   e_matrix = ee_m[:,num.newaxis]*Nmatrix

   average = sum(val_matrix)/sum(Nmatrix)
   e_average = sum(e_matrix)/sum(Nmatrix)
   
   # at this point, average is the original data, but with repeating data points
   # replaced with their average.  Now, we just pick out the unique x's and
   # the first of any repeating points:
   gids = num.concatenate([[1], num.greater(tt[1:] - tt[:-1], 0.)])
   tt = num.compress(gids, tt)
   mm = num.compress(gids, average)
   ee_m = num.compress(gids, e_average)

   # Now get rid of any data that's outside [tmin,tmax]
   gids = num.less_equal(tt, tmax)*num.greater_equal(tt, tmin)
   #tt = num.compress(gids,tt)
   tt = tt[gids]
   #mm = num.compress(gids,mm)
   mm = mm[gids]
   #ee_m = num.compress(gids,ee_m)
   ee_m = ee_m[gids]
   ee_m = num.where(num.less(ee_m, 0.001), 0.001, ee_m)

   # Next, add some anchors to the data to control the slopes
   if anchor_dist[0] > 0 and task != -1:
      if slopes[0] is not None:
         mm0 = mm[0] - slopes[0]*anchor_dist[0]
      else:
         mm0 = mm[0] - (mm[1] - mm[0])/(tt[1]-tt[0])*anchor_dist[0]
      tt = num.concatenate([[tt[0]-anchor_dist[0]], tt])
      mm = num.concatenate([[mm0], mm])
      ee_m = num.concatenate([[ee_m[0]], ee_m])
   if anchor_dist[1] > 0:
      if slopes[1] is not None:
         mm1 = mm[-1] + slopes[1]*anchor_dist[1]
      else:
         mm1 = mm[-1] + (mm[-1] - mm[-2])/(tt[-1]-tt[-2])*anchor_dist[1]
      tt = num.concatenate([tt,[tt[-1]+anchor_dist[1]]])
      mm = num.concatenate([mm,[mm1]])
      ee_m = num.concatenate([ee_m,[ee_m[-1]]])

   # Now convert to flux if requested:
   if fitflux:
      mm = num.power(10, -0.4*(mm - zpt))
      ee_m = mm*ee_m/1.087
   
   if knots is None and task==-1:
      # Use the minimal number
      knots = tmin + num.arange(2*k+3)*(tmax - tmin)/(2*k+2)
   
   # Okay, now make the spline representation
   tck,fp,ier,msg = splrep(tt, mm, 1.0/ee_m, k=k, s=s, t=knots, task=task,
         full_output=1)
   return(tck,fp,ier,msg)

def K2(x, tck):
   '''compute the square curvature of a spline2 at point x'''
   yp = num.power(spline2.evalsp(x, tck, 1),2)
   ypp = num.power(spline2.evalsp(x, tck, 2),2)
   return(ypp/num.power(1+yp,3))


def make_spline2(t, m, e_m, k=3, fitflux=0, zpt=0, tmin=-10, tmax=100,
      adaptive=0, max_curv_fac=10, **args):
   '''A wrapper around spline2 that makes sure the independent variable is
   monotonic and non-repeating.  Required arguments:  time (t), magnitudes
   (m) and errors (e_m).  k is the spline order (default 3)  If fitflux
   is nonzero, convert magnitudes to flux (using provided zpt).  tmin and tmax should
   be set to the limits of the spline.'''
   # first, we make sure that t is monotonically increasing with no repeated
   # elements
   sids = num.argsort(t)
   tt = num.take(t, sids)
   mm = num.take(m, sids)
   ee_m = num.take(e_m, sids)

   # here's some Numeric magic.  first, find where we have repeating x-values
   Nmatrix = num.equal(tt[:,num.newaxis], tt[num.newaxis,:])
   #val_matrix = mm[:,num.newaxis]*num.ones((len(mm),len(mm)))
   #e_matrix = ee_m[:,num.newaxis]*num.ones((len(mm),len(mm)))
   val_matrix = mm[:,num.newaxis]*Nmatrix
   e_matrix = ee_m[:,num.newaxis]*Nmatrix

   average = sum(val_matrix)/sum(Nmatrix)
   e_average = sum(e_matrix)/sum(Nmatrix)
   
   # at this point, average is the original data, but with repeating data points
   # replaced with their average.  Now, we just pick out the unique x's and
   # the first of any repeating points:
   gids = num.concatenate([[True], num.greater(tt[1:] - tt[:-1], 0.)])
   tt = tt[gids]           #num.compress(gids, tt)
   mm = average[gids]      # num.compress(gids, average)
   ee_m = e_average[gids]  #num.compress(gids, e_average)

   # Now get rid of any data that's outside [tmin,tmax]
   gids = num.less_equal(tt, tmax)*num.greater_equal(tt, tmin)
   tt = tt[gids]             #num.compress(gids,tt)
   mm = mm[gids]             # num.compress(gids,mm)
   ee_m = ee_m[gids]         #num.compress(gids,ee_m)
   ee_m = num.where(num.less(ee_m, 0.001), 0.001, ee_m)

   # Now convert to flux if requested:
   if fitflux:
      mm = num.power(10, -0.4*(mm - zpt))
      ee_m = mm*ee_m/1.087
   
   # Okay, now make the spline representation
   if not adaptive:
      tck = spline2.spline2(tt, mm, w=1.0/ee_m, degree=k, **args)
      fp = num.sum(num.power((mm - spline2.evalsp(tt, tck))/ee_m, 2))
      ier = 0
      msg = 'not much'
      return(tck,fp,ier,msg)
   
   # Do an adaptive (much slower) search for the best fit, subject
   #  to curvature constraints (right now, just a cap).
   if 'lopt' in args:  del args['lopt']
   Ks = []
   chisqs = []
   lopts = range(2, int(0.8*len(tt) - 1))
   for l in lopts:
      tck = spline2.spline2(tt, mm, w=1.0/ee_m, degree=k, lopt=l, **args)
      fp = num.sum(num.power((mm - spline2.evalsp(tt, tck))/ee_m, 2))
      K = quad(K2, tck[0][0], tck[0][-1], args=(tck,), epsrel=0.01)[0]
      chisqs.append(fp)
      Ks.append(K)

   chisqs = num.array(chisqs)
   Ks = num.array(Ks)
   chisqs = num.where(num.less(Ks, Ks.min()*max_curv_fac), chisqs, num.Inf)
   id = num.argmin(chisqs)

   ier = lopts[id]
   tck = spline2.spline2(tt, mm, w=1.0/ee_m, degree=k, lopt=ier, **args)
   fp = num.sum(num.power((mm - spline2.evalsp(tt, tck))/ee_m, 2))
   return(tck, fp, ier, "Optimized lopt = %d" % ier)

def fit_spline(t, m, e_m, knots=None, k=1, s=None, fitflux=0, zpt=0, tmin=-10, tmax=100,
      task=-1):
   '''Fit a spline to a set of lightcurve data (t,m,e_m).  The rest of the arguments are
   the same as make_spline() above.  The spline is then evaluated at one day intervals
   from t[0] to t[-1] and returned as a two-tuple (time,mag).'''

   tck,fp,ier,msg = make_spline(t, m, e_m, knots, k, s, fitflux, zpt, tmin, tmax,task)

   # Now do the interpolation
   evt = num.arange(t[0],t[-1],1, typecode=num.float64)
   evm = splev(evt, tck)

   if fitflux:
      evm = -2.5*num.log10(evm) + zpt

   return(evt,evm)

def interp_spline(t, m, e_m, eval_t, knots=None, k=1, s=None, fitflux=0, zpt=0, 
      tmin=-10, tmax=100, task=-1, tol=0.1):
   '''Same as fit_spline, except that you decide on the interpolated times.  Good
   for making color estimates when data in different bands in at different epochs.
   In this case, only the interpolated magnitudes are returned.'''
   tck,fp,ier,msg = make_spline(t,m,e_m,knots,k,s,fitflux,zpt,tmin,tmax,task)
   evm = splev(eval_t, tck)

   # Now, we scan t and eval_t and find where they are less than tol.  In these
   #  cases, we take the average of any matching times
   delta = num.absolute(t[num.newaxis,:] - eval_t[:,num.newaxis])
   cond = num.less(delta, tol)
   values = num.array([m]*len(eval_t))*cond
   N = num.sum(cond, axis=1)
   w = num.where(N > 0, N, 1)
   s = num.sum(values, axis=1)
   evm = num.where(N > 0, s/w, evm)

   if fitflux:
      evm = -2.5*num.log10(evm) + zpt

   return(evm)

def find_extr(tck, numpoints=1000):

   xs = num.arange(numpoints+1, typecode='d')/numpoints*(tck[0][-1]- tck[0][0]) + tck[0][0]
   f = lambda x:  splev(x, tck, 1)
   derivs = f(xs)
   gt0 = num.greater(derivs, 0)
   inds = num.nonzero(gt0[1:] - gt0[:-1])[0]
   if len(inds) == 0:
      return (None,None,None)
   ret = []
   for i in range(len(inds)):
      ret.append(scipy.optimize.brentq(f, xs[inds[i]], xs[inds[i]+1]))
   ret = num.array(ret)
   vals = splev(ret, tck)
   curvs = splev(ret, tck, 2)
   curvs = num.where(num.greater(curvs, 0), 1, curvs)
   curvs = num.where(num.less(curvs, 0), -1, curvs)
   if len(ret) == 1:
      curvs = num.array([curvs])
      vals = num.array([vals])
   
   return (ret,vals,curvs)
