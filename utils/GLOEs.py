'''GLOESpy

Algorithm:  Barry Madore
Python version:  Chris Burns

Requires:  Numeric, 
           fit_poly (http://www.ociw.edu/Code/python/burns-python-scripts/fit_poly.py)

A routine for generating a smooth 1D or 2D interpolation based on irregularly
sampled data.  The basic idea is simple:

- You've got a bunch of data (x, y +/- dy)
- consider a point where you want to interpolat (x0)
- generate a normalized Gaussian window function W(x0,sigma) centered at 
  x0 with a certain width sigma.
- weight the data like    W(x0,sigma)/dy
- Fit a polynomial at x0 and use that to interpolated the value y0.
- Repeat for each point x0.

For a 2D surface, you use an elliptical Gaussian window function and fit a
2D polynomial.'''
from numpy import *
from numpy import isinf,divide
import sys
from fit_poly import fit2Dpoly
from fit_poly import fitpoly
#from pygplot import *

debug=0

def smooth(x, y, weight, xeval, sigma=None, N=None):
   '''Given a set of x and y points, with associated internal errors in y given
   by erry, find a smothing esimate of the data using Barry's GLOEs 
   algorithm evaluated at the points in xeval.  The Gaussian window function 
   will have a sigma=sigma.  '''
   x,y,weight,xeval = map(asarray, [x,y,weight,xeval])

   if len(x) == 1:
      # only one point?  Just return it' value
      return xeval*0 + y[0], xeval*0+1.0/weight[0], xeval*0, xeval*0
   if len(x) == 2:
      # Two points?  do a linear interpolation
      sl = divide(y[1] - y[0], x[1] - x[0])
      if isinf(sl):
         raise ValueError, "slope is infinite"
      err2 = power(((xeval-x[0])/(x[1]-x[0])+1)/weight[0],2) + \
             power((xeval-x[0])/(x[1]-x[0])/weight[1],2)
      return y[0] + (xeval-x[0])*sl, sqrt(err2), xeval*0+sl, xeval*0

   # dists[i,j] is the distance from point xeval[i] to data point x[j]
   dists = x[NewAxis, :] - xeval[:, NewAxis]
   if sigma is not None:
      # Use the same weighting function everywhere
      sigmas = xeval*0+sigma
   else:
      FWHMs = []
      for i in range(len(dists)):
         lids = less(dists[i], 0)
         rids = 1-lids
         # Adapt the sigma to the sparcity of data
         sorts1 = sort(absolute(compress(lids, dists[i])))
         sorts2 = sort(absolute(compress(rids, dists[i])))
         if len(sorts1) >= N and len(sorts2) >= N:
            FWHMs.append(max(sorts1[N-1],sorts2[N-1]))
         elif len(sorts1) >= N:
            FWHMs.append(sorts1[N-1])
         elif len(sorts2) >= N:
            FWHMs.append(sorts2[N-1])
         else:
            raise RuntimeError, "Error:  N too large for this many points"
      sigmas = 10*array(FWHMs)

   # now weight these distances with a Gaussian
   arg = -0.5*power(dists, 2)/power(sigmas[:,NewAxis],2)
   arg = where(less(arg, -100), -100, arg)
   err_d = exp(arg)
   # normalize
   #norm = sum(err_d, axis=1)
   norm = maximum.reduce(err_d, axis=1)
   err_d = err_d/norm[:, NewAxis]
   weight = err_d*power(weight[NewAxis,:],2)

   # Lists of stuff we want to keep
   interps = []
   e_interps = []
   bs = []
   cs = []
   for i in range(len(xeval)):
      params,eparams = fitpoly(x, y, weight[i], x0=xeval[i], k=3)
      interps.append(params[0])
      e_interps.append(eparams[0])
      bs.append(params[1])
      cs.append(params[2])

   # Convert the lists arrays
   interps,e_interps,bs,cs = map(array, [interps,e_interps,bs,cs])
   return(interps, e_interps, bs, cs)


def smooth2d(x, y, z, wt, xeval, yeval, sigmax=None, sigmay=None,
      Nx=10, Ny=10, tol=0):
   '''Given a set of x, y, z points, with associated internal errors in z given
   by 1/wt, find a smothing esimate of the data using Barry's GLOEs 
   algorithm evaluated at the points in xeval, yeval.  The Gaussian window 
   function will have a sigma=sigmax in the x-direction and sigmay in the
   y-direction.'''

   x,y,z,wt,xeval,yeval = map(asarray, [x,y,z,wt,xeval,yeval])
   # We re-scale the problem to   0 < x 1 and 0 < y < 1
   x0 = min(x);  x1 = max(x)
   y0 = min(y);  y1 = max(y)

   # Rescaled data
   u = (x - x0)/(x1 - x0)
   v = (y - y0)/(y1 - y0)

   # Rescaled evaluation points
   si = (xeval - x0)/(x1 - x0)
   sj = (yeval - y0)/(y1 - y0)

   # rescaled sigmas:
   if sigmax is not None:  
      if len(shape(sigmax)) == 0:
         sigmax = sigmax + 0.0*si
      sigmax = sigmax/(x1-x0)
   if sigmay is not None:  
      if len(shape(sigmay)) == 0:
         sigmay = sigmay + 0.0*sj
      sigmay = sigmay/(y1-y0)


   interps = []
   e_interps = []
   fxs = []
   fys = []
   fxxs = []
   fyys = []
   fxys = []
   sigmaxs = []
   sigmays = []
   step = len(si)/79 + 1
   for i in range(len(si)):
      # [xy]dists[j] is the distance from point (si[i],sj[i]) to point u[j],v[j]
      xdists = u - si[i]
      ydists = v - sj[i]
      sids = None
 
      if sigmax is None:
         # Find the Nx-th unique x-distance sorted by euclidean distance
         dists2 = power(xdists,2) + power(ydists,2)
         sids = argsort(dists2)
         tempx = take(xdists,sids)
         gotem = []
         for dx in tempx:
            if dx not in gotem:
               gotem.append(dx)
               if len(gotem) == Nx:  break
         FWHMx = maximum.reduce(absolute(gotem))
         sigx = FWHMx/1.665
      else:
         sigx = sigmax[i]

      if sigmay is None:
         # Find the Nx-th unique x-distance sorted by euclidean distance
         if sids is None:
            dists2 = power(xdists,2) + power(ydists,2)
            sids = argsort(dists2)
         tempy = take(ydists,sids)
         gotem = []
         for dy in tempy:
            if dy not in gotem:
               gotem.append(dy)
               if len(gotem) == Ny:  break
         FWHMy = maximum.reduce(absolute(gotem))
         sigy = FWHMy/1.665
      else:
         sigy = sigmay[i]
 
      sigmaxs.append(sigx)
      sigmays.append(sigy)
      dists2 = -0.5*(power(xdists,2)/sigx**2 + power(ydists,2)/sigy**2)
      gids = greater(dists2, -100)
      err_d = exp(compress(gids, dists2))
      err_d = err_d/sum(err_d)
      weight = err_d*(compress(gids, wt))
      #print '-'*70
      #print i, xeval[i],yeval[i]
      #for j in range(len(weight)):
      #   print "   ",gids[j], compress(gids, x)[j], compress(gids, y)[j], compress(gids, xdists)[j], compress(gids, ydists)[j], sigx, sigy
      params,e_params = fit2Dpoly(compress(gids,u), compress(gids,v), compress(gids, z),
            weight, k=2, x0=si[i], y0=sj[i])
      interps.append(params[0])
      e_interps.append(e_params[0])
      fxs.append(params[1])
      fxxs.append(params[2])
      fys.append(params[3])
      fxys.append(params[4])
      fyys.append(params[5])

   interps,e_interps,fxs,fys,fxxs,fyys,fxys,sigmaxs,sigmays = map(array,
         [interps,e_interps,fxs,fys,fxxs,fyys,fxys,sigmaxs,sigmays])

   return(interps,e_interps,fxs,fys,fxxs,fyys,fxys,sigmaxs,sigmays)

class line:
   '''A class to represent a line of data.  The construtor takes x,y,dy data
   as arguments (and optionally the smoothing length).  Member functions 
   allow the user to interpolate on the line.'''
   
   def __init__(self, x,y,dy,sigma=None, N=None):
      self.xdata,self.ydata,self.dydata = map(asarray, [x,y,dy])

      if not (len(self.xdata) == len(self.ydata) == len(self.dydata)):
         raise AttributeError, "x,y,dy must all have same length"

      self.sigma = sigma
      self.N = N

   def chisq(self, sigma=None, N=None):
      '''Compute chisq for a supplied sigma or one stored in the instance.'''
      if sigma is None:
         if self.sigma is None:
            if N is None:
               if self.N is None:
                  raise AttributeError, "Must specify a smoothing length"
               else:
                  N = self.N
         else:
            sigma = self.sigma
      res = self.eval(self.xdata, sigma=sigma, N=N)
      rchisq = sum(power((res[0] - self.ydata)/self.dydata, 2))/len(self.xdata)
      return(rchisq)

   def find_sigma(self, low, high, ds):
      '''Given a range of sigmas (low,high) and interval (ds), compute rchisq
      for each and return the smoothing length that gives rchisq = 1 (or as
      close as possible).'''
      sigs = arange(low, high, ds)
      chisqs = array([self.chisq(sigma=sig) for sig in sigs])
      id = argsort(absolute(chisqs - 1))[0]
      return(sigs[id])

   def find_N(self, low, high):
      '''Given a range of N (low, high), compute rchisq for each and return
      the smoothing length that gives rchisq closest to 1.'''
      max_N = len(self.xdata)/2+1
      if low > max_N:
         raise ValueError, "Error, low must be less than %d" % (max_N)
      if high > max_N:
         high = max_N
      Ns = range(int(low), int(high)+1)
      chisqs = array([self.chisq(N=N) for N in Ns])
      id = argsort(absolute(chisqs - 1))[0]
      return(Ns[id])
      

   def eval(self, x, sigma=None, N=None):

      if not len(shape(x)):
         scalar = 1
         x = [x]
      else:
         scalar = 0
      x = asarray(x)
      if sigma is None:
         if self.sigma is None:
            if N is None:
               if self.N is None:
                  raise AttributeError, "Must specify a smoothing length"
               else:
                  N = self.N
         else:
            sigma = self.sigma
      
      res = smooth(self.xdata, self.ydata, 1.0/self.dydata, x, sigma=sigma, N=N)
      self.x = x
      self.y = res[0]
      self.dy = res[1]
      self.fx = res[2]
      self.fxx = res[3]

      if scalar:
         return(self.y[0], self.dy[0])
      else:
         return(self.y, self.dy)


class surface:
   '''A class to represent a surface of data.  The construtor takes x,y,z,dz
   data as arguments (and optionally the smoothing lengths).  Member functions 
   allow the user to interpolate on the surface.'''
   
   def __init__(self, x,y,z,dz,sigmax=None,sigmay=None, Nx=None, Ny=None):

      self.xdata,self.ydata,self.zdata,self.dzdata = map(asarray, [x,y,z,dz])

      if not (len(self.xdata) == len(self.ydata) == len(self.zdata) == \
            len(self.dzdata)):
         raise AttributeError, "x,y,dy must all have same length"

      self.sigmax = sigmax
      self.sigmay = sigmay

      self.Nx = Nx
      self.Ny = Ny

   def eval(self, x, y, sigmax=None, sigmay=None, Nx=None, Ny=None, tol=0):

      if not len(shape(x)):
         scalar = 1
         x = [x]
      else:
         scalar = 0

      if not len(shape(y)):
         scalar = 1
         y = [y]
      else:
         scalar = 0

      x = asarray(x)
      y = asarray(y)
      if not (len(x) == len(y)):
         raise RuntimeError, "x and y must have same length"

      if sigmax is None:
         if self.sigmax is None:
            if Nx is None:
               if self.Nx is None:
                  raise AttributeError, "Must specify a x smoothing length"
               else:
                  Nx = self.Nx
         else:
            sigmax = self.sigmax

      if sigmay is None:
         if self.sigmay is None:
            if Ny is None:
               if self.Ny is None:
                  raise AttributeError, "Must specify a y smoothing length"
               else:
                  Ny = self.Ny
         else:
            sigmay = self.sigmay

      res = smooth2d(self.xdata, self.ydata, self.zdata, 1.0/self.dzdata, x, y,
            sigmax=sigmax, sigmay=sigmay, Nx=Nx, Ny=Ny, tol=tol)
      self.x = x
      self.y = y
      self.z = res[0]
      self.dz = res[1]
      self.fx = res[2]
      self.fy = res[3]
      self.fxx = res[4]
      self.fyy = res[5]
      self.fxy = res[6]
      self.sigmaxs = res[7]
      self.sigmays = res[8]

      return(self.z, self.dz)
