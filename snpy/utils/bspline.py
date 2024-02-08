'''A basis-spline module shamelessly grabbed from Stack Overvlow:

   http://stackoverflow.com/questions/35236221/how-to-get-the-spline-basis-used-by-scipy-interpolate-splev

This is needed because the basis functions are hidden from python view in
the scipy.interpolate implementation. We want the basis functions themselves,
so we have to make our own.'''

from __future__ import print_function
import numpy as np

def memo(f):
   # Peter Norvig's
   """Memoize the return value for each call to f(args).
   Then when called again with same args, we can just look it up."""
   cache = {}

   def _f(*args):
      try:
         return cache[args]
      except KeyError:
         cache[args] = result = f(*args)
         return result
      except TypeError:
         # some element of args can't be a dict key
         return f(*args)
   _f.cache = cache
   return _f


def bspline_basis(knots, u, degree, zeroslope=False, gradient=False,
      extrap=False):
   '''Compute the basis functions for a b-spline interpolation.

   Args:
      knots (list or array length N): the knot positions.
      u (list or array length M): the evaluation points of the spline
      degree (int): The degree of the spline (e.g., 3 == cubic)
      zeroslope (bool): the slopes at the end points are constrained to be zero
      gradient (bool): If True, the curvature of the spline is zero at
                       the endpoints (so could be extrapolated linearly)
      extrap (bool):  If True and zeroslope or gradient, extrapolate beyond
                      the end of the knot points

   Returns:
      2-d array with shape (M,X) where X is:
         N+degree-1    if slopes=None and gradient=False
         N+degree-3    otherwise
   '''
   # Create knot vector and a range of samples on the curve
   try:
      knots = np.asarray(knots)
   except:
      raise ValueError("knots must be a list type or array")
   if len(np.shape(knots)) != 1:
      raise ValueError("knots must be a 1d array/list")
   if knots.shape[0] < degree:
      raise ValueError("You must have at least <degree> internal knots")
   if not np.all(np.greater_equal(knots[1:] - knots[:-1],0)):
      raise ValueError("knots must be strictly increasing")
   kv = np.concatenate([[knots[0]]*degree,
                         knots,
                         [knots[-1]]*degree])
   off = kv[0]
   scale = (kv[-1]-kv[0])
   kv = (kv - off)/scale    # map to 0 -> 1
   #u = np.linspace(0, c - degree, n)  # samples range
   u = (u - off)/scale

   # Cox - DeBoor recursive function to calculate basis
   @memo
   def coxDeBoor(k, d):
      # Test for end conditions
      if (d == 0):
         return np.where(np.greater_equal(u-kv[k],0)*\
                         np.less(u-kv[k+1],0),1.0, 0.0)
         #return ((u - kv[k] >= 0) & (u - kv[k + 1] < 0)).astype(int)
      denom1 = kv[k + d] - kv[k]
      term1 = 0
      if denom1 > 0:
         term1 = ((u - kv[k]) / denom1) * coxDeBoor(k, d - 1)

      denom2 = kv[k + d + 1] - kv[k + 1]
      term2 = 0
      if denom2 > 0:
         term2 = ((-(u - kv[k + d + 1]) / denom2) * coxDeBoor(k + 1, d - 1))

      return term1 + term2

   # Compute basis for each point
   b = np.column_stack([coxDeBoor(k, degree) for k in range(len(kv)-degree-1)])

   # Now deal with the endpoint
   zid = np.nonzero(np.greater(b[:,-1],0))[0][-1]+1
   b[zid,-1] = 1.0

   if zeroslope:
      '''contrain slopes at endpoints to be zero. This combines bases, so 
      two less degrees of freedom'''
      bb = np.zeros((b.shape[0],b.shape[1]-2))
      bb[:,0] = b[:,0] + b[:,1]
      bb[:,-1] = b[:,-1] + b[:,-2]
      bb[:,1:-1] = b[:,2:-2]
      if extrap:
         # simply add constant term to end bases
         bb[:,0] = bb[:,0] + np.less(u, kv[0])*1.0
         bb[:,-1] = bb[:,-1] + np.greater(u, kv[-1])*1.0
      return bb
   elif gradient:
      '''constrain zero curvature at end points. Again, combines bases, so
      two less degrees of fredom'''
      dl1 = knots[1]-knots[0]
      dl2 = knots[2]-knots[0]
      dlm1 = knots[-1]-knots[-2]
      dlm2 = knots[-1]-knots[-3]
      bb = np.zeros((b.shape[0],b.shape[1]-2))
      bb[:,0] = (1+dl1/dl2)*b[:,0]+b[:,1]
      bb[:,1] = b[:,2] - dl1/dl2*b[:,0]
      if degree == 3 and len(knots) == 3:
         # special case of mixed modes
         bb[:,1] = bb[:,1] - dlm1/dlm2*b[:,-1]
      else:
         bb[:,2:-2] = b[:,3:-3]
         bb[:,-2] = b[:,-3] - b[:,-1]*dlm1/dlm2
      bb[:,-1] = b[:,-2] + b[:,-1]*(1 + dlm1/dlm2)
      if extrap:
         step0 = np.less(u, knots[0])*1.0
         bb[:,0] = bb[:,0] + step0*(1 + dl1/dl2 - degree/dl2*(u-knots[0]))
         bb[:,1] = bb[:,1] + step0*(-dl1/dl2 + degree/dl2*(u-knots[0]))
         step1 = np.greater(u, knots[-1])*1.0
         bb[:,-2] = bb[:,-2] - step1*(dlm1/dlm2 + degree/dlm2*(u-knots[-1]))
         bb[:,-1] = bb[:,-1] + step1*(1 + dlm1/dlm2 + degree/dlm2*(u-knots[-1]))
      return bb

   return b

if __name__ == "__main__":
   # Control points
   cv = np.array([[50.,  25., 0.],
                  [59.,  12., 0.],
                  [50.,  10., 0.],
                  [57.,   2., 0.],
                  [40.,   4., 0.],
                  [40.,   14., 0.]])
   
   n = 10 ** 6
   degree = 3  # Curve degree
   points_scipy = scipy_bspline(cv, n, degree)
   
   b = bspline_basis(len(cv), n, degree)
   points_basis = np.dot(b, cv)  
   print(np.allclose(points_basis, points_scipy))
