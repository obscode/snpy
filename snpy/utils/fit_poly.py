'''Fit 1D and 2D polynomials of order k to a set of data points.  Actually, they're
McLaren series in 1 and 2 variables.

Written by Chris Burns
Obfuscated by Dan Kelson.
'''
from numpy  import *
from numpy.linalg import svd

def divz(x, y=1, repl=0.0,out=float32, tol=0):
   if len(shape(y)) or len(shape(x)):
      if len(shape(y)):  bad = less_equal(absolute(y), tol)
      else:  band = ones(x.shape)*(abs(y)<=tol)
      not_bad = 1-bad
      numer = (x*not_bad)
      denom = (y+bad)
      a = (repl*bad).astype(out)
      b = (numer/denom).astype(out)
      return(a + b).astype(out)
   else:
      if abs(y) <= tol:  return(repl)
      else:  return(x/y)

def fitsvd(A, b):
   decomp = svd(A, full_matrices=False)
   sol1 = transpose(decomp[2])
   sol2 = divz(1.0,decomp[1], tol=1e-10)
   sol3 = dot(transpose(decomp[0]),b)
   if any(sol3):
      solr = (sol2*sol3)
      soll = dot(sol1,solr)
      err = sqrt(sum(power(sol1*sol2, 2), axis=0))
   else:
      soll = zeros(sol3.shape)
      err = zeros(sol3.shape)
   return soll,err


fac = lambda N: (not N and 1) or (fac(N-1)*N)

def fitpoly(x, y, w, k=1, x0=0):
   '''Fit a 1D McLaren series of degree k to some data (x,y) with weights w, 
   centered on the point x=x0.  Returns (coeff,err) which are the 
   coefficients of the series:  f(x0), f'(x0), f''(x0), ... and the 
   associated errors.'''

   w2 = sqrt(w)
   A = [1.0/fac(i)*power(x-x0, i)*w2 for i in range(k+1)]
   A = transpose(array(A))
   b = y*w2

   soll,err = fitsvd(A,b)
   return(soll, err)

def poly(x, x0, soll):
   '''Compute the polynomial from the solution.'''
   y = sum(array([1.0/fac(i)*power(x*1.0-x0, i)*soll[i] for i in range(len(soll))]), axis=0)
   return(y)


def fit2Dpoly(x, y, z, w, k=1, x0=0, y0=0):
   '''Fit a 2D McLaren series of degree k to some data (x,y,z) with weights w,
   centered on the point (x=x0, y=y0).  Returns the tuple (coeff, err).  Coeff
   are the coefficients of the series in the following order.  Let f[i,j](x0,y0)
   be the ith partial derivative of f w.r.t. x and jth partial derivative of f
   w.r.t y evaluated at (x0,y0).  The Coeff is
      f[0,0](x0,y0), f[1,0](x0,y0), ..., f[k,0](x0,y0),
      f[0,1](x0,y0), f[1,1](x0,y0), ..., f[k,1](x0,y0),
      ...
      f[0,k](x0,y0), f[1,k](x0,y0), ..., f[k,k](x0,y0)'''

   xbasis = [power(x-x0,i)/fac(i) for i in range(k+1)]
   ybasis = [power(y-y0,j)/fac(j) for j in range(k+1)]
   A = [xbasis[i]*ybasis[j] for j in range(k+1) for i in range(k+1) if i+j <= k]
   soll,err = fitsvd(transpose(A)*w[::,NewAxis],z*w)
   return (soll, err)

def poly2D(x, y, x0, y0, soll):
   k = (sqrt(1+4*len(soll))-1)/2
   xbasis = [power(x-x0,i)/fac(i) for i in range(k+1)]
   ybasis = [power(y-x0,i)/fac(i) for i in range(k+1)]
   ret = sum([xbasis[i]*ybasis[j] for j in range(k+1) \
                                  for i in range(k+1) if i+j <= k])
   return(ret)

def legendre(x, y, u):
   '''Given the length of x and y, determine the legendre interpolating
   polynomial and evaluate at u'''
   if len(x) != len(y):
      raise ValueError("Error:  x and y must be of equal size")
   scalar = 0
   if len(shape(u)) == 0:
      scalar = 1
      u = array([u])

   Dx = -x[NewAxis,:] + u[:,NewAxis]
   Dxx = x[:,NewAxis] - x[NewAxis,:] + identity(len(x))

   ids = arange(len(x))
   numerator = array([product(compress(not_equal(ids, i), Dx), axis=1) \
         for i in range(len(x))])
   denominator = array([product(Dxx[:,i]) for i in range(len(x))])
   result = sum(numerator*y[:,NewAxis]/denominator[:,NewAxis])
   if scalar:
      return(result[0])
   else:
      return(result)
