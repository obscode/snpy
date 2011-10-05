## Automatically adapted for numpy.oldnumeric Feb 04, 2009 by ipython

import numpy.oldnumeric as num
import numpy.oldnumeric.mlab as MLab

def divz(x,y=1,repl=0.0,out=num.Float32):
   if len(num.shape(y)) or len(num.shape(x)):
      if len(num.shape(y)): bad = num.equal(y,0.0)
      else: bad = num.ones(x.shape)*(y==0.0)
      not_bad = 1-bad
      numer = (x*not_bad)
      denom = (y+bad)
      a = (repl*bad).astype(out)
      b = (numer/denom).astype(out)
      return (a + b).astype(out)
   else:
      if y == 0: return repl
      else: return x/y

def mean( x,axis=0):
  if len(x): return MLab.mean(num.asarray(x),axis)
  else: return 0.0

def average( x,axis=0):
  return mean(x,axis=axis)

def rms( x,y=None):
  if len(x) > 2:
    if y == None: y=num.average(x)
    return num.sqrt(mean(num.power(num.subtract(x,y),2)))
  else: return 0.0

def sigma(x, y=None):
   if len(x) > 1:
      if y == None:  y=num.average(x)
      N = len(x)
      return num.sqrt(sum(num.power(num.subtract(x,y),2))/(N-1))
   else: return -1.0

def median( x, sorted=0, nullval = 0.0, axis=0):
  x = num.asarray(x)
  n = x.shape[axis]
  median=0.0
  if n:
     if not sorted: x = num.sort(x,axis)
     if n%2:
        if axis == 0: median = x[(n-1)/2]
        else: median = num.take(x,[(n-1)/2],axis)
     else:
        median = num.add.reduce(num.take(x,[n/2-1,n/2],axis),axis)/2.0
     return median
  else: return nullval


def bwt( x, iter=3, sorted=0):
 x=num.asarray(x)
 ns = num.sqrt(len(x))
 if len(x)>0: M = median(x,sorted=sorted)
 else: M=0.0
 if len(x) <= 2:
  S = 0*M
 else:
  MAD = median(abs(x-M))
  if not MAD: MAD = rms(x-M)*0.6745
  if MAD:
     for i in range(iter):
        u6 = divz(x-M,MAD)/6.0
        omuu62 = num.power(1-num.power(u6,2),2)
        g6 = num.less(num.absolute(u6),1.0)
        M = M + num.add.reduce(g6*(x-M)*omuu62)/num.add.reduce(g6*omuu62)
        u9 = divz(x-M,MAD)/9.0
        omuu9 = 1-num.power(u9,2)
        omuu59 = 1-5*num.power(u9,2)
        g9 = num.less(abs(u9),1.0)
        S = ns * divz(num.sqrt(num.add.reduce( \
              g9*num.power(x-M,2)*num.power(omuu9,4))), 
              num.absolute(num.add.reduce(g9*omuu9*omuu59)))
        MAD = 0.6745*S
  else: M,S = mean(x),0.0
 return [M,S]

def cum_prob(prob):
   # Convert a probability density into a cumulative probability density
   shape = prob.shape
   sids = num.argsort(num.ravel(prob))[::-1]
   sprob = num.take(num.ravel(prob), sids)
   scprob = num.cumsum(sprob)
   cprob = 0.0*prob
   num.put(cprob, sids, scprob)
   cprob.shape = shape
   return(cprob)

