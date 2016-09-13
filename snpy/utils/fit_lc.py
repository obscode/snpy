'''
This module fits a parameterized function for a SNIa light-curve with one
or two peaks.  Taken from M. Stritzinger's PhD thesis, which was adapted
from Contardo, G., Leibundgut, B., & Vacca, W. D. 2000, A&A, 359, 876.
'''

from scipy.optimize import leastsq,brentq
from numpy import *


def Ialcn(t, par, n):
   m0,gamma,tau,theta,t0 = par[0:5]
   num = m0 + gamma*(t - t0)
   denom = 1 - exp((tau-t)/theta**2)
   for i in range(n):
      ti,gi,sigi = par[4+i*3:7+i*3]
      num -= gi**2*exp(-0.5*power(t-ti,2)/sigi**2)
   return(num/denom)

def d_Ialcn_dt(t, par, n):
   m0,gamma,tau,theta = par[0:4]
   num = gamma
   ex = exp((tau - t)/theta**2)
   for i in range(n):
      ti,gi,sigi = par[4+i*3:7+i*3]
      num -= gi*exp(-0.5*power(t-ti,2)/sigi**2)/sigi**2*(t-ti)
   num -= ex/theta**2*Ialcn(t, par,n)
   return (num/(1-ex))

def d_Ialcn_dp(par, t, y, dy, n):
   m0,gamma,tau,theta,t0 = par[0:5]
   ex = exp((tau - t)/theta**2)
   jac = zeros((len(par), len(t)), dtype=float32)
   jac[0,:] = 1.0/(1 - ex)
   jac[1,:] = jac[0,:]*(t - t0)
   jac[2,:] = Ialcn(t, par, n)/theta**2/(1-ex)*ex
   jac[3,:] = -2*jac[2,:]*(tau-t)/theta
   for i in range(n):
      ti,gi,sigi = par[4+i*3:7+i*3]
      Gi = gi*exp(-0.5*power(t-ti,2)/sigi**2)
      jac[4+i*3,:] = jac[0,:]*(Gi/sigi**2*(t - ti) - gamma)
      jac[5+i*3,:] = jac[0,:]*Gi/gi
      jac[6+i*3,:] = jac[0,:]*Gi/sigi**3*power(t-ti,2)
   return jac


def wrap_Ialcn(p, x, y, dy, n):
   return((Ialcn(x, p, n)-y)/dy)


def guess_parsn(t, n, mag, Tmax=None):
   p = [0]*(7+3*(n-1))
   if Tmax is None:
      id = argmin(mag)
      p[0] = mag[id]
      p[4] = t[id]
   else:
      p[0] = mag[argmin(absolute(t-Tmax))]
      p[4] = Tmax

   p[1] = (p[0] - mag[-1])/(p[4] - t[-1])
   p[2] = p[4] - 100
   p[3] = 3.0
   p[5] = -1.0
   p[6] = 10.0
   dt = (t[-1] - p[4])/n
   for i in range(n-1):
      p[7+i*3] = p[4]+(i+1)*dt
      p[8+i*3] = -1.0
      p[9+i*3] = 10.0
   return p


def fit_lc(t, mag, e_mag, ngauss=1, maxiter=10000, p0=None, Tmax=None):
   '''Fit a light-curve to the parameterized model.  t = time, mag = magnitudes,
   e_mag = error in magnitudes.  If ngauss=1, fit a single-peaked LC, if ngauss=2,
   fit a 2-peaked one.'''

   if p0 is None:
      p0 = guess_parsn(t, ngauss, mag, Tmax)
   par,cov,info,mesg,ier = leastsq(wrap_Ialcn, p0, 
      args=(t, mag, e_mag, ngauss), full_output=1, maxfev=maxiter,
      Dfun=d_Ialcn_dp, col_deriv=1)

   return(par,cov,info,mesg,ier)
