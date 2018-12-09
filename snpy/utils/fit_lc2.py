'''
This module fits a parameterized function for a SNIa light-curve with one
or two peaks.  Taken from M. Stritzinger's PhD thesis, which was adapted
from Contardo, G., Leibundgut, B., & Vacca, W. D. 2000, A&A, 359, 876.
'''

from scipy.optimize import leastsq,brentq
from numpy import *

def Ialc1(t, m0, gamma, t0, g0, sig0, theta):
   '''The parametric function with one maximum we want to fit.'''
   num = m0 + gamma*(t - t0)
   G0 = -g0**2*exp(-0.5*power(t - t0,2)/sig0**2)
   num += G0
   denom = where(less(t-t0, 0), cos((t - t0)/theta), 1.)
   return(num/denom)

def d_Ialc1_dt(t, m0, gamma, t0, g0, sig0, tau, theta):
   '''First derivative of Ialc1 wrt t.'''
   num = gamma
   G0 = -g0**2*exp(-0.5*power(t - t0,2)/sig0**2)
   ex = exp((tau - t)/theta**2)
   num -= G0/sig0**2*(t - t0)
   num -= ex/theta**2*Ialc1(t, m0, gamma, t0, g0, sig0, tau, theta)
   return(num/(1-ex))

def d_Ialc1_dp(t, p):
   m0,gamma,t0,g0,sig0,tau,theta = p
   G0 = g0*exp(-0.5*power(t - t0,2)/sig0**2)
   ex = exp((tau - t)/theta**2)
   jac = zeros((len(p), len(t)), dtype=float64)
   jac[0,:] = 1.0/(1 - ex)
   jac[1,:] = jac[0,:]*(t - t0)
   jac[2,:] = jac[0,:]*(G0/sig0**2*(t - t0) - gamma)
   jac[3,:] = jac[0,:]*G0/g0
   jac[4,:] = jac[0,:]*G0/sig0**3*power(t-t0,2)
   jac[5,:] = Ialc1(t, *p)/theta**2/(1-ex)*ex
   jac[6,:] = -2*jac[5,:]*(tau - t)/theta**3
   return jac


def Ialc2(t, m0, gamma, t0, g0, sig0, t1, g1, sig1, theta):
   '''The parametric function with two maxima we want to fit.'''
   num = m0 + gamma*(t - t0)
   G0 = -g0**2*exp(-0.5*power(t - t0,2)/sig0**2)
   G1 = -g1**2*exp(-0.5*power(t - t1,2)/sig1**2)
   num += G0 + G1
   denom = where(less(t-t0, 0), cos((t - t0)/theta), 1.)
   return(num/denom)

def d_Ialc2_dt(t, m0, gamma, t0, g0, sig0, t1, g1, sig1, tau, theta):
   '''First derivative of Ialc2 wrt t.'''
   num = gamma
   G0 = g0*exp(-0.5*power(t - t0,2)/sig0**2)
   G1 = g1*exp(-0.5*power(t - t1,2)/sig1**2)
   ex = exp((tau - t)/theta**2)
   num -= G0/sig0**2*(t - t0) + G1/sig1**2*(t - t1)
   num -= ex/theta**2*Ialc2(t, m0, gamma, t0, g0, sig0, t1, g1, sig1, tau, theta)
   return(num/(1-ex))

def d_Ialc2_dp(p,t,mag,e_mag):
   m0,gamma,t0,g0,sig0,t1,g1,sig1,tau,theta = p
   G0 = g0*exp(-0.5*power(t - t0,2)/sig0**2)
   G1 = g1*exp(-0.5*power(t - t1,2)/sig1**2)
   ex = exp((tau - t)/theta**2)
   jac = zeros((len(p), len(t)), dtype=float32)
   jac[0,:] = 1.0/(1 - ex)
   jac[1,:] = jac[0,:]*(t - t0)
   jac[2,:] = jac[0,:]*(G0/sig0**2*(t - t0) - gamma)
   jac[3,:] = jac[0,:]*G0/g0
   jac[4,:] = jac[0,:]*G0/sig0**3*power(t-t0,2)
   jac[5,:] = jac[0,:]*(G1/sig1**2*(t - t1))
   jac[6,:] = jac[0,:]*G1/g1
   jac[7,:] = jac[0,:]*G1/sig1**3*power(t-t1,2)
   jac[8,:] = Ialc2(t, *p)/theta**2/(1-ex)*ex
   jac[9,:] = -2*jac[8,:]*(tau - t)/theta
   return jac


def wrap_Ialc1(p, x, y, dy):
   return((y - Ialc1(x, *p))/dy)

def wrap_Ialc2(p, x, y, dy):
   return((y - Ialc2(x, *p))/dy)

def guess_pars1(t, mag, Tmax=None):
   '''Based on the input light-curve, guess the likely paramters.'''
   p = [0]*6

   if Tmax is None:
      id = argmin(mag)
      p[0] = mag[id]
      p[2] = t[id]
   else:
      p[0] = mag[argmin(absolute(t - Tmax))]
      p[2] = Tmax

   p[1] = (p[0] - mag[-1])/(p[2] - t[-1])
   p[3] = -1.0
   p[4] = 10.0
   p[5] = 3.
   return p

def guess_pars2(t, mag, Tmax=None):
   '''Based on the input light-curve, guess the likely paramters.'''
   p = [0]*9

   if Tmax is None:
      id = argmin(mag)
      p[0] = mag[id]
      p[2] = t[id]
   else:
      p[0] = mag[argmin(absolute(t - Tmax))]
      p[2] = Tmax
   p[1] = (p[0] - mag[-1])/(p[2] - t[-1])
   p[3] = -1.0
   p[4] = 10.0
   p[5] = p[2] + 20.
   p[6] = -1.0
   p[7] = 10.0
   p[8] = 3.
   return p

def fit_lc(t, mag, e_mag, ngauss=1, maxiter=10000, p0=None, Tmax=None):
   '''Fit a light-curve to the parameterized model.  t = time, mag = magnitudes,
   e_mag = error in magnitudes.  If ngauss=1, fit a single-peaked LC, if ngauss=2,
   fit a 2-peaked one.'''

   if ngauss==1:
      if p0 is None:
         p0 = guess_pars1(t, mag, Tmax)
      par,cov,info,mesg,ier = leastsq(wrap_Ialc1, p0, args=(t, mag, e_mag),
            full_output=1, maxfev=maxiter)
   else:
      if p0 is None:
         p0 = guess_pars2(t, mag, Tmax)
      par,cov,info,mesg,ier = leastsq(wrap_Ialc2, p0, args=(t, mag, e_mag),
            full_output=1, maxfev=maxiter)# , Dfun=d_Ialc2_dp, col_deriv=1)
   return(par,cov,info,mesg,ier)
