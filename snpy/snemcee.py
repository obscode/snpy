'''This is a module that tries to use emcee to solve SNooPy models.
The SN object should have all the necessary ingredients. All that is
left is to define prior probabilities.'''

import emcee
import numpy as np
from scipy.optimize import minimize
import types,os

gconst = -0.5*np.log(2*np.pi)

def builtin_priors(x, st):
   '''Some built-in priors that are simple strings that we parse.'''
   if st[0] not in ['U','G','E']:
      raise ValueError("I don't understand the prior code %s" % st)
   if st[0] == 'U':
      '''uniform prior:  U,l,u'''
      u,l = list(map(float, st.split(',')[1:]))
      if not u < x < l:
         return -np.inf
      return 0
   elif st[0] == 'G':
      '''Gaussian prior, G(mu,sigma) = 1/sigma/sqrt(2*pi)*exp((x-mu)^2/2/sigma^2)'''
      mu,sigma = list(map(float, st.split(',')[1:]))
      return gconst - 0.5*np.power(x - mu,2)/sigma**2 - np.log(sigma)
   elif st[0] == 'E':
      '''exponential prior, E(x,tau) = 1/tau*exp(-x/tau).'''
      tau = float(st.split(',')[1])
      return -np.log(tau) - x/tau

def vecgauss(x, mu, sig):
   '''A simple vector-based Gaussian with mean mu and std sig'''
   return np.sum(gconst - 0.5*np.power(x-mu,2)/np.power(sig,2) - np.log(sig))

def guess(varinfo, snobj):
   '''Get starting values from the fitter.'''
   p = np.zeros((varinfo['Nvar'],))
   ep = np.zeros((varinfo['Nvar'],))
   for var in varinfo['free']:
      if var in snobj.model.nparameters:
         p[varinfo[var]['index']] = snobj.model.nparameters[var]
         ep[varinfo[var]['index']] = snobj.model.enparameters[var]
      else:
         if snobj.model.parameters[var] is None:
            raise ValueError("model parameters not set, run initial fit() first")
         p[varinfo[var]['index']] = snobj.model.parameters[var]
         #ep[varinfo[var]['index']] = max(0.001,snobj.model.errors[var])
         ep[varinfo[var]['index']] = 1e-4
   return p,ep



def setup_varinfo(snobj, args):
   '''Given a sn object and its associated model, setup the varinfo.'''
   varinfo = {}
   varinfo['varlist'] = list(snobj.model.parameters.keys())
   # Nuissance parameters
   varinfo['varlist'] = varinfo['varlist'] + list(snobj.model.nparameters.keys())
   varinfo['fitflux'] = args.get('fitflux',True)
   i  = 0
   varinfo['free'] = []
   for var in varinfo['varlist']:
      varinfo[var] = {}
      if var in args:
         if type(args[var]) is bytes:
            varinfo[var]['fixed'] = False
            varinfo[var]['index'] = i
            varinfo['free'].append(var)
            i += 1
            varinfo[var]['prior'] = args[var]
            varinfo[var]['prior_type'] = 'builtin'
         elif type(args[var]) in np.ScalarType:
            varinfo[var]['value'] = args[var]
            varinfo[var]['fixed'] = True
         elif type(args[var]) is types.FunctionType:
            varinfo[var]['fixed'] = False
            varinfo[var]['index'] = i
            varinfo['free'].append(var)
            i += 1
            varinfo[var]['prior'] = args[var]
            varinfo[var]['prior_type'] = 'function'
      else:
         if var in snobj.model.nparameters:
            if snobj.model.enparameters[var] is not None:
               varinfo[var]['fixed'] = False
               varinfo[var]['prior_type'] = 'nuissance'
               varinfo[var]['value'] = snobj.model.nparameters[var]
               varinfo[var]['std'] = snobj.model.enparameters[var]
               if len(np.shape(varinfo[var]['value'])) == 1:
                  varinfo[var]['index'] = slice(i,i+varinfo[var]['value'].shape[0])
                  i += varinfo[var]['value'].shape[0]
               else:
                  varinfo[var]['index'] = i
                  i += 1
               varinfo['free'].append(var)
            else:
               varinfo[var]['fixed'] = True
               varinfo[var]['value'] = snobj.model.nparameters[var]
         else:
            varinfo[var]['fixed'] = False
            varinfo[var]['index'] = i
            varinfo['free'].append(var)
            i += 1
            varinfo[var]['prior_type'] = 'model'
   varinfo['Nvar'] = i
   return varinfo


def lnprior(p, varinfo, snobj):
   lp = 0
   for var in varinfo['free']:
      id = varinfo[var]['index']
      val = p[id]
      if varinfo[var]['prior_type'] == 'function':
         lp += varinfo[var]['prior'](val)
      elif varinfo[var]['prior_type'] == 'builtin':
         lp += builtin_priors(val, varinfo[var]['prior'])
      elif varinfo[var]['prior_type'] == 'model':
         lp += snobj.model.prior(var,val)
      elif varinfo[var]['prior_type'] == 'nuissance':
         lp += vecgauss(p[id], varinfo[var]['value'], varinfo[var]['std'])
   return lp

def lnlike(p, varinfo, snobj, bands):
   # first, assign all variables to the model:
   for id,var in enumerate(varinfo['varlist']):
      if varinfo[var]['fixed']:
         if var in snobj.model.parameters:
            snobj.model.parameters[var] = varinfo[var]['value']
         else:
            snobj.model.nparameters[var] = varinfo[var]['value']
      else:
         val = p[varinfo[var]['index']]
         if var in snobj.model.parameters:
            snobj.model.parameters[var] = val
         else:
            snobj.model.nparameters[var] = p[varinfo[var]['index']]
   lp = 0
   for band in bands:
      mod,err,mask = snobj.model.__call__(band, snobj.data[band].MJD)
      fitflux = varinfo['fitflux']
      if fitflux:
         if snobj.model.model_in_mags:
            f = np.power(10, -0.4*(mod - snobj.data[band].filter.zp))
            cov_f = np.power(f*err/1.0857,2)
         else:
            f = mod
            cov_f = np.power(err, 2)
      else:
         if snobj.model.model_in_mags:
            f = mod
            cov_f = np.power(err, 2)
         else:
            f = -2.5*log10(mod) + snobj.data[band].filter.zp
            cov_f = np.power(err/mod*1.0857,2)

      m = mask*snobj.data[band].mask
      if not np.any(m):
         # We're outside the support of the data
         return -np.inf
      N = sum(m)
      X = snobj.data[band].flux[m] - f[m]
      #if not np.all(m):
      #   ids = np.nonzero(-m)[0]
      #   thiscovar = np.delete(np.delete(snobj.bcovar[band],ids,axis=0), 
      #         ids, axis=1)
      #else:
      #   thiscovar = snobj.bcovar[band]
      #detcovar = np.linalg.det(thiscovar)
      #invcovar = np.linalg.inv(thiscovar)
      #lp = lp - 0.5*np.log(2*np.pi**N*detcovar) -\
      #      0.5*np.dot(X, np.dot(invcovar,X))
      denom = cov_f[m] + np.power(snobj.data[band].e_flux[m],2)
      lp = lp - 0.5*np.sum(np.power(X,2)/denom + \
            np.log(denom) + np.log(2*np.pi))
   return lp

def lnprob(p, varinfo, snobj, bands):

   lprior = lnprior(p, varinfo, snobj)
   if not np.isfinite(lprior):
      return -np.inf
   return lprior + lnlike(p, varinfo, snobj, bands)
            #raise RuntimeError, "Model must be in mags"
            #raise RuntimeError, "Model must be in mags"


def generateSampler(snobj, bands, nwalkers, threads=1, tracefile=None, **args):
   '''Generate an emcee sampler from the sn object [snobj] and its
   associated model (chosen with snobj.choose_model). You must set the
   number of walkers (see emcee documentation).  You can control
   the priors of the model by passing them as arguments. For example,
   using Tmax='G,1000,10' would use a Gaussian prior with mean 1000
   and standard deviation 10. You can also set any parameter to a
   constant value. Lastly, you can set a parameter equal to a function
   that takes a single argument and returns the log-probability as
   a prior. 
   
   This function returns:  sampler,p0
   where sampler is an emcee sampler, and p0 is [nwalkers] starting
   points.'''

   tp0 = None
   if tracefile is not None:
      if os.path.isfile(tracefile):
         tpars = []
         f = open(tracefile)
         line = f.readline()
         Nwalkers = 50
         while line[0] == '#':
            if line.find('Col') > 0:
               tpars.append(line.split()[-1])
            elif line.find('Nwalkers') >= 0:
               Nwalkers = int(line.split()[-1])
            line = f.readline()
         f.close()
         data = np.loadtxt(tracefile)
         Niter = data.shape[0]/Nwalkers
         endids = [(i+1)*Niter-1 for i in range(Nwalkers)]
         tp0 = [data[ids,:] for ids in endids]

   if not snobj.model._fbands:
      raise ValueError("You need to do an initial fit to the SN first")
   vinfo = setup_varinfo(snobj, args)
   p,ep = guess(vinfo, snobj)

   # Find the ML:
   #nll = lambda *args: -lnlike(*args)
   #result = minimize(nll, p, args=(vinfo,snobj,bands))
   #p = result["x"]
   ndim = p.shape[0]
   p0 = [];  fail=0
   while len(p0) < nwalkers and fail < 1000:
      pp = p + ep*np.random.randn(ndim)
      if not np.isinf(lnprior(pp, vinfo, snobj)):
         p0.append(pp)
      else:
         fail += 1
   if len(p0) < nwalkers:
      raise RuntimeError("Could not establish an initial set of MCMC walkers.\n" +\
            "Make sure your priors are consistent with your intial fit solution")
   if tp0 is not None:
      for i in range(len(p0)):
         for ii,par in enumerate(tpars):
            j = vinfo[par]['index']
            p0[i][j] = tp0[i][ii]
   sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(vinfo, snobj, bands),
         threads=threads)
   return sampler,vinfo,p0
