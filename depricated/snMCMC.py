'''This is a module that contains a class derived from the pymc MCMC
class. It takes a SN object and generates a MCMC sampler.'''
import pymc
import numpy as np

debug=False

class MCMC_generator(pymc.MCMC):

   def __init__(self, snobj, filters=None, inc_var=False, **args):
      '''Create an MCMC sampler based on a sn object. The specified filters
      are fit using the model that is currently selected. Uniform
      priors are assumed for the parameters unless overridden by assigning
      pymc Stochastics through **args.'''

      self.sn = snobj
      if filters is None:
         filters = list(self.sn.data.keys())


      self.model = snobj.model
      self.model.args = {}
      self.model._fbands = filters
      self.model.setup()
      params = []
      paramnames = list(self.model.parameters.keys())
      # First, setup stochastics for our parameters
      for param in paramnames:
         if param in args:
            params.append(args[param])
            del args[param]
            continue
         if param == 'dm15':
            params.append(pymc.Uniform('dm15', 0.7, 2.0))
         elif param == 'st':
            params.append(pymc.Uniform('st', 0.25, 1.22))
         elif param == 'Tmax':
            t0 = min([self.sn.data[f].MJD.min() for f in self.sn.data])
            t1 = max([self.sn.data[f].MJD.max() for f in self.sn.data])
            params.append(pymc.Uniform('Tmax', t0-30, t1+30))
         elif param == 'EBVhost':
            params.append(pymc.Uniform('EBVhost', 0, 10.))
         elif param == 'DM':
            params.append(pymc.Uniform('DM', 0, 100))
         elif param.find('max') > 0:
            params.append(pymc.Uniform(str(param), 10., 30.))
         else:
            raise AttributeError("Error, parameter %s not recognized. Update MCMC package" % (param))
         if self.model.parameters[param] is None:
            params[-1].value = self.model.guess(param)
         else:
            params[-1].value = self.model.parameters[param]
      params = pymc.Container(params)

      # now setup intrinsic variances for each filter
      if inc_var:
         vars = pymc.InverseGamma('taus', alpha=0.5, beta=0.1**2, 
               value=np.random.uniform(0,0.1**2,size=len(filters)))
      else:
         vars = np.array([0.0]*len(filters))

      # The data stochastic that maps parameters to observations
      @pymc.data
      @pymc.stochastic
      def model(params=params, vars=vars, paramnames=paramnames, filters=filters,
            value=1.0):
         # Set the parameters in the model
         for i,param in enumerate(paramnames):
            if debug:
               print("setting ",param, " to ",params[i])
            self.model.parameters[param] = params[i]

         logp = 0
         numpts = 0
         for i,f in enumerate(filters):
            mod,err,mask = self.model(f, self.sn.data[f].MJD)
            m = mask*self.sn.data[f].mask
            if not np.any(m):
               continue
            numpts += np.sum(m)
            tau = np.power(vars[i] + np.power(self.sn.data[f].e_mag,2),-1)
            logp += pymc.normal_like(self.sn.data[f].mag[m],mod[m],tau[m])
         #if numpts < len(paramnames):
         #   return -np.inf
         return logp

      pymc.MCMC.__init__(self, locals(), **args)

      # Setup the step methods
      # 1) params will be AdaptiveMetropolis, so we need to setup initial
      #    scales. If the model has been fit, use error, otherwise guess.
      def_scales = {'Tmax':0.5**2,
                    'st':0.001**2,
                    'dm15':0.001**2,
                    'max':0.01**2,
                    'DM':0.01**2,
                    'EBVhost':0.01**2}
      scales = {}
      for i,par in enumerate(self.paramnames):
         if par in self.model.errors and self.model.errors[par] > 0:
            scales[self.params[i]] = self.model.errors[par]
         else:
            if par in def_scales:
               scales[self.params[i]] = def_scales[par]
            elif par[0] == "T" and par[-3:] == "max":
               scales[self.params[i]] = def_scales['Tmax']
            elif par[-3:] == "max":
               scales[self.params[i]] = def_scales['max']
            else:
               scales[self.params[i]] = self.params[i].value/10.
      self.use_step_method(pymc.AdaptiveMetropolis, self.params, scales=scales,
            delay=1000, interval=1000)

      if inc_var:
         self.use_step_method(pymc.AdaptiveMetropolis, [self.vars], 
               scales={self.vars:self.vars.value*0+0.005**2})
