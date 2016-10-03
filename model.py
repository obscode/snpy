'''Model.py:  a module that defines the SN models to be fit by SNOOPY.

A base class (Model) is defined to handle most of the heavy-lifting and boiler
plate around scipy.optimize.leastsq. Defining a new model is done by sub-
classing Model and overriding the member functions.

New:  Add an optional [decline_param] to choose between a dm15 model and stretch
      (st)  model'''
import os,string
from snpy import ubertemp
from snpy import kcorr
from snpy.utils import redlaw
from numpy.linalg import cholesky
from scipy import stats
from scipy.optimize import leastsq
from scipy.optimize import brent
import scipy.interpolate
from numpy import *
#from numpy import median, bool, diag
from numpy.linalg import inv
import pickle

Ia_w,Ia_f = kcorr.get_SED(0, 'H3')
gconst = -0.5*log(2*pi)

debug = 0
base = os.path.dirname(globals()['__file__'])

it_sigma = 1.0

def _quniform_prior(p,pmin,pmax,sigma):
   '''A quasi-uniform prior (uniform with Gaussian tails).'''
   norm = 1.0/((pmax - pmin)+sqrt(2*pi)*sigma)
   if pmin < p < pmax:  return norm
   if p <= pmin:  return norm*exp(-0.5*(p-pmin)**2/sigma**2)
   if p >= pmax:  return norm*exp(-0.5*(p-pmax)**2/sigma**2)


class model:
   '''The base class for SNooPy light-curve models. It contains the parameters
   to be solved and its __call__ function returns the model for a filter's
   data. It has several convenience member functions for figuring out things
   like K-corrections and reddening coefficients.'''

   model_in_mags = 1

   def __init__(self, parent):
      '''Setup the model.
      Args:
         parent (snpy.sn instance): The parent sn class.
      '''
      self.parent = parent       # make a link to the SN object
      self.parent.model = self   # make a link to the model
      self.parameters = {}    # dictionary of parameters of the model
      self.nparameters = {}    # dictionary of nuissance parameters of the model
      self.errors = {}        # dictionary of errors on the parameters
      self.enparameters = {}    # dictionary of nuissance parameter errors
      self.fixed = []         # list of parameters that are held fixed
      self.do_kcorr = 1       # Do we perform initial k-corrections?
      self.rbs = []           # List of rest-frame filters this model supports
      self._fbands = []
      self.args = {}
      self.MWRobs = {}

   def guess(self, param):
      '''A function that is run and picks initial guesses for
      the parameter.
      
      Args:
         param (str):  the name of the parameter.
      Returns:
         Value: initial guess for the parameter.'''
      raise NotImplementedError('Derived class must overide')

   def prior(self, param, value):
      '''An optional prior on the parameters.  To be filled in
      based on implementation.  This should return the probability
      given the current value of a parameter.

      Args:
         param (str): The name of the parameter
         value (float): The current value of the parameter.
      Returns:
         prior (float):  the probability of parameter.'''
      return 1.0

   def setup(self):
      '''A function that does any initial setup (estimating k-correcions,
      etc)'''
      raise NotImplementedError('Derived class must overide')

   def kcorr(self, band, t):
      '''Convenience function to call back to the parent and get the 
      k-corrections.
      
      Args:
         band (str):  the name of the filter
         t (float array): The epochs (t - T(Bmax)) of the observations
      Returns:
         2-tuple (K,mask)
         
         * K (float array): K-corrections
         * mask (bool array): True where K-corrections are valid
      '''

      # If k-corrections are there, use them
      if band in self.parent.ks_tck:   
         tck = self.parent.ks_tck[band]
         K = scipy.interpolate.splev(t+self.Tmax, tck)
         k0 = scipy.interpolate.splev(tck[0][0], tck)
         k1 = scipy.interpolate.splev(tck[0][-1], tck)
         K = where(less(t+self.Tmax, tck[0][0]), k0, K)
         K = where(greater(t+self.Tmax, tck[0][-1]), k1, K)
         mids = argmin(absolute(t[:,newaxis]-self.parent.data[band].MJD[newaxis,:]+\
               self.Tmax), axis=1)
         # mask based on original mask and limits of Hsiao spectrum
         mask2 = self.parent.ks_mask[band][mids]*\
               greater_equal(t, -19)*less_equal(t, 70)
      else:
         K = 0*t
         mask2 = ones(t.shape, dtype=bool)
      return K,mask2

   def MWR(self, band, t):
      '''Determine the best :math:`R_\lambda` for the foreground MW extinction.
      
      Args:
         band (str):  The name of the filter
         t (float array): The epoch (t - T(Bmax)) of the observations.
         
      Returns:
         float array: :math:`R_\lambda` for the filter at times t.
      '''

      if band in self.parent.Robs:
         if type(self.parent.Robs[band]) is type(()):
            R = scipy.interpolate.splev(t+self.parent.Tmax, 
                  self.parent.Robs[band])
         else:
            R = self.parent.Robs[band]
      else:
         if band in self.MWRobs:
            R = self.MWRobs[band]
         else:
            self.MWRobs[band] = \
                  self.parent.data[band].filter.R(self.parent.Rv_gal)
            R = self.MWRobs[band]
      return R

   def _extra_error(self, parameters):
      return 0

   def fit(self, bands, epsfcn=0, **args):
      '''Fit the model with currently fixed and free parameters againts the
      set of bands [bands].  All other arguments are passed directly to
      the model() member function as optional arguments.  After running,
      the free parameters will be set to their fit values and self.C will
      contain the covariance matrix.
      
      Args:
         bands (list of str): The filters to fit
         epsfcn (float): see scipy.optmize.leastsq.
         
      Returns:
         None
         
      Effects:
         The member variables self.parameters and self.eparameters will
         be set with best-fit values. self.C will be updated with the 
         values from the covariance matrix.
      '''
      self.args = args.copy()
      if debug:  print "model.fit() called with args = ", args

      for b in bands:
         if b not in self.parent.data:
            raise ValueError, "band %s not defined in this SN object" % (b)

      self._fbands = bands

      # Setup the MW reddening
      self.setup()

      # build up the fixed variables
      self.fixed = []
      for key in args.keys():
         if key in self.parameters:
            # check that it is a valid floag
            try:
               testvalue = 1.0*args[key]
            except TypeError:
               raise ValueError, "You are trying to hold %s fixed, but with an illegal value or type: %s" % (key, str(args[key]))
            self.parameters[key] = args[key]
            self.fixed.append(key)
            del args[key]
      # build up list of free variables:
      self._free = [p for p in self.parameters.keys() if p not in self.fixed]
      pars = [self.parameters[p] for p in self._free]
      # Make initial parameter guesses if any are None
      for i in range(len(pars)):
         if pars[i] is None:
            pars[i] = self.guess(self._free[i])
            if debug:  print "Initial guess for %s = %f" % (self._free[i],pars[i])

      # get the error matrix or vector
      error = {}
      for band in bands:
         error[band] = self.parent.data[band].get_covar(flux=1)

      pars,C,self.info,self.mesg,self.ier = \
            leastsq(self._wrap_model, pars, (bands,error), full_output=1)
      if self.ier > 4:  print self.mesg

      self.chisquare = sum(power(self.info['fvec'], 2))
      self.dof = len(self.info['fvec']) - len(self._free)
      if self.dof < 1:
         print "Warning!  less than 1 degree of freedom!"
         self.rchisquare = self.chisquare
      else:
         self.rchisquare = self.chisquare/self.dof
      
      # update errors and covariance matrix
      if C is None:
         raise RuntimeError, "Error:  Covariance Matrix is singular.  " + \
               "Either two or more parameters are degenerate or the model" + \
               "has become insensitive to one or more parameters."
      C = C*self.rchisquare    # the trick for underestimated errors.
      self.C = {}
      for p in self.parameters:  self.errors[p] = self._extra_error(p)
      for i in range(len(self._free)):
         if C[i,i] < 0:
            print "Error: covariance matrix has negative diagonal element."
            print "       Error for %s not computed" % (self._free[i])
            self.errors[self._free[i]] = 0
         else:
            self.errors[self._free[i]] += sqrt(C[i,i])
         if self._free[i] not in self.C:  self.C[self._free[i]] = {}
         for j in range(len(self._free)):
            self.C[self._free[i]][self._free[j]] = C[i,j]

   def covar(self, band, t):
      return zeros((t.shape[0],t.shape[0]))

   def systematics(self):
      '''compute the systematic errors associated with the model paramters.
      
      Args:
         None
      
      Returns:
         dict:  a dictionary of systematic errors keyed by parameter
                name.  It therefore depends on the model being used.  Also
                see the specific model for any extra arguments.  If None
                is returned as a value, no systematic has been estimated
                for it.'''
      raise NotImplementedError('Derived class must overide')

   def get_max(self, bands, restframe=0, deredden=0):
      '''Get the maxima of the light-curves, given the current state of
      the model.

      Args:
         bands (list of str):  list of filters to find maxima for
         restframe (bool):  If True, apply K-corrections to maximum magnitudes
         deredden (bool):  If True, apply all known extinction corrections
                           to maximum magnitudes. Note: this will always
                           correct for Milky-Way extinction, but in some
                           models, host-galaxy extinction may be removed as
                           well.

      Returns: 
         4-tuple:  (Tmax,Mmax,eMmax,restbands):

         * Tmax (list of floats): Time of maxima for each filter
         * Mmax (list of floats): Maximum magnitudes for each filter
         * eMmax (list of floats): errors in maximum magnitudes
         * restbands (list of str): The rest-bands used to fit each observed
                                    filter
   '''
      raise NotImplementedError('Derived class must overide')

   def _wrap_model(self, pars, bands, error):
      resids_list = []
      sum_w = 0
      for i in range(len(pars)):
         self.parameters[self._free[i]] = pars[i]

      if debug:  
         print ">>> _wrap_model called with pars:"
         for i in range(len(self._free)):
            print "   %s:  %f" % (self._free[i],pars[i])
      for band in bands:
         if debug: print">>> calling model member function"
         mod,err,mask = self.__call__(band, self.parent.data[band].MJD)
         if not sometrue(mask):
            msg = "All weights for filter %s are zero." % band
            msg += " The fitter is in a part of parameter space where the model"
            msg += " is not valid or there is no useful data."
            raise RuntimeError, msg
         if self.model_in_mags:
            f = power(10, -0.4*(mod - self.parent.data[band].filter.zp))
            #ef = sqrt(power(err*mod/1.0857,2) + \
            #          power(self.parent.data[band].e_flux,2))
            cov_f = power(f*err/1.0857,2)
         else:
            f = mod
            cov_f = power(err,2)
         m = mask*self.parent.data[band].mask
         if len(shape(error[band])) == 1:
            # simple weight array
            W = m*1.0/sqrt(error[band]+cov_f)
            sum_w += sum(W)
            resids_list.append((f - self.parent.data[band].flux)*W)
         else:
            W = cholesky(inv(error[band]+diag(cov_f)))
            W = W*m[newaxis,:]*m[:,newaxis]
            sum_w += sum(diagonal(W))
            resids_list.append(dot(W, f - self.parent.data[band].flux))

      res = concatenate(resids_list)

      # now apply any priors:  chi2 = chi2 + 2*ln(p) -N/2log(2pi)-sum(log(sigi))
      if debug:  print "  weighted resids = ", resids_list
      #N = len(W[m])
      #extra = log(2*pi)-2/N*sum(log(W[m])) - 2*log(self.prior())/N
      #res = sqrt(power(res,2) + extra)
      return(res)

   def __getattr__(self, attr):
      if 'parameters' in self.__dict__:
         if attr in self.__dict__['parameters']:
            return self.__dict__['parameters'][attr]
      if attr in self.__dict__:
         return self.__dict__[attr]
      else:
         raise AttributeError, "Attribute %s not found" % (attr)

   def __setattr__(self, attr, value):
      if 'parameters' in self.__dict__:
         if attr in self.__dict__['parameters']:
            self.__dict__['parameters'][attr] = value
         else:
            self.__dict__[attr] = value
      else:
         self.__dict__[attr] = value



class EBV_model(model):
   '''This model fits any number of lightcurves with CSP uBVgriYJHK templates
   or Prieto BsVsRsIs templates.  The parameters you can fit:

   - dm15 (decline rate)
   - Tmax (time of peak B maximum)
   - DM   (distance modulus)
   - EBVhost  (host galaxy extinction)

   The model is constructed by assuming a peak B absolute magnitude  and B-X
   colors based on the current value of dm15.  The colors are from Folatelli
   et al. (2010), as are the calibration of Bmax vs dm15.  For the latter,
   there are 6 calibrations, based on the sample used to make the fit.  The
   default is 6 (best observed, excluding heavily extincted SNe), but you can
   choose a different calibration by setting that argument in the fit() call.
   Aside from the instrinsic colors, a global extinction parameter EBVhost
   is applied to each light-curve, as well as Milky way extinction from 
   the SN object's EBVgal.  The value of R_V for the host galaxy is
   not a parameter, but is controled by the choice of calibration in order to
   remain consistent with Folatelli et al. (2009).  The R_V for the galactic
   extinction is taken from the SN object (default 3.1).'''

   def __init__(self, parent, stype='dm15'):

      if stype != 'dm15':
         raise ValueError, "This model only supports the dm15 parameter"
      model.__init__(self, parent)
      self.rbs = ['u','B','V','g','r','i','Y','J','H','K','Bs','Vs','Rs','Is',
            'J_K','H_K']
      self.parameters = {'DM':None, 'dm15':None, 'EBVhost':None, 'Tmax':None}
      self.errors = {'DM':0, 'dm15':0, 'EBVhost':0, 'Tmax':0}
      self.template = ubertemp.template()
      # R_V as a function of which calibration fit number (see Folatelli et
      #  al. (2009) table 9
      self.Rv_host = {1:0, 2:3.10, 3:1.50, 4:1.46, 5:1.46, 6:1.01}
      self.dRv_host = {1:0, 2:0, 3:0.11, 4:0.10, 5:0.33, 6:0.31}
      self.M0 = {1:-19.07,2:-19.39,3:-19.15,4:-19.13,5:-19.16,6:-19.11}
      self.dM0 = {1:0.01, 2:0.02, 3:0.02, 4:0.01, 5:0.03, 6:0.02}
      self.b = {1:1.03, 2:0.98, 3:0.94, 4:1.05, 5:0.94, 6:1.08}
      self.db = {1:0.25, 2:0.41, 3:0.11, 4:0.11, 5:0.12, 6:0.11}
      # B-X pseudo-colors from Folatelli et al. (2009) table 3
      self.colors = {'u':-0.32, 'B':0.0, 'V':-0.02, 'g':0.05, 'r':-0.09, 'i':-0.63,
            'Y':-0.69, 'J':-0.65, 'H':-0.79, 'K':-0.61, 'J_K':-0.65, 'H_K':-0.79}
      self.dcolors = {'u':0.04, 'B':0, 'V':0.01, 'g':0.02, 'r':0.02, 'i':0.02,
            'Y':0.03, 'J':0.02, 'H':0.03, 'K':0.05, 'J_K':0.02, 'H_K':0.03}
      self.color_slopes = {'u':-0.47, 'B':0.0, 'V':0.12, 'g':0.05, 'r':0.29,
                           'i':0.39, 'Y':0.63, 'J':0.67, 'H':0.66, 'K':0.26,
                           'J_K':0.67, 'H_K':0.66}
      self.dcolor_slopes = {'u':0.25, 'B':0, 'V':0.05, 'g':0.06, 'r':0.07,
                            'i':0.08, 'Y':0.17, 'J':0.10, 'H':0.11, 'K':0.18,
                            'J_K':0.10, 'H_K':0.11}
      self.do_Robs = 0
      self.Robs = {}


   def setup(self):
      if 'EBVhost' not in self.args:
         if len(self._fbands) < 2:
            raise RuntimeError, "Error:  to solve for EBVhost, you need to fit more than one filter"

      self.calibration = self.args.get('calibration',6)
      self.gen = self.args.get('gen',2)
      for band in self._fbands:
         #cal = self.args.get('cal',6)
         cal = self.calibration
         self.Robs[band] = kcorr.R_obs(band, self.parent.z, 0, 0.01, 0,
               self.Rv_host[cal], self.parent.Rv_gal, self.parent.k_version,
               redlaw=self.parent.redlaw)
      
   def guess(self, param):
      s = self.parent
      if param == 'Tmax':
         Tmaxs = []
         for f in s.data:
            Tmaxs.append(s.data[f].MJD[argmin(s.data[f].mag)])
         return median(Tmaxs)

      if param == 'DM':
         # Quick DM based on Ho = 72
         if s.z < 1e-10:
            raise ValueError, "SN redshift is too close to zero.  Set it properly."
         return 43.11 + 5*log10(s.z)

      if param == 'dm15':
         # choose just the average dm15:
         return(1.1)

      return(0.0)

   def __call__(self, band, t, extrap=False):
      self.template.mktemplate(self.dm15)
      if len(shape(t)) == 0:
         t = array([t])
      t = t - self.Tmax
      rband = self.parent.restbands[band]

      # Now build the lc model
      temp,etemp,mask = self.template.eval(rband, t, self.parent.z, 
            gen=self.gen, extrap=extrap)
      K,mask2 = self.kcorr(band, t)
      temp = temp + K

      # Apply reddening correction:
      # Figure out the reddening law
      if self.do_Robs:
         self.Robs[band] = kcorr.R_obs(band, self.parent.z, t, self.EBVhost,
               self.parent.EBVgal, self.Rv_host[self.calibration], 
               self.parent.Rv_gal, self.parent.k_version, 
               redlaw=self.parent.redlaw)
         temp = temp + self.Robs[band]*(self.EBVhost + self.parent.EBVgal)
      else:
         # Apply Robs*EBVgal:
         R = self.MWR(band, t)
         temp = temp + self.Robs[band]*self.EBVhost + R*self.parent.EBVgal
      temp = temp + self.DM + self.MMax(rband, self.calibration)
      return temp,etemp,mask*mask2

   def get_max(self, bands, restframe=0, deredden=0):
      Tmaxs = []
      Mmaxs = []
      eMmaxs = []
      rbands = []
      self.template.mktemplate(self.dm15)
      for band in bands:
         rband = self.parent.restbands[band]
         # find where the template truly peaks:
         x0 = brent(lambda x: self.template.eval(rband, x, gen=self.gen)[0], brack=(0.,5.))
         Tmaxs.append(x0 + self.Tmax)
         mmax = self.DM + self.MMax(rband, self.calibration)
         if not restframe and band in self.parent.ks_tck:
            # add the K-correction
            mmax = mmax + scipy.interpolate.splev(Tmaxs[-1], self.parent.ks_tck[band])
         if not deredden:
            if self.do_Robs:
               Robs = kcorr.R_obs(band, self.parent.z, x0, self.EBVhost, 
                     self.parent.EBVgal, self.Rv_host[self.calibration], 
                     self.parent.Rv_gal, 'H3', redlaw=self.parent.redlaw)
               mmax = mmax + Robs*(self.EBVhost + self.parent.EBVgal)
            else:
               # Apply Robs*EBVgal:
               if band in self.parent.Robs:
                  if type(self.parent.Robs[band]) is type(()):
                     R = scipy.interpolate.splev(Tmaxs[-1], self.parent.Robs[band])
                  else:
                     R = self.parent.Robs[band]
               else:
                  EBV = max(self.parent.EBVgal, 0.01)
                  R = kcorr.R_obs(band, self.parent.z, 0, 0, EBV,
                                  self.Rv_host[self.calibration], 
                                  self.parent.Rv_gal,
                                  self.parent.k_version,
                                  redlaw=self.parent.redlaw)
               EBV = max(self.EBVhost, 0.01)
               Robs = kcorr.R_obs(band, self.parent.z, 0, EBV,
                     0, self.Rv_host[self.calibration], self.parent.Rv_gal,
                     self.parent.k_version, redlaw=self.parent.redlaw)
               mmax = mmax + Robs*self.EBVhost + R*self.parent.EBVgal
         Mmaxs.append(mmax)
         eMmaxs.append(self.errors['DM'])
         rbands.append(rband)
      return(Tmaxs, Mmaxs, eMmaxs, rbands)

   def MMax(self, band, calibration=6):
      '''Given self.dm15, return the absolute magnitude at maximum for the given
      filter [band].  The calibration paramter allows you to choose which
      fit (1-6) in Folatelli et al. (2009), table 9'''
      if band == 'Bs':
         return -19.319 + (self.dm15-1.1)*0.634
      elif band == 'Vs':
         return -19.246 + (self.dm15-1.1)*0.606
      elif band == 'Rs':
         return -19.248 + (self.dm15-1.1)*0.566
      elif band == 'Is':
         return -18.981 + (self.dm15-1.1)*0.524
      elif band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
         return self.M0[calibration] + (self.dm15-1.1)*self.b[calibration] -\
                self.colors[band] - self.color_slopes[band]*(self.dm15 -1.1)
      else:
         return -19.0

   def systematics(self, calibration=6, include_Ho=False):
      '''Returns the systematic errors in the paramters as a dictionary.  
      If no estimate is available, return None for that paramter.'''
      systs = dict.fromkeys(self.parameters.keys())
      # DM contains systematics for Ho, plus average of calibration
      #  uncertainties
      syst_DM = []
      syst_EBV = []
      weights = []
      for band in self._fbands:
         rb = self.parent.restbands[band]
         # make a call to get the weights
         mod,err,mask = self.__call__(band, self.parent.data[band].MJD)
         weights.append(sum(where(mask, power(err,-2), 0)))
         ddm15 = self.dm15 - 1.1
         Robs = kcorr.R_obs(band, self.parent.z, 0, self.EBVhost,
               self.parent.EBVgal, self.Rv_host[calibration], 
               self.parent.Rv_gal, self.parent.k_version, 
               redlaw=self.parent.redlaw)
         dRobs = Robs*self.dRv_host[calibration]/self.Rv_host[calibration]
         syst_DM.append(power(ddm15*self.db[calibration],2)+\
                        power(ddm15*self.dcolor_slopes[rb],2)+\
                        power(self.EBVhost*dRobs, 2)+\
                        power(self.dM0[calibration], 2)+\
                        power(self.dcolors[rb]*Robs, 2) +\
                        power(0.06*Robs, 2) +\
                        #power(2.17*velerr/(3e5*self.parent.z),2) +\
                        power(0.06,2))
         syst_EBV.append(power(0.06*Robs,2))
      syst_DM = array(syst_DM)
      weights = array(weights)
      systs['DM'] = sum(weights*syst_DM)/sum(weights)
      if include_Ho:
         systs['DM'] += power(2.17*0.1,2)     # assume 10% error in Ho
      systs['DM'] = sqrt(systs['DM'])
      systs['EBVhost'] = 0.06
      return(systs)

def read_table(file):
   f = open(file, 'r')
   lines = f.readlines()
   f.close()
   lines = map(string.strip, lines)
   lines = [line for line in lines if line[0] != "#"]
   lines = map(string.split, lines)
   a = {}; b = {};  c={};  Rv = {};  sig = {}
   ea = {}; eb = {};  ec={};  eRv = {}
   for line in lines:
      id = int(line[0])
      if id not in a:
         a[id] = {};  b[id] = {};  c[id] = {};  sig[id] = {}
         ea[id] = {};  eb[id] = {};  ec[id] = {}
      f = line[1]
      a[id][f] = float(line[2]); ea[id][f] = float(line[3])
      b[id][f] = float(line[4]); eb[id][f] = float(line[5])
      c[id][f] = float(line[6]); ec[id][f] = float(line[7])
      Rv[id] = float(line[8]); eRv[id] = float(line[9])
      sig[id][f] = float(line[10])
   return(a,ea,b,eb,c,ec,Rv,eRv,sig)



class EBV_model2(model):
   '''This model fits any number of lightcurves with CSP uBVgriYJHK templates
   or Prieto BsVsRsIs templates.  The parameters you can fit:

   - dm15 or st (decline rate or stretch)
   - Tmax (time of peak B maximum)
   - DM   (distance modulus)
   - EBVhost  (host galaxy extinction)

   The model is constructed by assuming a peak absolute magnitudes in tall
   the filters, as derived in Burns et al. 2011.  Calibrations were determined
   using MCMC modeling on all filters at once, determining M_X and b_X for
   each filter, and one value for R_V.  There are 2 calibrations, based on the 
   sample used to make the fit and the prior used on the extinction.  Default
   is 0, where the 2 red SNe are excluded and the blue sub-sample is used to
   anchor the colors.  Value of 1 is for the sample where the two red
   SNe were included.  A global extinction parameter EBVhost
   is applied to each light-curve, as well as Milky way extinction from 
   the SN object's EBVgal.  The value of R_V for the host galaxy is
   not a parameter, but is controled by the choice of calibration
   R_V for the galactic extinction is taken from the SN object (default 3.1).'''

   def __init__(self, parent, stype='st'):

      if stype not in ['dm15','st']:
         raise ValueError, "This model only supports the dm15 and st parameter"
      model.__init__(self, parent)
      self.rbs = ['u','B','V','g','r','i','Y','J','H']
      self.parameters = {'DM':None, stype:None, 'EBVhost':None, 'Tmax':None}
      self.errors = {'DM':0, stype:0, 'EBVhost':0, 'Tmax':0}
      if stype == 'dm15':
         self.template = ubertemp.template()
      else:
         self.template = ubertemp.stemplate()
      self.stype = stype
      
      if stype in ['st']:
         self.a,self.ea,self.b,self.eb,self.c,self.ec,self.Rv_host, self.eRv_host,self.sigSN = read_table(os.path.join(base,'st_calibration2.dat'))
      else:
         self.a,self.ea,self.b,self.eb,self.c,self.ec,self.Rv_host, self.eRv_host,self.sigSN = read_table(os.path.join(base,'dm15_calibration2.dat'))

      self.do_Robs = 0
      self.Robs = {}


   def setup(self):
      # check to see if we have more than one filter when solving for EBV
      if 'EBVhost' not in self.args:
         if len(self._fbands) < 2:
            raise RuntimeError, "Error:  to solve for EBVhost, you need to fit more than one filter"

      self.calibration = self.args.get('calibration',0)
      self.gen = 2

      for band in self._fbands:
         self.Robs[band] = kcorr.R_obs(band, self.parent.z, 0, 0.01, 0,
               self.Rv_host[self.calibration], self.parent.Rv_gal, 
               self.parent.k_version, redlaw=self.parent.redlaw)
      
   def guess(self, param):
      s = self.parent
      if param == 'Tmax':
         Tmaxs = []
         for f in s.data:
            Tmaxs.append(s.data[f].MJD[argmin(s.data[f].mag)])
         return median(Tmaxs)

      if param == 'DM':
         # Quick DM based on Ho = 72
         if s.z < 0.0001:
            return 0
         return 43.11 + 5*log10(s.z)

      if param == 'dm15':
         # choose just the average dm15:
         return(1.1)

      if param == 'st':
         return (1.0)

      return(0.0)

   def __call__(self, band, t, extrap=False):
      self.template.mktemplate(self.parameters[self.stype])
      t = t - self.Tmax
      rband = self.parent.restbands[band]

      # Now build the lc model
      temp,etemp,mask = self.template.eval(rband, t, self.parent.z, 
            gen=self.gen, extrap=extrap)
      K,mask2 = self.kcorr(band, t)
      temp = temp + K

      # Apply reddening correction:
      # Figure out the reddening law
      if self.do_Robs:
         self.Robs[band] = kcorr.R_obs(band, self.parent.z, t, self.EBVhost,
               self.parent.EBVgal, self.Rv_host[self.calibration], 
               self.parent.Rv_gal, self.parent.k_version,
               redlaw=self.parent.redlaw)
         temp = temp + self.Robs[band]*(self.EBVhost + self.parent.EBVgal)
      else:
         # Apply Robs*EBVgal:
         R = self.MWR(band, t)
         temp = temp + self.Robs[band]*self.EBVhost + R*self.parent.EBVgal
      temp = temp + self.DM + self.MMax(rband, self.calibration)
      return temp,etemp,mask*mask2

   def get_max(self, bands, restframe=0, deredden=0):
      Tmaxs = []
      Mmaxs = []
      eMmaxs = []
      rbands = []
      self.template.mktemplate(self.parameters[self.stype])
      for band in bands:
         rband = self.parent.restbands[band]
         # find where the template truly peaks:
         x0 = brent(lambda x: self.template.eval(rband, x, gen=self.gen)[0], brack=(0.,5.))
         Tmaxs.append(x0 + self.Tmax)
         mmax = self.DM + self.MMax(rband, self.calibration)
         if not restframe and band in self.parent.ks_tck:
            # add the K-correction
            mmax = mmax + scipy.interpolate.splev(Tmaxs[-1], self.parent.ks_tck[band])
         if not deredden:
            if self.do_Robs:
               Robs = kcorr.R_obs(band, self.parent.z, x0, self.EBVhost, 
                     self.parent.EBVgal, self.Rv_host[self.calibration], 
                     self.parent.Rv_gal, 'H3', redlaw=self.parent.redlaw)
               mmax = mmax + Robs*(self.EBVhost + self.parent.EBVgal)
            else:
               # Apply Robs*EBVgal:
               if band in self.parent.Robs:
                  if type(self.parent.Robs[band]) is type(()):
                     R = scipy.interpolate.splev(Tmaxs[-1], self.parent.Robs[band])
                  else:
                     R = self.parent.Robs[band]
               else:
                  EBV = max(self.parent.EBVgal, 0.01)
                  R = kcorr.R_obs(band, self.parent.z, 0, 0, EBV,
                                  self.Rv_host[self.calibration], 
                                  self.parent.Rv_gal, self.parent.k_version,
                                  redlaw=self.parent.redlaw)
               EBV = max(self.EBVhost, 0.01)
               Robs = kcorr.R_obs(band, self.parent.z, 0, EBV,
                     0, self.Rv_host[self.calibration], self.parent.Rv_gal,
                     self.parent.k_version, redlaw=self.parent.redlaw)
               mmax = mmax + Robs*self.EBVhost + R*self.parent.EBVgal
         Mmaxs.append(mmax)
         eMmaxs.append(self.errors['DM'])
         rbands.append(rband)
      return(Tmaxs, Mmaxs, eMmaxs, rbands)

   def MMax(self, band, calibration=1):
      '''Given self.dm15, return the absolute magnitude at maximum for the given
      filter [band].  The calibration paramter allows you to choose which
      fit (1-6) in Folatelli et al. (2009), table 9'''
      if band not in self.a[calibration]:
         raise ValueError, "Error, filter %s cannot be fit with this calibration"
      if self.stype in ['st']:
         delta = self.st - 1.0
      else:
         delta = self.dm15 - 1.1
      return self.a[calibration][band] + self.b[calibration][band]*delta +\
             self.c[calibration][band]*delta**2

   def systematics(self, calibration=1, include_Ho=False):
      '''Returns the systematic errors in the paramters as a dictionary.  
      If no estimate is available, return None for that paramter.'''
      systs = dict.fromkeys(self.parameters.keys())
      # DM contains systematics for Ho, plus average of calibration
      #  uncertainties
      syst_DM = []
      weights = []
      for band in self._fbands:
         rb = self.parent.restbands[band]
         # make a call to get the weights
         mod,err,mask = self.__call__(band, self.parent.data[band].MJD)
         weights.append(sum(where(mask, power(err,-2), 0)))
         if self.stype in ['st']:
            dst = self.st - 1.0
         else:
            dst = self.dm15 - 1.1
         Robs = kcorr.R_obs(band, self.parent.z, 0, self.EBVhost,
               self.parent.EBVgal, self.Rv_host[calibration], 
               self.parent.Rv_gal, self.parent.k_version,
               redlaw=self.parent.redlaw)
         dRobs = Robs*self.eRv_host[calibration]/self.Rv_host[calibration]
         syst_DM.append(power(dst*self.eb[calibration][rb],2)+\
                        power(self.EBVhost*dRobs, 2)+\
                        power(self.ea[calibration][rb], 2)+\
                        #power(2.17*velerr/(3e5*self.parent.z),2) +\
                        power(self.sigSN[calibration][rb],2))
      syst_DM = array(syst_DM)
      weights = array(weights)
      systs['DM'] = sum(weights*syst_DM)/sum(weights)
      if include_Ho:
         systs['DM'] += power(2.17*0.1,2)     # assume 10% error in Ho
      systs['DM'] = sqrt(systs['DM'])
      return(systs)


class max_model(model):
   '''This model is very similar to the EBV_model, but instead of having an
   extinction parameter (EBVhost) that controls all the colors, we simply
   fit a peak magnitude for each band independently.  There are therefore
   a variable number of parameters, based on the number of bands you fit:
   - dm15 (decline rate parameter)
   - Tmax (time of maximum in B)
   - Xmax (peak magnitude in restband X).  N of these for N bands.
   The lightcurves are constructed by offsetting each template by Xmax vertically
   and Tmax horizontally, as a function of dm15.  Each band is also offset by
   R_X*EBVgal, where EBVgal is taken from the parent SN object (as is the value
   of R_V).'''

   def __init__(self, parent, stype='dm15'):

      if stype not in ['dm15','st']:
         raise ValueError, "This model only supports dm15 and st as shape parameters"
      self.stype = stype
      model.__init__(self, parent)
      self.rbs = ['u','B','V','g','r','i','Y','J','H','K']
      self.M0s = {'u':-18.82, 'B':-19.11, 'V':-19.12, 'g':-19.16, 'r':-19.03,
                  'i':-18.50, 'Y':-18.45, 'J':-18.44, 'H':-18.38, 'K':-18.42,
                  'J_K':-18.44, 'H_K':-18.38,
                  'Bs':-19.319, 'Vs':-19.246, 'Rs':-19.248, 'Is':-18.981}
      self.parameters = {stype:None, 'Tmax':None}
      self.errors = {stype:0, 'Tmax':0}
      if stype == 'dm15':
         self.template = ubertemp.template()
         self.rbs = self.rbs + ['Bs','Vs','Rs','Is']
      else:
         self.template = ubertemp.stemplate()
      self.do_Robs = 1
      self.R_obs = {}

   def setup(self):
      '''Since we have a variable number of parameters, we need to 
      do this dynamically before the fitting is done.'''
      self.N_bands = len(self._fbands)
      self._rbs = [self.parent.restbands[band] for band in self._fbands]
      # start with the set we always have:
      shape = self.stype
      pars = {shape:self.parameters[shape], 'Tmax':self.parameters['Tmax']}
      errs = {shape:self.errors[shape], 'Tmax':self.errors['Tmax']}
      # now build up maxs, but use previously fit values if they exist.
      for band in self._rbs:
         if band+"max" not in pars:
            if band+"max" in self.parameters:
               pars[band+"max"] = self.parameters[band+"max"]
               errs[band+"max"] = self.errors[band+"max"]
            else:
               pars[band+"max"] = None
               errs[band+"max"] = 0
      self.parameters = pars
      self.errors = errs
      self.gen = self.args.get('gen',2)

   def guess(self, param):
      s = self.parent
      if param == 'Tmax':
         Tmaxs = []
         for f in s.data:
            Tmaxs.append(s.data[f].MJD[argmin(s.data[f].mag)])
         return median(Tmaxs)

      if param == 'st':
         return 1.0

      if param.find('max') > 0:
         filt = param.replace('max','')
         if filt in self.M0s:
            M0 = self.M0s[param.replace('max','')]
         else:
            M0 = -19
         if s.z < 1e-6:
            return M0
         return 43.11 + 5*log10(s.z) + M0
      
      if param == 'dm15':
         # choose just the average dm15:
         return(1.1)

      return(0.0)

   #def prior(self):
   #   if self.stype == 'dm15':
   #      return _quniform_prior(self.parameters['dm15'], 0.4, 2.5,0.01)
   #   elif self.stype == 'st':
   #      return _quniform_prior(self.parameters['st'], 0.2, 1.3, 0.01)


   def __call__(self, band, t, extrap=False):
      if debug:  print ">>>   Now in max_model"
      if debug:  print ">>>> setting shape parameter to ", self.parameters[self.stype]
      self.template.mktemplate(self.parameters[self.stype])
      t = t - self.Tmax
      rband = self.parent.restbands[band]

      # Now build the lc model
      if debug:  
         print ">>>> calling template.eval with"
         print "rband = %s, t =" % (rband), t

      temp,etemp,mask = self.template.eval(rband, t, self.parent.z, 
            gen=self.gen, extrap=extrap)
      if debug:  print ">>>>  done."
      K,mask2 = self.kcorr(band, t)
      temp = temp + K

      # Apply Max
      temp = temp + self.parameters[rband+'max'] 

      # Apply Robs*EBVgal:
      R = self.MWR(band, t)
      temp = temp + R*self.parent.EBVgal
      return temp,etemp,mask*mask2

   def get_max(self, bands, restframe=0, deredden=0):
      Tmaxs = []
      Mmaxs = []
      eMmaxs = []
      rbands = []
      self.template.mktemplate(self.parameters[self.stype])
      for band in bands:
         rband = self.parent.restbands[band]
         # find where the template truly peaks:
         x0 = brent(lambda x: self.template.eval(rband, x, gen=self.gen)[0], brack=(0.,5.))
         Tmaxs.append(x0 + self.Tmax)
         if rband+"max" not in self.parameters:
            raise ValueError, "Trying to find max of %s, but haven't solved for %s" %\
                  (band, rband+"max")
         mmax = self.parameters[rband+'max']
         if not restframe and band in self.parent.ks_tck:
            # add the K-correction
            mmax = mmax + scipy.interpolate.splev(Tmaxs[-1], self.parent.ks_tck[band])
         if not deredden:
            if band in self.parent.Robs:
               if type(self.parent.Robs[band]) is type(()):
                  R = scipy.interpolate.splev(x0+self.Tmax, self.parent.Robs[band])
               else:
                  R = self.parent.Robs[band]
            else:
               R = kcorr.R_obs(band, self.parent.z, int(floor(x0)), 0.0, 
                     self.parent.EBVgal, self.parent.Rv_gal, 
                     version=self.parent.k_version, redlaw=self.parent.redlaw)

            mmax = mmax + R*self.parent.EBVgal
         Mmaxs.append(mmax)
         eMmaxs.append(self.errors[rband+'max'])
         rbands.append(rband)
      return(Tmaxs, Mmaxs, eMmaxs, rbands)

   def _extra_error(self, parameter):
      '''Return the extra error we computed in the SNooPy paper when
      observations start after Tmax.'''
      days = array([0.,5.,10.,15.,20.])
      errors = {'dm15':array([0.00, 0.00, 0.03, 0.03, 0.03]),
                'Tmax':array([0.07, 0.16, 0.21, 0.21, 0.46]),
                'umax':array([0.00, 0.03, 0.05, 0.06, 0.06]),
                'Bmax':array([0.00, 0.02, 0.02, 0.03, 0.04]),
                'Vmax':array([0.00, 0.00, 0.01, 0.03, 0.03]),
                'gmax':array([0.00, 0.01, 0.01, 0.01, 0.02]),
                'rmax':array([0.00, 0.01, 0.02, 0.02, 0.02]),
                'imax':array([0.01, 0.03, 0.05, 0.05, 0.05]),
                'Ymax':array([0.01, 0.03, 0.03, 0.03, 0.04]),
                'Jmax':array([0.02, 0.04, 0.05, 0.05, 0.05]),
                'Hmax':array([0.01, 0.02, 0.02, 0.02, 0.02]),
                'Kmax':array([0.00, 0.00, 0.00, 0.00, 0.00])}
      errors['st'] = errors['dm15']*13.74/30    # conversion factor from dm15 to s
      id = parameter.find('max')
      if id > 0:
         f = parameter[0:id]
         for band in self.parent.data:
            if self.parent.restbands[band] == f:
               break
         t0 = (self.parent.data[band].MJD - self.Tmax).min()/(1+self.parent.z)
         if t0 < 0:  return 0
         id = argmax(absolute(t0-days))
         return errors[parameter][id]

      t0s = []
      for band in self.parent.data:
         rband = self.parent.restbands[band]
         if rband+'max' not in self.parameters:
            continue
         t0 = (self.parent.data[band].MJD - self.Tmax).min()/(1+self.parent.z)
         t0s.append(t0)
      t0 = median(t0s)
      if t0 < 0:
         return 0
      id = argmax(absolute(t0 - days))
      return errors[parameter][id]

class max_model2(model):
   '''Same as max_model, but here we let Tmax for each filter be a free parameter.'''

   def __init__(self, parent, stype = 'dm15'):

      if stype not in ['dm15','st']:
         raise ValueError, "This model only supports dm15 and st as shape parameters"
      self.stype = stype
      model.__init__(self, parent)
      self.rbs = ['u','B','V','g','r','i','Y','J','H','K',
            'H_K','J_K']
      if self.stype == 'dm15':
         self.rbs = self.rbs + ['Bs','Vs','Rs','Is']
      self.M0s = {'u':-18.82, 'B':-19.11, 'V':-19.12, 'g':-19.16, 'r':-19.03,
                  'i':-18.50, 'Y':-18.45, 'J':-18.44, 'H':-18.38, 'K':-18.42,
                  'J_K':-18.44, 'H_K':-18.38,
                  'Bs':-19.319, 'Vs':-19.246, 'Rs':-19.248, 'Is':-18.981,
                  'UVM2':-19, 'UVW1':-19, 'UVW2':-19}
      self.parameters = {stype:None}
      self.errors = {stype:0}
      if self.stype == 'dm15':
         self.template = ubertemp.template()
      else:
         self.template = ubertemp.stemplate()
      self.do_Robs = 1
      self.R_obs = {}

   def setup(self):
      '''Since we have a variable number of parameters, we need to 
      do this dynamically before the fitting is done.'''
      self.N_bands = len(self._fbands)
      self._rbs = [self.parent.restbands[band] for band in self._fbands]
      # start with the set we always have:
      pars = {self.stype:self.parameters[self.stype]}
      errs = {self.stype:self.errors[self.stype]}
      # now build up maxs, but use previously fit values if they exist.
      for band in self._rbs:
         if band+"max" not in pars:
            if band+"max" in self.parameters:
               pars[band+"max"] = self.parameters[band+"max"]
               errs[band+"max"] = self.errors[band+"max"]
            else:
               pars[band+"max"] = None
               errs[band+"max"] = 0
         if 'T'+band+"max" not in pars:
            if 'T'+band+"max" in self.parameters:
               pars['T'+band+"max"] = self.parameters['T'+band+"max"]
               errs['T'+band+"max"] = self.errors['T'+band+"max"]
            else:
               pars['T'+band+"max"] = None
               errs['T'+band+"max"] = 0
      self.parameters = pars
      self.errors = errs
      self.gen = self.args.get('gen',2)

   def guess(self, param):
      s = self.parent

      if param == 'st':
         return(1.)

      if param[0] == 'T':
         rfilt = param[1:].replace('max','')
         for f in self._fbands:
            if self.parent.restbands[f] == rfilt:  filt=f
         return s.data[filt].MJD[argmin(s.data[filt].mag)]

      elif param.find('max') > 0:
         M0 = self.M0s[param.replace('max','')]
         if s.z < 1e-6:
            return M0
         return 43.11 + 5*log10(s.z) + M0
      
      if param == 'dm15':
         # choose just the average dm15:
         return(1.1)

      return(0.0)

   def __call__(self, band, t, extrap=False):
      self.template.mktemplate(self.parameters[self.stype])
      rband = self.parent.restbands[band]
      Tmax = self.parameters['T'+rband+'max']
      t = t - Tmax

      # Now build the lc model
      temp,etemp,mask = self.template.eval(rband, t, self.parent.z, 
            gen=self.gen, extrap=extrap)
      # If k-corrections are there, use them
      K,mask2 = self.kcorr(band, t)
      temp = temp + K

      # Apply Max
      temp = temp + self.parameters[rband+'max'] 

      # Apply Robs*EBVgal:
      R = self.MWR(band, t)
      temp = temp + R*self.parent.EBVgal
      return temp,etemp,mask*mask2

   def get_max(self, bands, restframe=0, deredden=0):
      Tmaxs = []
      Mmaxs = []
      eMmaxs = []
      rbands = []
      self.template.mktemplate(self.parameters[self.stype])
      for band in bands:
         rband = self.parent.restbands[band]
         # find where the template truly peaks:
         x0 = brent(lambda x: self.template.eval(rband, x, gen=self.gen)[0], brack=(0.,5.))
         if rband+"max" not in self.parameters:
            raise ValueError, "Trying to find max of %s, but haven't solved for %s" %\
                  (band, rband+"max")
         Tmaxs.append(x0 + self.parameters['T'+rband+'max'])
         mmax = self.parameters[rband+'max']
         if not restframe and band in self.parent.ks_tck:
            # add the K-correction
            mmax = mmax + scipy.interpolate.splev(Tmaxs[-1], self.parent.ks_tck[band])
         if not deredden:
            if band in self.parent.Robs:
               if type(self.parent.Robs[band]) is type(()):
                  R = scipy.interpolate.splev(x0+self.parameters['T'+rband+'max'], 
                        self.parent.Robs[band])
               else:
                  R = self.parent.Robs[band]
            else:
               R = kcorr.R_obs(band, self.parent.z, int(floor(x0)), 0.0, 
                     self.parent.EBVgal, self.parent.Rv_gal, 
                     version=self.parent.k_version, redlaw=self.parent.redlaw)

            mmax = mmax + R*self.parent.EBVgal
         Mmaxs.append(mmax)
         eMmaxs.append(self.errors[rband+'max'])
         rbands.append(rband)
      return(Tmaxs, Mmaxs, eMmaxs, rbands)

   def __getattr__(self, attr):
      if 'parameters' in self.__dict__:
         if attr in self.__dict__['parameters']:
            return self.__dict__['parameters'][attr]
         if attr == 'Tmax' and 'TBmax' in self.__dict__['parameters']:
            return self.__dict__['parameters']['TBmax']
      if attr in self.__dict__:
         return self.__dict__[attr]
      else:
         raise AttributeError, "Attribute %s not found" % (attr)



class Rv_model(model):
   '''This model fits any number of lightcurves with CSP uBVgriYJHK templates
   or Prieto BsVsRsIs templates.  The parameters you can fit:

   - dm15 (decline rate)
   - Tmax (time of peak B maximum)
   - Bmax (maximum B magnitude)
   - EBVhost  (host galaxy extinction)
   - Rv (host galaxy reddening law)

   The model is constructed by assuming a peak observed magnitude in B, 
   value of dm15, and Tmax.  This will fit the B-lightcurve.  To fit
   the others, we use the Burns et al. (2011) calibrations to compute
   B - X  intrinsic colors, add extinction consistent with the current
   value of Rv and EBVhost to get predicted maximum magnitudes.'''

   def __init__(self, parent, stype='dm15'):

      model.__init__(self, parent)
      self.rbs = ['u','B','V','g','r','i','Y','J','H']
      if stype == 'dm15':
         self.rbs = self.rbs + ['Bs','Vs','Rs','Is']
      self.parameters = {'Bmax':None, 'dm15':None, 'EBVhost':None, 
            'Rv':None, 'Tmax':None}
      self.errors = {'Bmax':0, 'dm15':0, 'EBVhost':0, 'Tmax':0, 'Rv':0}
      self.template = ubertemp.template()
      # R_V as a function of which calibration fit number (see Folatelli et
      #  al. (2009) table 9
      self.M0 = {'u':[-18.64,-18.62],
                 'B':[-19.02,-19.02],
                 'V':[-19.01,-19.00],
                 'g':[-19.07,-19.07],
                 'r':[-18.93,-18.92],
                 'i':[-18.35,-18.32],
                 'Y':[-18.35,-18.33],
                 'J':[-18.44,-18.43],
                 'H':[-18.26,-18.25]}
      self.eM0 = {'u':[0.08,0.08],
                  'B':[0.06,0.06],
                  'V':[0.04,0.04],
                  'g':[0.06,0.05],
                  'r':[0.04,0.03],
                  'i':[0.03,0.02],
                  'Y':[0.02,0.02],
                  'J':[0.02,0.02],
                  'H':[0.02,0.02]}
      self.b = {'u':[0.58,0.46],
                'B':[0.32,0.32],
                'V':[0.33,0.32],
                'g':[0.31,0.33],
                'r':[0.26,0.26],
                'i':[0.14,0.11],
                'Y':[0.10,0.10],
                'J':[0.10,0.11],
                'H':[0.12,0.11]}
      self.eb = {'u':[0.32, 0.31],
                 'B':[0.26, 0.24],
                 'V':[0.18, 0.15],
                 'g':[0.24, 0.22],
                 'r':[0.15, 0.12],
                 'i':[0.11, 0.09],
                 'Y':[0.07, 0.07],
                 'J':[0.07, 0.07],
                 'H':[0.06, 0.07]}
      self.do_Robs = 0
      self.Robs = {}

   def setup(self):
      # check to see if we have more than one filter when solving for EBV
      if 'EBVhost' not in self.args:
         if len(self._fbands) < 2:
            raise RuntimeError, "Error:  to solve for EBVhost, you need to fit more than one filter"

      self.calibration = self.args.get('calibration',0)
      self.gen = self.args.get('gen',2)
      
   def guess(self, param):
      s = self.parent
      if param == 'Tmax':
         Tmaxs = []
         for f in s.data:
            Tmaxs.append(s.data[f].MJD[argmin(s.data[f].mag)])
         return median(Tmaxs)

      if param == 'Bmax':
         if s.z < 0.0001:
            return self.M0['B'][self.calibration]
         return 43.11 + 5*log10(s.z) + self.M0['B'][self.calibration]

      if param == 'dm15':
         # choose just the average dm15:
         return(1.1)

      if param == 'EBVhost':
         return(0.0)

      if param == 'Rv':
         return(2.0)

      return(0.0)

   def __call__(self, band, t, extrap=False):
      self.template.mktemplate(self.dm15)
      t = t - self.Tmax
      rband = self.parent.restbands[band]

      # Now build the lc model
      temp,etemp,mask = self.template.eval(rband, t, self.parent.z, 
            gen=self.gen, extrap=extrap)
      K,mask2 = self.kcorr(band, t)
      temp = temp + K

      # Apply reddening correction:
      # Figure out the reddening law
      # Milky-Way reddening:
      if band in self.parent.Robs:
         if type(self.parent.Robs[band]) is type(()):
            Rmw = scipy.interpolate.splev(t+self.Tmax, self.parent.Robs[band])
         else:
            Rmw = self.parent.Robs[band]
      # Host reddening:
      Ahost = kcorr.A_obs(rband, 0,  0, self.EBVhost, 0, self.Rv, \
            self.parent.Rv_gal, 'H3')
      temp = temp + Ahost + Rmw*self.parent.EBVgal
      temp = temp + self.Bmax + self.MMax(rband, self.calibration) - \
            self.MMax('B', self.calibration)
      return temp,etemp,mask*mask2

   def get_max(self, bands, restframe=0, deredden=0):
      Tmaxs = []
      Mmaxs = []
      eMmaxs = []
      rbands = []
      self.template.mktemplate(self.dm15)
      for band in bands:
         rband = self.parent.restbands[band]
         # find where the template truly peaks:
         x0 = brent(lambda x: self.template.eval(rband, x, gen=self.gen)[0], brack=(0.,5.))
         Tmaxs.append(x0 + self.Tmax)
         mmax = self.Bmax + self.MMax(rband, self.calibration) - \
               self.MMax('B', self.calibation)
         if not restframe and band in self.parent.ks_tck:
            # add the K-correction
            mmax = mmax + scipy.interpolate.splev(Tmaxs[-1], self.parent.ks_tck[band])
         if not deredden:
            if band in self.parent.Robs:
               if type(self.parent.Robs[band]) is type(()):
                  Rmw = scipy.interpolate.splev(x0+self.Tmax, 
                        self.parent.Robs[band])
               else:
                  Rmw = self.parent.Robs[band]
            Rhost = kcorr.R_obs(rband, 0, x0, self.EBVhost, 
                     self.parent.EBVgal, self.Rv_host, self.parent.Rv_gal, 'H3',
                     redlaw=self.parent.redlaw)
            mmax = mmax + Rhost*self.EBVhost + Rmw*self.parent.EBVgal
         Mmaxs.append(mmax)
         eMmaxs.append(self.errors['Bmax'])
         rbands.append(rband)
      return(Tmaxs, Mmaxs, eMmaxs, rbands)

   def MMax(self, band, calibration=1):
      '''Given self.dm15, return the absolute magnitude at maximum for the given
      filter [band].  The calibration paramter allows you to choose which
      fit (1-6) in Folatelli et al. (2009), table 9'''
      if band == 'Bs':
         return -19.319 + (self.dm15-1.1)*0.634
      elif band == 'Vs':
         return -19.246 + (self.dm15-1.1)*0.606
      elif band == 'Rs':
         return -19.248 + (self.dm15-1.1)*0.566
      elif band == 'Is':
         return -18.981 + (self.dm15-1.1)*0.524
      if band in ['u','B','V','g','r','i','Y','J','H']:
         return self.M0[band][calibration] + (self.dm15-1.1)*self.b[band][calibration]
      else:
         return -19.0

   def systematics(self, calibration=1, include_Ho=False):
      '''Returns the systematic errors in the paramters as a dictionary.  
      If no estimate is available, return None for that paramter.'''
      systs = dict.fromkeys(self.parameters.keys())
      # DM contains systematics for Ho, plus average of calibration
      #  uncertainties
      return(systs)


class color_model(model):
   '''This model fits any number of lightcurves with CSP uBVgriYJHK templates
   or Prieto BsVsRsIs templates.  The model is an observed B maximum (so B
   must be one of the restbands) and N-1 colors constructed from an intrinsic
   color-st/dm15 relation and extinction computed from E(B-V) and Rv. The
   parameters are:

   - st (decline rate)
   - Tmax (time of peak B maximum)
   - Bmax (maximum B magnitude)
   - EBVhost  (host galaxy extinction)
   - Rv (host galaxy reddening law)

   The assumed colors are take from Burns et al. (2014).
   '''

   def __init__(self, parent, stype='st', normfilter='B'):

      if stype == 'dm15':
         raise AttributeError, "This model only supports st parameter"
      model.__init__(self, parent)
      self.rbs = ['u','B','V','g','r','i','Y','J','H']
      self.normfilter = normfilter
      self.parameters = {'st':None, 'EBVhost':None, 
            'Rv':None, 'Tmax':None}
      self.parameters[self.normfilter+'max'] = None
      self.nparameters = {'a':None, 'b':None, 'c':None}
      self.enparameters = {'a':None, 'b':None, 'c':None}
      self.errors = {'st':0, 'EBVhost':0, 'Tmax':0, 'Rv':0}
      self.errors[self.normfilter+'max'] = 0
      self.template = ubertemp.stemplate()

   def __setstate__(self, d):
      # Needed to update older pickles. Cruft cruft cruft
      if 'normfilter' not in d:
         d['normfilter'] = 'B'
         d['ncolor'] = -1
      self.__dict__.update(d)

   def setup(self):
      # check to see if we have more than one filter when solving for EBV
      if 'EBVhost' not in self.args:
         if len(self._fbands) < 2:
            raise RuntimeError, "Error:  to solve for EBVhost, you need to fit more than one filter"

      # Make sure that at least one filter has restband for norm_filter
      nf = []
      for f in self.parent.data:
         if self.parent.restbands[f] == self.normfilter: nf.append(f)
      if len(nf) == 0:
         raise RuntimeError, "Error, to use color-model, you must have a B rest-frame observation"
      self.nf = nf

      cfile = os.path.join(base,"color_priors.pickle")
      if not os.path.isfile(cfile):
         raise ValueError, "Calibration file not found: %s" % \
               (self.calibration+".pickle")
      self.redlaw = self.args.get('redlaw', 'ccm')
      self.rvprior = self.args.get('rvprior', 'uniform')
      f = open(cfile)
      cdata = pickle.load(f)
      f.close()
      if self.redlaw not in cdata or \
            self.rvprior not in cdata[self.redlaw]:
         raise ValueError, "Intrinsic colors for %s prior and %s reddening law not found" %\
               (self.rvprior,self.redlaw)

      d = cdata[self.redlaw][self.rvprior]
      self.colors = d['colors']
      #self.ref_band = self.colors[0][0]
      self.Xfilters = [color[1] for color in self.colors]
      if self.normfilter == 'B':
         self.ncolor = -1
      else:
         self.ncolor = self.Xfilters.index(self.normfilter)
      self.Nf = len(self.Xfilters)
      self.nparameters['a'] = d['mean'][0:self.Nf]
      self.nparameters['b'] = d['mean'][self.Nf:self.Nf*2]
      self.nparameters['c'] = d['mean'][2*self.Nf:self.Nf*3]
      self.enparameters['a'] = sqrt(d['vars'])
      # Now we setup a covariance matrix for the data for each
      # band.  Right now, only including intrinsic color error
      self.covar = d['covar']
      self.gen = self.args.get('gen',2)
      if self.rvprior == 'bin':
         self.mu_i = d['mu_i']
         self.tau_i = d['tau_i']
         self.bins = d['bins']
      elif self.rvprior == 'mix':
         self.mu_i = d['mu_i']
         self.tau_i = d['tau_i']
         self.pi_i = d['pi_i']

   def get_covar(self, band, t):
      '''Return a covariance matrix to handle systematics'''
      cov_f = ones((t.shape[0],t.shape[0]))
      if band == "B":
         return cov_f*0
      #else:
      #   return cov_f*self.evar[self.Xfilters.index(band)]

   def guess(self, param):
      s = self.parent
      if param == 'Tmax':
         Tmaxs = []
         for f in s.data:
            Tmaxs.append(s.data[f].MJD[argmin(s.data[f].mag)])
         return median(Tmaxs)

      if param == self.normfilter+'max':
         return median([s.data[f].mag.min() for f in self.nf])

      if param == 'st':
         # choose just the average dm15:
         return(1.0)

      if param == 'EBVhost':
         return(0.0)

      if param == 'Rv':
         return(2.0)

      return(0.0)

   def prior(self, param, value):
      if param == "Rv":
         if value < 0.5:
            return -inf
         if self.rvprior == 'uniform':
            return 0.0
         elif self.rvprior == 'bin':
            id = searchsorted(self.bins, self.parameters['EBVhost'])
            return gconst - 0.5*power(value-self.mu_i[id],2)*self.tau_i[id] + \
                  0.5*log(self.tau_i[id])
         elif self.rvprior == 'mix':
            p = sum(pi_i*sqrt(tau_i/2/pi)*exp(-0.5*power(value-self.mu_i,2)*self.tau_i))
            return log(p)

      return 0.0


   def __call__(self, band, t, extrap=False):
      self.template.mktemplate(self.st)
      t = t - self.Tmax
      rband = self.parent.restbands[band]

      # Now build the lc model
      temp,etemp,mask = self.template.eval(rband, t, self.parent.z, 
            gen=self.gen, extrap=extrap)
      K,mask2 = self.kcorr(band, t)
      temp = temp + K

      # Apply reddening correction:
      # Figure out the reddening law
      # Milky-Way reddening:
      #mwRB = redlaw.R_lambda('B', 3.1, self.parent.EBVgal, redlaw=self.redlaw)

      temp = temp + self.parameters[self.normfilter+'max']
      a = self.nparameters['a']
      b = self.nparameters['b']
      c = self.nparameters['c']
      # Colors and reddening differential:
      if self.ncolor >= 0:
         temp = temp + a[self.ncolor] + b[self.ncolor]*(self.st - 1) + \
               c[self.ncolor]*(self.st - 1)**2
      if rband != 'B':
         cid = self.Xfilters.index(rband)
         temp = temp - a[cid] - b[cid]*(self.st-1) - c[cid]*(self.st -1)**2
      # Host reddening
      #RB = redlaw.R_lambda('B', self.Rv, self.EBVhost)
      RX = redlaw.R_lambda(rband, self.Rv, self.EBVhost, redlaw=self.redlaw)
      mwRX = redlaw.R_lambda(rband, 3.1, self.parent.EBVgal, redlaw=self.redlaw)
      #temp = temp - (RB - RX)*self.EBVhost + mwRX*self.parent.EBVgal
      temp = temp + RX*self.EBVhost +  mwRX*self.parent.EBVgal

      # added dispersion in intrinsic color
      
      return temp,etemp,mask*mask2

   def get_max(self, bands, restframe=0, deredden=0):
      Tmaxs = []
      Mmaxs = []
      eMmaxs = []
      rbands = []
      self.template.mktemplate(self.dm15)
      for band in bands:
         rband = self.parent.restbands[band]
         # find where the template truly peaks:
         x0 = brent(lambda x: self.template.eval(rband, x, gen=self.gen)[0], brack=(0.,5.))
         Tmaxs.append(x0 + self.Tmax)
         mmax = self.Bmax + self.MMax(rband, self.calibration) - \
               self.MMax('B', self.calibation)
         if not restframe and band in self.parent.ks_tck:
            # add the K-correction
            mmax = mmax + scipy.interpolate.splev(Tmaxs[-1], self.parent.ks_tck[band])
         if not deredden:
            if band in self.parent.Robs:
               if type(self.parent.Robs[band]) is type(()):
                  Rmw = scipy.interpolate.splev(x0+self.Tmax, 
                        self.parent.Robs[band])
               else:
                  Rmw = self.parent.Robs[band]
            Rhost = kcorr.R_obs(rband, 0, x0, self.EBVhost, 
                     self.parent.EBVgal, self.Rv_host, self.parent.Rv_gal, 'H3',
                     redlaw=self.parent.redlaw)
            mmax = mmax + Rhost*self.EBVhost + Rmw*self.parent.EBVgal
         Mmaxs.append(mmax)
         eMmaxs.append(self.errors['Bmax'])
         rbands.append(rband)
      return(Tmaxs, Mmaxs, eMmaxs, rbands)

   def systematics(self, calibration=1, include_Ho=False):
      '''Returns the systematic errors in the paramters as a dictionary.  
      If no estimate is available, return None for that parameter.'''
      systs = dict.fromkeys(self.parameters.keys())
      # DM contains systematics for Ho, plus average of calibration
      #  uncertainties
      return(systs)
