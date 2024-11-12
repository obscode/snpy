#!/usr/bin/env python
from __future__ import print_function, absolute_import
import six
import sys,string,os
import warnings
have_sql = 1
from . import sqlmod
if 'SQLSERVER' in os.environ:
   sqlmod.setSQL(os.environ['SQLSERVER'])
have_sql = sqlmod.have_sql

try:
   from . import snemcee
except ImportError:
   snemcee = None

try:
   import corner as triangle
except:
   triangle=None
   
import types
import time
from . import plot_sne_mpl as plotmod
from .lc import lc           # the light-curve class
from numpy import *       # Vectors
from . import ubertemp             # a template class that contains these two
from . import kcorr                # Code for generating k-corrections
from . import bolometric
from .utils import IRSA_dust_getval as dust_getval

from .utils import fit_poly  # polynomial fitter
import scipy                # Scientific python routines
linalg = scipy.linalg       # Several linear algebra routines
from scipy.interpolate import interp1d
from .utils import fit_spline # My Spline fitting routines
from .filters import fset    # filter definitions.
from .filters import standards # standard SEDs
from . import mangle_spectrum      # SN SED mangling routines
import pickle
from . import model
from .utils.fit1dcurve import list_types,regularize
from .version import __version__

# Some useful functions in other modules which the interactive user may want:
getSED = kcorr.get_SED
Robs = kcorr.R_obs
Ia_w,Ia_f = getSED(0, 'H3')
Vega = standards.Vega.VegaB
BD17 = standards.Smith.bd17

def check_file(o):
   '''Return true if o is an open file.'''
   if six.PY2:
      return type(o) is types.FileType
   else:
      from io import IOBase
      return isinstance(o, IOBase)


def myupdate(d1, d2):
   '''Update keys in d1 based on keys in d2, but only if they exist in d1'''
   for key in d1:
      if key in d2:
         d1[key] = d2[key]

class dict_def:
   '''A class that acts like a dictionary, but if you ask for a key
   that is not in the dict, it returns the key instead of raising
   an exception.'''
   def __init__(self, parent, dict={}):
      self.dict = dict
      self.parent = parent

   def __getitem__(self, key):
      if key in self.dict:
         return self.dict[key]
      else:
         return key

   def __setitem__(self, key, value):
      self.dict[key] = value

   def __delitem__(self, key):
      self.dict.__delitem__(key)

   def __contains__(self, key):
      return self.dict.__contains__(key)

   def __iter__(self):
      return self.dict.__iter__()

   def __str__(self):
      ret = ""
      for key in list(self.parent.data.keys()):
         ret += "%s -> %s, " % (key, self.__getitem__(key))
      return ret
   def __repr__(self):
      return self.__str__()

   def keys(self):
      return list(self.dict.keys())

   def values(self):
      vals = []
      for key in list(self.parent.data.keys()):
         val = self.__getitem__(key)
         if val not in vals: vals.append(val)
      return vals

class sn(object):
   '''This class is the heart of SNooPy.  Create a supernova object by 
   calling the constructor with the name of the superova as the argument.  
   e.g::
   
      In[1]:  s = sn('SN1999T')
   
   if the supernova is in the SQL database, its lightcurve data will be loaded 
   into its member data.  Once the object is created, use its member data 
   and functions to do your work.  Of course, you can have multiple 
   supernovae defined at the same time.

   Args:
      name (str): A SN name, or a filename containing data.
      source (sqlbase): An instance of sqlmod.sqlbase class (for database access)
      ra (float): degrees Right-ascention (J2000), of object
      dec (float): degrees Declination (J2000) of object
      z (float): heliocentric redshift of object
   '''

   def __init__(self, name, source=None, ra=None, dec=None, z=0,
         get_spectra=False):
      '''Create the object.  Only required parameter is the [name].  If this 
      is a new object, you can also specify [ra], [dec], and [z].'''
      self.__dict__['data'] = {}        # the photometric data, one for each band.
      self.__dict__['model'] = model.EBV_model(self)
      self.template_bands = ubertemp.template_bands

      self.Version = __version__    # A version-stamp for when we upgrade
      self.name = name
      self.z = z              # Redshift of the SN
      self.ra = ra              # Coordinates
      self.decl = dec
      self.filter_order = None  # The order in which to plot the filters
      self.xrange = None 
      self.yrange = None        # Impose any global plotting ranges?
      self.Rv_gal = 3.1         #  galaxy
      self.EBVgal = 0.0
      self.fit_mag = False     # fit in magnitude space?

      self.restbands = dict_def(self, {})   # the band to which we are fitting 
                                            #  for each band
      self.ks = {}          # k-corrections
      self.ks_mask = {}     # mask for k-corrections
      self.ks_tck = {}      # spline rep of k-corrections
      self.Robs = {}          # The observed R based on Rv and Ia spectrum

      self.p = None
      self.replot = 1          # Do we replot every time the fit finishes?
      self.quiet = 1           # Have copious output?

      if source is None:
         if ra is None or dec is None:
            if have_sql:
               self.sql = sqlmod.default_sql
               self.read_sql(self.name, spectra=get_spectra)
               self._sql_read_time = time.gmtime()
            else:
               print("Warning:  ra and/or decl not specified and no source specified.")
               print("   setting ra, decl and z = 0") 
               self.ra = 0;  self.decl = 0;
         else:
            self.ra = ra
            self.decl = dec
      else:
         self.sql = source
         self.read_sql(self.name, spectra=get_spectra)

      #self.summary()
      self.getEBVgal()
      self.get_restbands()     # based on z, assign rest-frame BVRI filters to 
                               # data
      self.k_version = 'H3'
      self.k_extrapolate = False   # Extrapolate SED beyond the ends?

   def __getattr__(self, name):
      if 'data' in self.__dict__:
         if name in self.data:
            return(self.data[name])
      if name == 'zcmb':
         return self.get_zcmb()

      if name == 'Tmax':
         if 'model' in self.__dict__:
            if 'Tmax' in self.__dict__['model'].parameters:
               if self.__dict__['model'].parameters['Tmax'] is not None:
                  return self.__dict__['model'].parameters['Tmax']
         for f in self.data:
            if self.restbands[f] == 'B':
               if self.data[f].Tmax is not None:
                  return self.data[f].Tmax
         return 0.0

      if name == 'dm15':
         if 'model' in self.__dict__:
            if 'dm15' in self.__dict__['model'].parameters:
               if self.__dict__['model'].parameters['dm15'] is not None:
                  return self.__dict__['model'].parameters['dm15']
         for f in self.data:
            if self.restbands[f] == 'B':
               if self.data[f].dm15 is not None:
                  return self.data[f].dm15
         return None
         
      if name == 'parameters':
         if 'model' in self.__dict__:
            return self.__dict__['model'].parameters
         else:
            raise AttributeError("Error, model not defined, so no paramters")
      if name == 'errors':
         if 'model' in self.__dict__:
            return self.__dict__['model'].errors
         else:
            raise AttributeError("Error, model not defined, so no errors")
      if name in self.parameters:
            return self.parameters[name]
      if name.replace('e_','') in self.errors:
         return self.errors[name.replace('e_','')]
      if name == 'dm15':
         if 'B' in self.data:
            return getattr(self.data['B'], 'dm15', None)
         else:
            return None
      if name == 'st':
         return None

      if name == 'redlaw':
         return 'ccm'
      raise AttributeError("Error:  attribute %s not defined" % (name))

   def __setattr__(self, name, value):
      if 'model' in self.__dict__:
         if name in self.__dict__['model'].parameters:
            self.__dict__['model'].parameters[name] = value
            return
      #if name == 'Rv_host'
      self.__dict__[name] = value

   def __eq__(self, other):
      '''Check if the data in this instance is equal to the data in another.
      Only checks raw data (light-curves).'''
      res = set(self.data.keys()) == set(other.data.keys())
      if not res:
         return False
      res = res and (abs(self.ra - other.ra) < 1e-7)
      res = res and (abs(self.decl - other.decl) < 1e-7)
      res = res and (abs(self.z - other.z) < 1e-9)
      for f in list(self.data.keys()):
         res = res and (self.data[f] == other.data[f])
      return res

   def __ne__(self, other):
      return not self.__eq__(other)

   def allbands(self):
      '''Returns all filters except those beginning with a "_"'''
      return [b for b in list(self.data.keys()) if b[0] != "_"]
      

   def choose_model(self, name, stype='st', set_restbands=True, **kwargs):
      '''A convenience function for selecting a model from the model module.
      [name] is the model to use.  The model will be used when self.fit() is 
      called and will contain all the parameters and errors. Refer to
      :class:`.model` for models and their parameters.

      Args:
         name (str): The name of the model (default: ``EBV_model``)
         stype (str): the template parameter (``dm15`` or ``st``)
         kwargs (dict): Any other arguments are sent to model's constructor

      Returns:
         None
      '''
      models = []
      for item in model.__dict__:
         obj = model.__dict__[item]
         try:
            if issubclass(obj, model.model):
               models.append(item)
         except:
            pass
      if name not in models:
         st = "Not a valid model.  Choose one of:  "+str(models)
         raise ValueError(st)

      self.model = model.__dict__[name](self, stype=stype, **kwargs)
      if not self.model.external_fitter:
         self.template_bands = [b for b in self.model.rbs \
            if b not in ['Bs','Vs','Rs','Is']]
      if set_restbands and not self.model.external_fitter:
         self.set_restbands()
     
   def get_mag_table(self, bands=None, dt=1.0, outfile=None):
      '''This routine returns a table of the photometry, where the data from
      different filters are grouped according intervals of dt.  The
      desired filters can be specified, otherwise all filters are
      returned. When data is missing, a value of 99.9 is inserted. Multiple
      observations are averaged in flux. The times are takent to be the
      average of all times in each bin.
      
      Args: 
         bands (list):  filters to include in the table
         dt (flaot): controls how to group by time: data are grouped
               into bins of width ``dt`` days.
         outfile (str or open file): optinal file name for output
      
      Returns:
         dict: numpy arrays keyed by:

         - ``MJD``: the epoch of observation
         - [band]: the magnitude in filter [band]
         - e_[band]: the error in [band]

      Raises:
         TypeError: the outfile is an incorrect type.
      '''

      if bands is None:  bands = list(self.data.keys())

      ret_data = {}
      # First, we make a list of observation dates from all the bands.
      times = [self.data[band].MJD for band in bands]
      times = sort(concatenate(times))

      # Create bin boundaries
      tbins = arange(times.min()+dt*0.5, times.max()+dt, dt) 
      bids = searchsorted(tbins, times)
      # Compute mean time in each bin. Empty bins will be assigned 'nan'
      with warnings.catch_warnings():
         warnings.filterwarnings('ignore', message=".*Mean of empty slice.*")
         warnings.filterwarnings('ignore', message=".*invalid value*")
         times = array([mean(times[bids == i]) for i in range(len(tbins))])
         gids = ~isnan(times)
      # Now we deal with edge cases: due to sparse sampling, you can get
      # bins separatedma/def 

      ret_data['MJD'] = times[gids]
      # Now loop through the bands and see where we need to fill in data
      for band in bands:
         bids = searchsorted(tbins, self.data[band].MJD)
         temp1 = 0.0*tbins + 99.9
         temp2 = 0.0*tbins + 9.9
         for i in range(len(tbins)):
            ggids = equal(bids, i)
            if sum(ggids) > 0:
               fs = self.data[band].flux[ggids]
               ws = power(self.data[band].e_flux[ggids],-2)
               af = sum(fs*ws)/sum(ws)
               ef = sqrt(1.0/sum(ws))
               temp1[i] = -2.5*log10(af) + self.data[band].filter.zp
               temp2[i] = 1.087*ef/af
         ret_data[band] = temp1[gids]
         ret_data["e_"+band] = temp2[gids]

      if outfile is not None:
         if type(outfile) is str:
            fp = open(outfile, 'w')
         elif check_file(outfile):
            fp = outfile
         else:
            raise TypeError("outfile must be a file name or file handle")
         JDlen = len(str(int(mean(ret_data['MJD'])))) + 3
         title = "MJD" + " "*(JDlen+2)
         for b in bands:  title += "%5s +/-   " % b
         print(title, file=fp)
         format = ("%%%d.2f  " % (JDlen-3)) + "%5.2f %4.2f  "*len(bands)
         for i in range(len(ret_data['MJD'])):
            data = [ret_data['MJD'][i]]
            for b in bands:  data += [ret_data[b][i], ret_data['e_'+b][i]]
            print(format % tuple(data), file=fp)
         fp.close()
         return
      else:
         return(ret_data)

   def color_stretch(self, Bband, Vband, interpolate=True, plot=False, 
         Nboot=10, method=None, **args):
      '''Given two filters, compute a [Bband]-[Vband] color curve, 
      fit with interpolator, and measure the time of maximum (reddest) color.
      If Tmax is defined for this SN, comptute the color-stretch (sBV).
      
      Args:
         Bband (str):  the observed filter corresponding to B-band
         Vband (str):  the observed filter corresponding to V-band
         interpolate (bool):  If true, and a model or interpolator exists,
                              interpolate missing data to compute colors
         Nboot (int):  Number of bootstrap iterations to get an error.
                      default: 10
         plot (bool): Make a plot of the fit
         method (str): interpolation method to use (default: GP or spline)
         args (dict):  arguments to be sent to fit1dcurve.

      Returns:
         2-tuple:  sBV, e_sBV

         sBV:  color stretch ((time of B-V maximum) - Tmax)/30
         e_sBV:  uncertainty (through bootstrapping)

         Values set to None if no maximum is found
      '''
      from .utils.fit1dcurve import Interpolator,functions

      if method is None:
         if 'gp' in functions:
            method='gp'
         else:
            method='spline2'

      if Bband not in self.data:
         raise ValueError("Bband filter {} not found".format(Bband))
      if Vband not in self.data:
         raise ValueError("Vband filter {} not found".format(Vband))

      # See if we can figure out a Tmax
      if self.data[Bband].Tmax is not None and self.data[Bband].Tmax > 0:
         Tmax = self.data[Bband].Tmax
      elif getattr(self, 'Tmax', None) is not None:
         Tmax = self.Tmax
      else:
         print("Warning:  no Tmax found")
         Tmax = 0

      t,BV,eBV,flags = self.get_color(Bband, Vband, interp=interpolate)
      gids = less(flags, 2) 

      interp = Interpolator(method, t[gids], BV[gids], eBV[gids], **args)

      tmaxs,maxs,dirs = interp.find_extrema()

      if plot:
         ts = linspace(t[gids].min(), t[gids].max(), 100)
         ys = interp(ts)[0]

      nids = less(dirs, 0)   # concave-down (i.e., maxima)
      if not any(nids):
         return None,None

      idx = argmax(maxs[nids])  # take the true maximum

      tBVmax = tmaxs[nids][idx]

      # Now bootstrap
      tBVmaxs = []
      if plot: yys = []
      for i in range(Nboot):
         interp.draw()
         btmaxs,bmaxs,bdirs = interp.find_extrema()
         if plot: yys.append(interp(ts)[0])
         if len(bdirs) == len(dirs) and all(bdirs == dirs):
            # same maxima
            tBVmaxs.append(btmaxs[nids][idx])
         else:
            bnids = less(bdirs, 0)
            if any(bnids):
               bidx = argmax(bmaxs[bnids])
               if absolute(btmaxs[bnids][bidx] - tBVmax) < 10.0:
                  tBVmaxs.append(btmaxs[bnids][bidx])
      # Do the stats
      if len(tBVmaxs) < Nboot:
         print("Warning, {}% of iterations failed".format(
            (Nboot - len(tBVmaxs))/Nboot*100))

      stBVmax = std(tBVmaxs)
      sBV = (tBVmax - Tmax)/(1+self.z)/30.0
      e_sBV = stBVmax/(1+self.z)/30.0
      if plot:
         eys = std(array(yys), axis=0)
         plotmod.plot_sBV(self, t[gids], BV[gids], eBV[gids], tBVmax, stBVmax, 
               ts, ys, eys, Tmax)

      return(sBV,e_sBV)


   def lira(self, Bband, Vband, interpolate=0, tmin=30, tmax=90, plot=0,
         dokcorr=True, kcorr=None, B14=False, deredden=True):
      '''Use the Lira Law to derive a color excess.  [Bband] and [Vband] 
      should be whichever observed bands corresponds to restframe B and V,
      respectively.  The color excess is estimated to be the median of the 
      offset between the Lira
      line and the data.  The uncertainty is 1.49 times the median absolute
      deviation of the offset data from the Lira line.
      
      Args:
         Bband (str): the observed filter corresponding to B-band
         Vband (str): the observed filter corresponding to V-band
         interpolate (bool): If true and a model (or interpolator) exists for
                             the observed filters, use it to interpolate
                             missing B or V data
         tmin/tmax (flaot): range over which to fit Lira Law
         plot (bool):  If True, produce a plot with the fit.
         dokcoor (bool):  If True, k-correct the data before fitting
         B14 (bool):  If True, use Burns et al. (2014) rather than Lira (1996)
         deredden (bool): If True, remove Milky-Way reddening before computing
                          Lira reddening and plotting
         
      Returns:
         4-tuple:  (EBV, error, slope, eslope)

         EBV: the E(B-V) color-excess
         error: undertainty based on fit
         slope: the late-time slope of the B-V color curve
         eslope: the error in the late-time slope
         
      '''
      if kcorr is not None:
         warnings.warn("Use of kcorr argument is deprecated. Use dokcorr "
                       " instead", stacklevel=2)
         dokcorr=kcorr

      # find V-maximum
      t_maxes,maxes,e_maxes,restbands = self.get_rest_max([Vband])
      Tmax = t_maxes[0]

      t,BV,eBV,flag = self.get_color(Bband, Vband, dokcorr=dokcorr, 
            interp=interpolate)

      if deredden:
         A_B = fset[Bband].R()*self.EBVgal
         A_V = fset[Vband].R()*self.EBVgal
         BV = BV - (A_B - A_V)

      # find all points that have data in both bands
      gids = equal(flag, 0)
      # If we're allowed to interpolate, add flag=1
      if interpolate:
         gids = gids + equal(flag, 1)

      # Now apply a time criterion
      gids = gids*greater_equal(t-Tmax, tmin)*less_equal(t-Tmax, tmax)

      # Now check that we actually HAVE some data left
      if not any(gids):
         raise RuntimeError("Sorry, no data available between t=%f and t=%f" % (tmin,tmax)) 
      
      # extract the data we want and convert to Vmax epochs
      t2 = compress(gids, (t-Tmax)/(1+self.z))
      BV2 = compress(gids, BV)
      eBV2 = compress(gids, eBV)
      
      # Next, solve for a linear fit (as diagnostic)
      w = power(eBV2,-2)
      c,ec = fit_poly.fitpoly(t2, BV2, w=w, k=1, x0=55.0)
      rchisq = sum(power(BV2 - c[0] - c[1]*(t2-55.),2)*w)/(len(BV2) - 2)
      ec = ec*sqrt(rchisq)

      if B14:
         st = getattr(self, 'st', 1.0)
         lira_BV = 0.78 - (0.0097 - 0.004*(st-1))*(t2 - 45.)
      else:
         lira_BV = 0.732 - 0.0095*(t2 - 55.0)
      w = power(eBV2,-2)
      lira_EBV = sum((BV2 - lira_BV)*w)/sum(w)
      e_lira_EBV = power(sum(w), -0.5)
      if B14:
         e_lira_EBV = sqrt(power(e_lira_EBV,2) + 0.04**2)

      print("Vmax occurred at %f" % (t_maxes[0]))
      print("Slope of (B-V) vs. t-Tvmax was %f(%f)" % (c[1], ec[1]))
      print("median E(B-V) = %f    1.49*mad(E(B-V)) = %f" % (lira_EBV, e_lira_EBV))
      if absolute(c[1] + 0.0118) > 3*ec[1]:
         print("WARNING:  fit slope differs from Lira Law by more than three sigma")

      if plot:
         plotmod.plot_lira(t, t2, t_maxes, BV, eBV, BV2, tmin, tmax, c)

      return (lira_EBV, e_lira_EBV, c[1], ec[1])

   def get_rest_max(self, bands, deredden=0):
      return self.get_max(bands, deredden=deredden, restframe=1)

   def get_max(self, bands, restframe=0, deredden=0, use_model=0):
      '''Get the  maximum magnitude in [bands] based on the currently
      defined model or spline fits.
      
      Args:
         bands (list): List of filters to fine maximum
         restframe (bool): If True, apply k-corrections (default: False)
         deredden (bool):  If True, de-redden using Milky-Way color excess
                           and any model color excess, if defined 
                           (default: False)
         use_model (bool): If True and both a model and interpolator are 
                           defined for a filter, use the model. 
                           (default: False, i.e., use interpolator)

      Returns:
         4-tuple:  (Tmax, Mmax, e_Mmax, rband)
            Tmax:  array of time-of-maximum, one for each of [bands]
            Mamx:  array of maximum magnitudes
            e_Mmax: error
            rband: list of rest-bands (if model was used to interpolate)

      '''
      if type(bands) is str:
         bands = [bands]
         scalar = True
      else:
         scalar = False
      model_bands = [b for b in bands if b in self.model._fbands]
      lc_model_bands = [b for b in bands if self.data[b].Mmax is not None]
      if use_model:
         lc_model_bands = [b for b in lc_model_bands if b not in model_bands]
      for band in bands:
         if band not in model_bands and band not in lc_model_bands:
            raise ValueError("Error:  filter %s has not been fit " % band + \
                  "with a light-curve yet, so I cannot compute it's maximum")
      N = len(bands)
      result = (zeros(N, dtype=float32), zeros(N, dtype=float32),
            zeros(N, dtype=float32), [""]*N)
      if len(model_bands) > 0:
         mod_result = self.model.get_max(model_bands, restframe=restframe,
               deredden=deredden)
      for i in range(N):
         b = bands[i]
         if b in model_bands:
            mid = model_bands.index(b)
            for j in range(4): result[j][i] = mod_result[j][mid]
         if b in lc_model_bands:
            if b in model_bands:
               print("Warning:  both model and spline fits present, using " +\
                     "spline values")
            if restframe:
               print("Warning:  can't k-correct spline fits, you're getting " +\
                     "observed maxima!")
            if deredden:
               if band in self.Robs:
                  if type(self.Robs[band]) is type(()):
                     R = scipy.interpolate.splev(self.data[b].Tmax, self.Robs[band])
                  else:
                     R = self.parent.Robs[band]
               else:
                  R = fset[band].R(wave=Ia_w, flux=Ia_f)
            else:
               R = 0

            result[0][i] = self.data[b].Tmax
            result[1][i] = self.data[b].Mmax - R*self.EBVgal
            result[2][i] = self.data[b].e_Mmax
            result[3][i] = b
      if scalar:
         return (result[0][0],result[1][0],result[2][0],result[3][0])
      return result

   def scorr(self, bands=None, SED='H3', merge=False):
      '''Compute the S-corrections for the named filters. The underlying
      SED is taken as the same that was used to compute k-corretions. If
      k-corrections have not been computed, we simply used the SED
      specified with no mangling.

      Args:
         bands (list or None): List of filters to S-correct or all if None
                               (default: None)
         SED (string):  If k-corrections have not been computed, we fall
                        back on using the SED spedified here. See 
                        kcorr.get_SED().
         merge (bool):  If True, then after computing the s-corrections,
                        merge all photometry with the same rest-band into
                        the same LC instances.

         Returns:
            None
         
         Effects:
            Upon successful completion, the following member variables will
            be populated:

            * self.Ss:       dictionary of S-corrections indexed by filter
            * self.Ss_mask:  dictionary indicating valid S-corrections
            * if merge, self.data is updated to merge filters to a common set
              and keep the old ones with a '_' prefix.
      '''
      if bands is None:
         bands = self.allbands()

      self.Ss = {}
      self.Ss_mask = {}
      for band in bands:
         if self.restbands[band] == band:
            # No S-correction
            self.Ss[band] = self.data[band].mag*0
            self.Ss_mask[band] = less(self.Ss[band],1)
            continue
         st = getattr(self, 'ks_s', 1.0)
         mopt = getattr(self, 'ks_mopts', {})
         if band not in mopt:
            # No k-corrections, we'll just use SED
            x = self.data[band].MJD
            # days since Bmax in the frame of the SN
            if SED != 'H3+L':
                days = (x - self.Tmax)/(1+self.z)/st
            days = days.tolist()
            self.Ss[band],self.Ss_mask[band] = list(map(array,kcorr.kcorr(days, 
               self.restbands[band], band, self.z, st, self.EBVgal, 0.0,
               version=SED, Scorr=True)))
            self.Ss_mask[band] = self.Ss_mask[band].astype(bool)
         else:
            self.Ss[band] = []
            self.Ss_mask[band] = []
            for i in range(len(self.data[band].MJD)):
               wave,flux,dum1,dum2 = self.get_mangled_SED(band,i)
               if wave is None:
                  S,f = 0,0
               else:
                  S,f = kcorr.S(wave,flux,fset[band],fset[self.restbands[band]],
                        z=self.z)
               self.Ss[band].append(S)
               self.Ss_mask[band].append(f)
            self.Ss[band] = array(self.Ss[band])
            self.Ss_mask[band] = array(self.Ss_mask[band])
      if merge:
         self.merge_filters(bands, method='scorr')

   def merge_filters(self, bands=None, method='scorr'):
      '''For all filters defined in self.data, merge those with the same
      restband (as defined by self.restbands) into a single filter
      instance.

      Args:
         bands (list)   filters to consider merging. If None, all filters
                        defined in self.data.
         method (str):  If 'scorr' (default), s-correct before merging.
                        If 'kcorr', do a cross-band k-correction before merging.
                        Any other string:  no correction is made.
         Returns:
            None
         
         Effects:
            Upon successful completion, all filters in bands will be merged
            into filter instances defined by self.restbands. Either K- or 
            S-corrections will be applied, so self.Ss and self.ks will be
            nullified.'''
      if bands is None:
         bands = self.allbands()
      for b in bands:
         if method == 'scorr' and (getattr(self, 'Ss', None) is None or\
               b not in self.Ss):
            raise ValueError("Error:  s-corrections for {} not found".format(b))
         elif method == 'kcorr' and b not in self.ks:
            raise ValueError("Error:  k-corrections for {} not found".format(b))

      # merge all filters with the same restband
      rbs = list(self.restbands.values())
      for rb in rbs:
         # list of bands to merge
         bs = [b for b in bands if self.restbands[b] == rb]
         if len(bs) == 0:
            continue
         lcs = [self.data[b] for b in bs]
         MJDs = concatenate([lc.MJD for lc in lcs])
         if method=='kcorr':
            mags = concatenate([self.data[b].magnitude - self.ks[b] \
                  for b in bs])
            masks = concatenate([self.data[b].mask & self.ks_mask[b] \
                  for b in bs])
         elif method=='scorr':
            mags = concatenate([self.data[b].magnitude - self.Ss[b] \
                  for b in bs])
            masks = concatenate([self.data[b].mask & self.Ss_mask[b] \
                  for b in bs])
         else:
            mags = concatenate([self.data[b].magnitude for b in bs])
            masks = ~isnan(mags)

         emags = concatenate([lc.e_mag for lc in lcs])
         SNRs = concatenate([lc.SNR for lc in lcs])
         sids = concatenate([lc.sids for lc in lcs])

         # Book-keeping
         for b in bs:
            if method == 'kcorr':
               del self.ks[b]
               del self.ks_mask[b]
               if b in self.ks_tck: del self.ks_tck[b]
            elif method == 'scorr':
               del self.Ss[b]
               del self.Ss_mask[b]
            del self.data[b]
            if self.filter_order is not None and  b in self.filter_order \
                  and b != rb:
               del self.filter_order[self.filter_order.index(b)]
               if len(self.filter_order) == 0:
                  self.filter_order = None
            if b in self.restbands: del self.restbands[b]

         #del self.ks[rb]
         #del self.ks_mask[rb]
         #del self.ks_tck[rb]
         #if getattr(self, 'Ss', None) is not None and rb in self.Ss:
         #   del self.Ss[rb]
         #   del self.Ss_mask[rb]

         # Create the LC
         self.data[rb] = lc(self, rb, MJDs, mags, emags, 
               SNR=SNRs, sids=sids)
         self.data[rb].mask[:] = masks[:]
         self.data[rb].time_sort()
         
    
   def kcorr(self, bands=None, mbands=None, mangle=1, interp=1, use_model=0, 
         min_filter_sep=400, use_stretch=1, **mopts):
      '''Compute the k-corrections for the named filters.
      In order to get the best k-corrections possible,
      we warp the SNIa SED (defined by self.k_version) to match the observed
      photometry. Not all bands will be observed on the same day (or some
      data may be less than reliable), so there are several arguments that
      control how the warping is done.

      Args:
         bands (list or None): List of filters to k-correct or all if None
                               (default: None)
         mbands (list of None): List of filters to use for mangling the SED.
                               (default: None:  same as bands)
         mangle (bool): If True, mangle (color-match) the SED to observed
                        colors. (default: True)
         interp (bool): If True, interpolate missing colors. (default: True)
         use_model (bool):  If True, use a model to interpolate colors, if
                            one exists (default: False)
         min_filter_sep (float): Filters whose effective wavelength are closer
                                 than this are rejected. (Default: 400 A)
         use_stretch (bool): If True, stretch the SED in time to match the
                             stretch/dm15 of the object. (Default: True)
         mopts (dict): Any additional arguments are sent to the function
                       mangle_spectrum.mangle_spectrum2()

      Returns:
         None

      Effects:
         Upon successful completion, the following member variables will be
         populated:

         * self.ks:      dictionary (indexed by filter) of k-corrections
         * self.ks_mask: dictionary indicating valid k-corrections
         * self.ks_tck: dictionary of spline coefficients for the k-corrections
                         (useful for interpolating the k-corrections).
         * self.mopts:  If mangling was used, contains the parameters of the 
                        mangling function.
      '''
      if use_stretch and self.k_version != '91bg':
         dm15 = getattr(self, 'dm15', None)
         st = getattr(self, 'st', None)
         if dm15 is None and st is None:
            raise AttributeError("Before you can k-correct with stretch, you"+\
                  " need to solve for dm15 or st, using either a model or LC fit")
         if dm15 is None:
            s = st
         else:
            if dm15 > 1.7:
               print("Warning:  dm15 > 1.7.  Using this stretch on the Hsiao SED")
               print("  is not recommended.  I'm setting stretch to 1.0.  You")
               print("  might consider using the 91bg SED template.")
               s = kcorr.dm152s(1.7)
            elif dm15 < 0.7:
               s = kcorr.dm152s(0.7)
            else:
               s = kcorr.dm152s(dm15)
      elif use_stretch and self.k_version == '91bg':
         print("Warning:  you asked for stretching the template SED, but")
         print("you have selected the 91bg template.  Setting stretch to 1.0.")
         s = 1.0
      else:
         s = 1.0
      self.ks_s = s

      if bands is None:  bands = list(self.data.keys())
      if mbands is None:  mbands = [b for b in bands]
      # Check the simple case:
      if not mangle:
         for band in bands:
            x = self.data[band].MJD
            # days since Bmax in the frame of the SN
            if self.k_version != 'H3+L':
                days = (x - self.Tmax)/(1+self.z)/s
            else:
                days = (x - self.Tmax)/(1+self.z)
            days = days.tolist()
            kextrap = getattr(self, 'k_extrapolate', False)
            self.ks[band],self.ks_mask[band] = list(map(array,kcorr.kcorr(days, 
               self.restbands[band], band, self.z, self.ks_s, self.EBVgal, 0.0,
               version=self.k_version, extrapolate=kextrap)))
            self.ks_mask[band] = self.ks_mask[band].astype(bool)
            #self.ks_tck[band] = scipy.interpolate.splrep(x, self.ks[band], k=1, s=0)
            if len(x) > 1:
               self.ks_tck[band] = fit_spline.make_spline(x, self.ks[band], x*0+1,
                                k=1, s=0, task=0, tmin=x.min(), anchor_dist=[0,0],
                                tmax=x.max())[0]
         return

      # Now see if we need to eliminate filters
      eff_waves = array([fset[band].eff_wave(Ia_w,Ia_f) for band in mbands])
      sids = argsort(eff_waves)
      eff_waves = eff_waves[sids]
      mbands = [mbands[sids[i]] for i in range(len(sids))]
      dwaves = eff_waves[1:] - eff_waves[0:-1]
      while any(less(dwaves, min_filter_sep)):
         bids = less(dwaves, min_filter_sep)
         mbands = [mbands[i] for i in range(len(bids)) if not bids[i]] + mbands[-1:]
         eff_waves = array([fset[band].eff_wave(Ia_w,Ia_f) for band in mbands])
         dwaves = eff_waves[1:] - eff_waves[0:-1]
      if not self.quiet:  print("Mangling based on filters:", mbands)

      restbands = [self.restbands[band] for band in bands]
      # now get the interpolated magnitudes all along the extent of the
      #  lightcurves.
      mags = []
      masks = []
      res = self.get_mag_table(bands)
      for band in mbands:
         bids = greater(res[band], 90)
         # find where we need to interpolate:
         if use_model:
            ev,eev,ma = self.model(band, res['MJD'])
            mags.append(ev)
            masks.append(ma)
         elif interp:
            ev,ma = self.data[band].eval(res['MJD'], t_tol=-1)
            #mags.append(where(bids, ev, res[band]))
            #masks.append(where(bids, ma, 1))
            mags.append(ev)
            masks.append(ma)
         else:
            mags.append(where(bids, 0.0, res[band]))
            masks.append(1-bids)

      mags = transpose(array(mags))
      masks = transpose(array(masks))
      # don't forget to convert to rest-frame epochs!
      if self.Tmax is None:
         raise AttributeError("Error.  self.Tmax must be set in oder to compute K-correctsions")
      t = res['MJD'] - self.Tmax
      if not any(greater_equal(t, -19)*less(t, 70)):
         raise RuntimeError("Error:  your epochs are all outside -20 < t < 70.  Check self.Tmax")
      kextrap = getattr(self, 'k_extrapolate', False)
      if  self.k_version != 'H3+L':
          use_t = t/(1+self.z)/s
      else:
          use_t = t/(1+self.z)
      kcorrs,mask,Rts,m_opts = kcorr.kcorr_mangle(use_t, bands, 
            mags, masks, restbands, self.z, sBV=s,
            colorfilts=mbands, version=self.k_version, full_output=1, 
            extrapolate=kextrap, **mopts)
      mask = greater(mask, 0)
      kcorrs = array(kcorrs)
      Rts = array(Rts)
      
      # At this point, we have k-corrections for all dates in res['MDJ']:
      #   kcorrs[i,j]  is kcorr for bands[j] on date res['MJD'][i]
      #   But there may be two observations separated by less than a day,
      #   in which case, they share the same k-correction.  So figure that out
      self.ks_mopts = {}
      for i in range(len(bands)):
         b = bands[i]
         self.ks_tck[b] = fit_spline.make_spline(res['MJD'], kcorrs[:,i],
                          res['MJD']*0+1, k=1, s=0, task=0, 
                          tmin=res['MJD'].min(), tmax = res['MJD'].max())[0]
         self.ks[b] = scipy.interpolate.splev(self.data[b].MJD, self.ks_tck[b])
         self.ks_mask[b] = array([mask[argmin(absolute(res['MJD'] - self.data[b].MJD[j])),i] \
               for j in range(len(self.data[b].MJD))]).astype(bool)
         self.ks_mopts[b] = [m_opts[argmin(absolute(res['MJD'] - self.data[b].MJD[j]))] \
               for j in range(len(self.data[b].MJD))]

         self.Robs[b] = fit_spline.make_spline(res['MJD'], Rts[:,i],
                          res['MJD']*0+1, k=1, s=0, task=0, 
                          tmin=res['MJD'].min(), tmax=res['MJD'].max())[0]


   def get_mangled_SED(self, band, i, normalize=True):
      '''After the mangle_kcorr function has been run, you can use this
      function to retrieve the mangled SED that was used to compute the
      k-correction for the [i]'th day in [band]'s light-curve.
      
      Args:
         band (str):  the refernece filter
         i (int):  index of filter [band]'s photometry
         normalize (bool): If True, normalize the SED to observed flux.

      Returns:
         4-tuple: (wave, mflux, oflux, mfunc)

         * wave: wavelength of SED in \AA
         * mflux: mangled SED flux
         * oflux: original (un-mangled) SED flux
         * mfunc: mangling function evaluated at [wave]
      '''
      
      if 'ks_mopts' not in self.__dict__:
         raise AttributeError("Mangling info not found... try running self.kcorr()")
      if self.ks_mopts[band][i] is None:
         return(None,None,None,None)
      if  self.k_version != 'H3+L':
         epoch = int(self.data[band].t[i]/(1+self.z)/self.ks_s)
      else:
         epoch = self.data[band].t[i]/(1+self.z)
      kextrap = getattr(self, 'k_extrapolate', False)
      wave,flux = kcorr.get_SED(int(epoch), version=self.k_version, 
            sBV=self.ks_s,extrapolate=kextrap)
      if self.ks_mopts[band][i]:
         man_flux = mangle_spectrum.apply_mangle(wave,flux, 
               **self.ks_mopts[band][i])[0]
      else:
         man_flux = flux
      if not normalize:
         return(wave*(1+self.z),man_flux,flux,man_flux/flux)

      # Now compute the normalizing factor
      num = power(10,-0.4*(self.data[band].magnitude[i] -\
            self.data[band].filter.zp))
      denom = fset[band].response(wave*(1+self.z), man_flux)
      return(wave*(1+self.z),man_flux*num/denom,flux,man_flux/flux)
   
   #def fit_color(self, band1, band2, **args):
   #   '''Fit a 1dcurve to the band1-band2 color. Can be an interactive fit
   #  as well.

   #  Args:
   #     band1,band2 (str):  Filters comprising the color
   #     interp (bool): If True, interpolate missing data in either filter
   #     use_model (bool): If True and both a model and interpolator are
   #                       defined for the filter, use the model to interpolate
   #     model_float (bool): If True, then re-fit the model to each filter
   #                         independently (so each has independent Tmax)
   #     kcorr (bool): If True, return k-corrected color.
   #     method (str):  Which interpolator to use. Use :meth`.list_types`
   #                    to find out what is available.
   #  '''
   #  interp = args.get('interp', 1)
   #  use_model = args.get('use_model', 0)
   #  model_float = args.get('model_flot', 0)
   #  kcorr = args.get('kcorr', 0)
   #  method = args.get('method', 'spline2')
   #  t,c,ec,mask = self.get_color(band1, band2, interp, use_model, model_float,
   #        kcorr)

   #  if self.Tmax > 0:
   #     evt = arange(int(self.t[0]), int(self.t[-1]+1)*1.0 + self.Tmax)
   #  else:
   #     evt = arange(self.MJD[0], self.MJD[-1])*1.0

   #  tt,cc,ecc = fit1dcurve.regularize(t, c, ec)
   #  if getattr(self, 'cinterp', None) is None:
   #     self.cinterp = {}
   #  self.cinterp[(band1,band2)] = fit1dcurve.Interpolator(method,
   #        t, c, ec, less_equal(mask,2))
   #  if interactive:

   def interp_table(self, bands, use_model=False, model_float=False, 
         dokcorr=False):
      '''For a given set of bands, construct a table of contemporaneous
      photometry, interpolating any missing data.

      Args:
         bands (list): Filters with which to construct the table
         use_model (bool): If True, use the model to interpolate and
                           allow extrapolations based on the model
         model_float (bool): If True, re-fit the model to each filter
                           independently (so each as independent Tmax)
         dokcorr (bool): If True, return k-corrected values
      
      Returns:
         4-tuple:  (MJD, mags, emags, flags)

         MJD:  the epoch of observations. This is the same for all filters,
               so is a single 1D array
         mags: list of 1D arrays of magnitudes, one for each input band
         emags: list of 1D arrays of errors, one for each intput band
         flags: list of 1D arrays, indicating how the value was derived:
                * 0 - observation was measured
                * 1 - value was based on an interpolation. Pretty safe
                * 2 - value was based on extrapolation of the model.  Meh.
                * 4 - value was obtained by extrapolation. Not safe to use.
                * 8 - k-corrections are not valid
      '''
      # First, get a table of all photometry:
      data = self.get_mag_table(bands)

      ms = []
      ems = []
      flags = []
      for band in bands:
         bids = greater(data[band], 90)   # where we need to fill in
         if dokcorr:
            if band not in self.ks_tck:
               raise RuntimeError("No k-corrections defined for %s.  Either set dokcorr=0 or run self.kcorr() first" % band)
            k = scipy.interpolate.splev(data['MJD'], self.ks_tck[band])
            kflag = (less(data['MJD'],self.ks_tck[band][0][0]) &\
                    greater(data['MJD'], self.ks_tck[band][0][-1]))*8

         else:
            k = 0
            kflag = 0

         if not any(bids):
            # They all good, carry on.
            ms.append(data[band]-k)
            ems.append(data['e_'+band])
            flags.append(bids*1 + kflag)
            continue

         # Need to interpolate from here on in
         if use_model:
            if model_float:
               temp,etemp,mask = self.model(band, self.data[band].MJD)
               weight = self.data[band].e_flux**2
               weight = power(weight, -1)*mask*self.data[band].mask
               # weighted mean offset
               offset = sum((self.data[band].mag - temp)*weight)/\
                              sum(weight)
               # error in mean offset  (MAD)
               doffset = median(absolute(self.data[band].mag - \
                                                     temp-offset))
            else:
               offset = 0
               doffset = 0
            temp,etemp,mask = self.model(band, data['MJD'])
            ms.append(where(bids, temp+offset, data[band]) - k)
            ems.append(where(bids, doffset, data['e_'+band]))
            # IDs where valid interpolation is done
            iids = mask & greater_equal(data['MJD'],self.data[band].MJD.min())\
                   & less_equal(data['MJD'], self.data[band].MJD.max())
            # IDs where valid extrapolation is done
            eids = mask & (less(data['MJD'],self.data[band].MJD.min())\
                   | greater(data['MJD'], self.data[band].MJD.max()))
            flags.append(where(bids, iids*1+eids*2+(~mask)*4, bids)+kflag)
         else:
            interp = getattr(self.data[band], 'interp', None)
            if interp is None:
               raise ValueError("Error: interpolator missing for %s" % band)
            temp,mask = self.data[band].interp(data['MJD'])
            etemp = self.data[band].interp.error(data['MJD'])
            ms.append(where(bids, temp, data[band])-k)
            ems.append(where(bids, etemp, data['e_'+band]))
            iids = mask
            eids = ~mask
            flags.append(where(bids, iids*1 + eids*4, bids)+kflag)

      return (data['MJD'], ms, ems, flags)

   def get_color(self, band1, band2, interp=True, use_model=False, 
         model_float=False, dokcorr=False, kcorr=None):
      '''return the observed SN color of [band1] - [band2].  
      
      Args:
         band1,band2 (str):  Filters comprising the color
         interp (bool): If True, interpolate missing data in either filter
         use_model (bool): If True and both a model and interpolator are
                           defined for the filter, use the model to interpolate
         model_float (bool): If True, then re-fit the model to each filter
                             independently (so each has independent Tmax)
         dokcorr (bool): If True, return k-corrected color.

      Returns:
         4-tuple:  (MJD, color, error, flag)

         MJD (float array):  
            epoch
         color (float array):  
            observed color
         error (float array):  
            uncertainty in color
         flag (int array):  
            binary flag. Each of the following conditions are logically-or'ed:

            * 0 - both bands measured at given epoch
            * 1 - only one band measured
            * 2 - extrapolation (based on template) needed, so reasonably safe
            * 4 - interpolation invalid (extrapolation or what have you)
            * 8 - One or both k-corrections invalid
      '''

      if kcorr is not None:
         warnings.warn("Use of kcorr argument is deprecated. Use dokcorr "
                       " instead", stacklevel=2)
         dokcorr=kcorr

      it1 = getattr(self.data[band1], 'interp', None)
      it2 = getattr(self.data[band2], 'interp', None)
      mod1 = band1 in self.model._fbands
      mod2 = band2 in self.model._fbands

      if not interp:
         # easy
         data = self.get_mag_table([band1, band2])
         gids = less(data[band1], 90) & less(data[band2], 90)
         mjd = data['MJD'][gids]
         col = (data[band1] - data[band2])[gids]
         ecol = sqrt(data['e_'+band1]**2 + data['e_'+band2]**2)[gids]
         flags = gids[gids]*0
         if dokcorr:
            if band1 not in self.ks:
               raise ValueError("band %s has no k-corrections, use self.kcorr()" % band1)
            if band2 not in self.ks:
               raise ValueError("band %s has no k-corrections, use self.kcorr()" % band2)
            # Now, we need to find indexes into each dataset that correspond
            #  to data['MJD']. 
            ids1 = argmin(absolute(self.data[band1].MJD[:,newaxis] - \
                  mjd[newaxis,:]), axis=0)
            ids2 = argmin(absolute(self.data[band2].MJD[:,newaxis] - \
                  mjd[newaxis,:]), axis=0)
            col = col - self.ks[band1][ids1] + self.ks[band2][ids2]
            flags = flags + (logical_not(self.ks_mask[band1][ids1])*\
                             logical_not(self.ks_mask[band2][ids2]))*8
         return (mjd,col,ecol,flags)
      elif use_model and mod1 and mod2:
         MJD,ms,ems,flags = self.interp_table([band1,band2], use_model=True,
            model_float=model_float, dokcorr=dokcorr)
      elif not use_model and it1 and it2:
         MJD,ms,ems,flags = self.interp_table([band1,band2], use_model=False,
            model_float=model_float, dokcorr=dokcorr)
      elif not use_model and mod1 and mod2:
         MJD,ms,ems,flags = self.interp_table([band1,band2], use_model=True,
            model_float=model_float, dokcorr=dokcorr)
      else:
         raise RuntimeError("You asked for interpolation, but there are no"\
               " models or interpolators defined")

      # We really don't care which one flags conditions
      flags = bitwise_or(flags[0],flags[1])

      return(MJD, ms[0]-ms[1], sqrt(ems[0]**2 + ems[1]**2), flags)

   def getEBVgal(self, calibration='SF11'):
      '''Gets the value of E(B-V) due to galactic extinction.  The ra and decl
      member variables must be set beforehand.
      
      Args:
         calibration (str):  Which MW extionction calibraiton ('SF11' or 
                             'SFD98')

      Returns:
         None

      Effects:
         self.EBVgal is set to Milky-Way color excess.
      '''
      if self.ra is not None and self.decl is not None:
         self.EBVgal,mask = dust_getval.get_dust_RADEC(self.ra, self.decl,
               calibration=calibration)
         self.EBVgal = self.EBVgal[0]
      else:
         print("Error:  need ra and dec to be defined, E(B-V)_gal not computed")

   def get_zcmb(self):
      '''Gets the CMB redshift from NED calculator and stores it locally.'''
      zcmb = getattr(self, '_zcmb', None)
      if zcmb is not None:
         return zcmb
      else:
         from .utils import zCMB
         self._zcmb = zCMB.z_cmb(self.z, self.ra, self.decl)
         return self._zcmb

   def get_distmod(self, cosmo='LambdaCDM', **kwargs):
      '''Gets the distance modulus based on a given cosmology. Requires the
      astropy module to work. 
      
      Args:
         cosmo (str): The cosmology to use. The available cosmologies are 
                      those in astropy.cosmology and kwargs can be set if the 
                      cosmology takes them (see astropy.cosmology)
                      Default:  LambdaCDM with Ho=74.3, O_m = 0.27, O_de = 0.73
         kwargs (dict): Extra arguments for given cosmology
         
      Returns:
         float: mu, the distance modulus in magnitudes.
      '''
      try:
         from astropy import cosmology
      except:
         raise ImportError("Sorry, in order to compute distance modulus, " +\
            "need the astropy module. Try 'pip install --no-deps astropy'")
      try:
         c = getattr(cosmology, cosmo)
      except:
         raise AttributeError("Unknown cosmology specified. See astropy.cosmology")
      if not isinstance(c, cosmology.core.Cosmology):
         kwargs.setdefault('H0', 74.3)
         kwargs.setdefault('Om0', 0.27)
         kwargs.setdefault('Ode0', 0.73)
         cos = c(**kwargs)
      else:
         cos = c

      return cos.distmod(self.zcmb).value

   def get_lb(self):
      '''Computes the galactic coordinates of this objects.
      
      Returns:
         2-tuple: (l,b)
         
           l:  galactic longitude (in degrees)
           b:  galactic latitude (in degrees)
      '''
      from .utils import radec2gal
      return radec2gal.radec2gal(self.ra, self.decl)

   def summary(self, out=sys.stdout):
      '''Get a quick summary of the data for this SN, along with fitted 
      parameters (if such exist).
      
      Args:
         out (str or open file): where to write the summary

      Returns:
         None
      '''
      print('-'*80, file=out)
      print("SN ",self.name, file=out)
      if self.z:  
         print("z = %.4f         " % (self.z), end=' ', file=out)
      if getattr(self, '_zcmb', None) is not None:
         print("zcmb = %.4f         " % (self.zcmb), end=' ', file=out)
      if self.ra:  print("ra=%9.5f        " % (self.ra), end=' ', file=out)
      if self.decl:  print("dec=%9.5f" % (self.decl), end=' ', file=out)
      print("", file=out)
      print("Data in the following bands:", end=' ', file=out)
      for band in self.data:  print(band + ", ", end=' ', file=out)
      print("", file=out)

      print("Fit results (if any):", file=out)
      #for band in self.restbands:
      #   print >> out, "   Observed %s fit to restbad %s" % (band, self.restbands[band])
      systs = self.systematics()
      for param in self.parameters:
         if self.parameters[param] is not None:
            line = "   %s = %.3f  +/-  %.3f" % (param, self.parameters[param],
                                                       self.errors[param])
            if param in systs and systs[param] is not None:
               line += "  +/- %.3f (sys)" % systs[param]
            print(line, file=out)

   def dump_lc(self, epoch=0, tmin=-10, tmax=70, k_correct=0, 
               s_correct=False, mw_correct=0):
      '''Outputs several files that contain the lc information: the data,
      uncertainties, and the models themselves.
      
      Args:
         epoch (bool): If True, output times relative to Tmax
         tmin/tmax (float): the time range over which to output the model
                            (default:  -10 days to 70 days after Tmax)
         k_correct (bool): If True, k-correct the data/models (default: False)
         s_correct (bool): If True, s-correct the data/models (default: False)
         mw_correct (bool): If True, de-reddent the MW extintcion 
                            (default: False)

      Returns:
         None

      Effects:
         This function will create several files with the following template 
         names:

            {SN}_lc_{filter}_data.dat
         
            {SN}_lc_{filter}_model.dat
         
         which will contain the photometric data and model for each filter (if
         that filter was fit with a model).  In the \*_data.dat, there is an
         extra flag column that indicates if the k-corrections are valid (0) or
         invalid (1).
      '''
      base = self.name + "_lc_"
      if not epoch:
         toff = 0
      else:
         toff = self.Tmax
      for filt in list(self.data.keys()):
         f = open(base+filt+"_data.dat", 'w')
         print("#  column 1:  time", file=f)
         print("#  column 2:  oberved magnitude", file=f)
         print("#  column 3:  error in observed magnitude", file=f)
         print("#  column 4:  Flag:  0=OK  1=Invalid S/K-correction", file=f)
         if s_correct:
            print("# NOTE  Data have been S-corrected from %s to %s" %\
                  (filt, self.restbands[filt]), file=f)
         if mw_correct:
            Ia_w,Ia_f = kcorr.get_SED(0, 'H3')
            Alamb = fset[filt].R(wave=Ia_w, flux=Ia_f, Rv=3.1)*self.EBVgal
         else:
            Alamb = 0
         for i in range(len(self.data[filt].mag)):
            if k_correct:
               flag = 0
               if filt not in self.ks:
                  flag = 1
                  ks = 0
               else:
                  flag = (not self.ks_mask[filt][i])
                  ks = self.ks[filt][i]
               print("%.2f  %.3f  %.3f  %d" % \
                     (self.data[filt].MJD[i]-toff, 
                     self.data[filt].mag[i] - ks - Alamb, 
                     self.data[filt].e_mag[i], flag), file=f)
            elif s_correct:
               flag = 0
               if filt not in self.Ss:
                  flag = 1
                  Ss = 0
               else:
                  flag = (not self.Ss_mask[filt][i])
                  Ss = self.Ss[filt][i]
               print("%.2f  %.3f  %.3f  %d" % \
                     (self.data[filt].MJD[i]-toff, 
                     self.data[filt].mag[i] + Ss - Alamb, 
                     self.data[filt].e_mag[i], flag), file=f)

            else:
               flag = (not self.data[filt].mask[i])
               print("%.2f  %.3f  %.3f  %d" % \
                     (self.data[filt].MJD[i]-toff, 
                     self.data[filt].mag[i] - Alamb
                     , self.data[filt].e_mag[i], flag), file=f)
         f.close()
         if filt in self.model._fbands:
            ts = arange(tmin, tmax+1, 1.0)
            ms,e_ms,mask = self.model(filt, ts+self.Tmax)
            if k_correct and filt in self.ks_tck:
               ks = scipy.interpolate.splev(ts + self.Tmax, self.ks_tck[filt])
               # mask out valid k-corrections
               mids = argmin(absolute(ts[:,newaxis]-\
                     self.data[filt].MJD[newaxis,:]+\
                     self.Tmax), axis=1)
               ks_mask = self.ks_mask[filt][mids]*greater_equal(ts, -19)*\
                     less_equal(ts, 70)
               mask = mask*ks_mask
               ms = ms - ks
            ms = ms[mask]
            e_ms = e_ms[mask]
            ts = ts[mask]
            f = open(base+filt+"_model.dat", 'w')
            print("# column 1: time", file=f)
            print("# column 2:  model magnitude", file=f)
            print("# column 3:  model magnitude error", file=f)
            for i in range(len(ts)):
               print("%.1f %.3f %.3f" % (ts[i]+self.Tmax-toff, ms[i],
                  e_ms[i]), file=f)
            f.close()
         if self.data[filt].interp is not None:
            f = open(base+filt+"_smooth.dat", 'w')
            x0,x1 = self.data[filt].interp.domain()
            ts = arange(x0, x1+1, 1.0)
            m,mask = self.data[filt].interp(ts)
            print("# column 1:  time", file=f)
            print("# column 2:  splined magnitude", file=f)
            for i in range(len(ts)):
               if not mask[i]:  continue
               print("%.1f  %.3f" % (ts[i]-toff, m[i]), file=f)
            f.close()

   def to_txt(self, filename, mask=True):
      '''Outputs the photometry to a .txt format the SNooPy can later import.
      
      Args:
         filename (str): Filename to output.
         mask (bool):    Mask out flagged data? default: True

      Returns:
         None
      '''
      f = open(filename, 'w')
      print("%s %f %f %f" % (self.name, self.z, self.ra, self.decl), file=f)
      for filt in list(self.data.keys()):
         print("filter %s" % (filt), file=f)
         for i in range(len(self.data[filt].mag)):
            # mask out bad data
            if mask and not self.data[filt].mask[i]: continue
            print("%.2f  %.3f  %.3f" % \
               (self.data[filt].MJD[i], self.data[filt].mag[i],
                self.data[filt].e_mag[i]), file=f)
      f.close()

   def to_json(self, filename, reference=None, bibcode=None, source=1):
      from . import get_osc
      '''Output a file with JSON format for this object. The format should
      conform to the Open Supernova Catalog and be suitable for upload.

      Args:
         filename (str):  Output filename
         reference (str): If specified, a string-based reference for this
                          data (e.g., 'Smith et al. 2000')
         bibcode (str): If specified, include the bibcode for this data
                       (e.g., '2017AJ....154..211K')
      Returns:
         None'''

      st = get_osc.to_osc(self, ref=reference, bibcode=bibcode, source=source)
      fout = open(filename, 'w')
      fout.write(st)
      fout.close()

   def to_mlcs(self, bands=None, sninfo='sn.info', outfile=None, Tmax=None, 
         spec_sample='hsiao', prior='rv19', vector='rv19-early-smix',
         stock=True):
      '''Output the LC data to formats suitable for fitting with MLCS2k2.
      
      Args:
         bands (list of str): List of bands to export. If None, all filters
                              that have a MLCS2k2 passband are exported.
         sinfo (str): File to update with SN info (default: sn.info)
         outfile (str): The light-curve data file to output. If None,
                         SN name + '.dat'
         Tmax (float or None): Time of maximum (initial guess for MLCS). If 
                               None, try to get from fits or splines.
         spec_sample (str): The spetroscopic sample to use for k-corrections.
                            (See MLCS documentation)
         prior (str):  Reddening prior for MLCS (see MLCS documentation)
         vector (str):  Vector to use for MLCS (see MLCS documentation)
         stock (bool):  If True, assume MLCS2k2 doesn't have CSP mods.

      Returns:
         trans (dict):  a translation between SNooPy filters and MLCS filters
                        that were exported.

      Effects:
         Creates (or appends to) an output file named sn.info with the SN
         information needed by MLCS. Also creates a file named `outfile`
         with the light-curve data. Both are used to fit in MLCS
      '''
      from .mlcs_utils import snpy_to_mlcs,snpy_to_mlcs0
      if bands is None:
         bands = self.data.keys()
      if Tmax is None:
         if self.Tmax > 1.0:
            Tmax = self.Tmax
         else:
            # Try to get max from closest filter to B
            allbands = list(self.data.keys())
            diffs = [abs(self.data[band].filter.ave_wave-fset['B'].ave_wave) \
                  for band in allbands]
            B = allbands[argmin(diffs)]
            Tmax = self.data[B].MJD[argmin(self.data[B].mag)]

      if outfile is None:
         outfile = self.name+".dat"

      # First sn.info 
      fin = open(sninfo, 'r')
      lines = fin.readlines()
      names = [line.split()[0] for line in lines]
      if self.name in names:
         idx = names.index(self.name)
      else:
         idx = None
      fin.close()

      l = "%-23s %6.4f %10.4f    %6.4f     %s    %s   %s\n" % \
            (self.name, self.z, Tmax, self.EBVgal, spec_sample, prior, vector)
      fout = open(sninfo,'w')
      if idx is not None:
         lines[idx] = l
      else:
         lines.append(l)
      for line in lines:
         fout.write(line)
      fout.close()

      fmt = "%-11s %10.4f     %5.3f    %4.3f\n"
      fout = open(outfile, 'w')
      trans = {}
      if stock:
         s2s = snpy_to_mlcs0
      else:
         s2s = snpy_to_mlcs
      for filt in bands:
         if filt not in s2s:
            print("Warning: MLCS2k2 passband not found for {}".format(filt))
            continue

         trans[s2s[filt]] = filt
         for i in range(len(self.data[filt].MJD)):
            fout.write(fmt % (s2s[filt], self.data[filt].MJD[i], 
               self.data[filt].mag[i], self.data[filt].e_mag[i]))
      fout.close()
      return trans

   def to_salt(self, bands=None, outfile=None, stock=True, flux=False):
      '''Output a LC file that can be fed into SALT.
      
      Args:
         outfile (str or None): Output file that SALT2 will input. If
                                None, default to {SN}.list
         bands (list of str): List of bands to export. All if None.
         stock (bool): Do we use the stock SALT2 filter/magsys (True)
                       or the improved CSP definitaions (False)?
         flux (bool):  If True, output in flux units rather than 
                       magnitudes.
      
      Returns:
         dict:  mapping from SNooPy filters to SALT filters
      
      Effects: 
         Creates output file with all SN info needed by SALT2 to fit.
      '''
      from .salt_utils import snpy_to_salt,snpy_to_salt0

      if outfile is None:
         outfile = self.name+".list"
      if bands is None:
         bands = self.data.keys()

      fout = open(outfile, 'w')
      fout.write('@SN %s\n' % self.name)
      fout.write('@RA %f\n' % self.ra)
      fout.write('@DEC %f\n' % self.decl)
      fout.write('@Z_HELIO %f\n' % self.z)
      fout.write('@MWEBV %f\n' % self.EBVgal)
      if flux:
         fout.write('#Date :\n#Flux :\n#Fluxerr :\n#ZP :\n#Filter :\n#MagSys :\n#end\n')
         fmt = "{:.2f} {:.6f} {:.6f} {:.6f} {}::{} {}\n"
      else:
         fout.write('#Date :\n#Mag :\n#Magerr :\n#Filter :\n#MagSys :\n#end\n')
         fmt = "{:.2f} {:.6f} {:.6f} {}::{} {}\n"
      trans = {}
      if stock:
         s2s = snpy_to_salt0
      else:
         s2s = snpy_to_salt
      for filt in bands:
         if filt not in s2s:
            print("Warning: SALT2 filter not found for {}".format(filt))
            continue
         inst,f,magsys = s2s[filt]
         trans["{}::{} {}".format(inst,f,magsys)] = filt
         l = self.data[filt]
         for i in range(len(l.mag)):
            if flux:
               fout.write(fmt.format(
                  l.MJD[i],l.flux[i],l.e_flux[i],l.filter.zp,inst,f,magsys))
            else:
               fout.write(fmt.format(
                  l.MJD[i], l.mag[i], l.e_mag[i], inst,f,magsys))
      fout.close()
      return trans

   def update_sql(self, attributes=None, dokcorr=1):
      '''Updates the current information in the SQL database, creating a new SN
      if needed.   
      
      Args:
         attributes (list or None): attributes of the SN to update. If None,
                                    then all attributes are updated.
         dokcorr (bool):  If True, also update the k-corrections in the DB.

      Returns:
         None

      Effects:
         SQL database is updated.
      '''
      if have_sql:
         N = self.sql.connect(self.name)
         if N == 0:
            self.sql.create_SN(self.ra, self.decl, self.z)
            data = {}
            for f in self.data:
               if dokcorr and data[f].K is not None:
                  data[f] = [self.data[f].MJD, self.data[f].magnitude, self.data[f].e_mag,
                          self.data[f].K]
               else:
                  data[f] = [self.data[f].MJD, self.data[f].magnitude, self.data[f].e_mag,
                          None]
            self.sql.create_SN_photometry(data)
         elif dokcorr:
            for f in self.data:
               self.sql.update_photometry(f, self.data[f].MJD, "K", self.data[f].K)
         attr_list = ['z','ra','decl'] + list(self.parameters.keys())
         for attr in attr_list:
            try:
               self.sql.set_SN_parameter(attr, self.__getattr__[attr])
            except:
               pass
         self.sql.close()

   def read_sql(self, name, spectra=False):
      '''Get the data from the SQL server for supernova.
      
      Args:
          name (str): The name to retrieve from the SQL database
          
      Returns:
         int:  -1 for failure
         
      Effects:
         If successful, the SN objects is updated with data from the 
         SQL database.
      '''
      if have_sql:
         N = self.sql.connect(name)
         if N == 0:
            print("%s not found in database, starting from scratch..." % (name))
            self.sql.close()
            return
         #try:
         self.z = self.sql.get_SN_parameter('z')
         self.ra = self.sql.get_SN_parameter('ra')
         self.decl = self.sql.get_SN_parameter('decl')
         data,source = self.sql.get_SN_photometry()
         self.sources = [source]
         for filt in data:
            d = data[filt]
            if 'K' in d:
               K = d['K']
            else:
               K = None
            if 'SNR' in d:
               SNR = d['SNR']
            else:
               SNR = None
            self.data[filt] = lc(self, filt, d['t'], d['m'], d['em'], 
                  K=K, SNR=SNR, sids=zeros(d['t'].shape, dtype=int))

         if spectra:
            from .filters import spectrum
            from .specobj import timespec
            specs = []
            MJDs = []
            ts,ws,fs = self.sql.get_SN_spectra()
            for i in range(ts.shape[0]):
               MJDs.append(ts[i] + self.sql.JD_OFFSET)
               specs.append(spectrum(wave=ws[i], flux=fs[i], fluxed=True,
                  name="Spectrum MJD={:.1f}".format(MJDs[-1])))
            self.sdata = timespec(self, MJDs, specs)
            self.sdata.sort_time()

         #finally:
         #   self.sql.close()

   def set_restbands(self):
      return self.get_restbands()

   def get_restbands(self):
      '''Automatically populates the restbands member data with one of 
      the filters supported by the currently selected model. The filter with
      the closest effective wavelength to the observed filter is selected.


      Returns:
         none

      Effects:
         self.restbands is updated with valid filter set.
      '''
      for band in self.data:
         self.restbands[band] = self.closest_band(band)

   def lc_offsets(self, min_off=0.5):
      '''Find offsets such that the lcs, when plotted, won't overlap.
      
      Args:
         min_off (float):  the minimum offset between the light-curves.

      Returns:
         list:  list of offsets, in the order in which filters are
                plotted (controled by self.filter_order)
      '''

      if self.filter_order is None:
         bands = list(self.data.keys())
         eff_wavs = []
         for filt in bands:
            eff_wavs.append(fset[filt].ave_wave)
         eff_wavs = asarray(eff_wavs)
         ids = argsort(eff_wavs)
         self.filter_order = [bands[i] for i in ids]

      offs = [0]
      filt = self.filter_order[0]

      x,y,ey = regularize(self.data[filt].MJD, self.data[filt].mag,
            self.data[filt].e_mag)
      if len(x) < 3:
         mn = y.mean()
         f = lambda x:  mn
      else:
         f = interp1d(x,y, bounds_error=False, 
               fill_value=self.data[filt].mag.max())
      for filt in self.filter_order[1:]:
         deltas = self.data[filt].mag + offs[-1] - f(self.data[filt].MJD)
         off = - deltas.max() - 0.5
         offs.append(offs[-1]+off)
         x,y,ey = regularize(self.data[filt].MJD, self.data[filt].mag,
            self.data[filt].e_mag)
         if len(x) > 1:
            f = interp1d(x,y+offs[-1], bounds_error=False,
                           fill_value=self.data[filt].mag.max()+offs[-1])
         else:
            f = lambda x:  self.data[filt].mag.mean()+offs[-1]
      return offs



   def save(self, filename):
      '''Save this SN instance to a pickle file, which can be loaded again
      using the get_sn() function.
      
      Args:
         filename (str):  output filename
         
      Returns:
         None
      '''
      f = open(filename, 'wb')
      pickle.dump(self, f)
      f.close()
   
   def fit(self, bands=None, mangle=1, dokcorr=1, reset_kcorrs=1, 
         k_stretch=True, margs={}, kcorr=None, **args):
      '''Fit the N light curves with the currently set model (see 
      self.choose_model()).  The parameters that can be varried or held 
      fixed depend on the model being used (try help(self.model)
      for this info).  If one of these parameters is specified with a 
      value as an argument, it is held fixed.  Otherwise it is varied.  
      If you set a parameter to None, it will be automoatically chosen by 
      self.model.guess().
      
      Args:
         bands (list or None):  List of observed filters to fit. If None
                                (default), fit all filters with valid
                                rest-bands.
         mangle (bool):  If True, mangle the Ia SED to fit observed colors
                         before computing k-corrections.
         dokcorr (bool):  If True, compute k-corrections as part of the fit.
         reset_kcorrs (bool):  If True, zero-out k-corrections before fitting.
         k_stretch (bool):  If True, stretch the Ia SED in time to match
                            dm15/st of the object.
         margs (dict): A set of extra arguments to send to 
                       kcorr.mangle_spectrum.mangle_spectrum2()
         args (dict): Any extra arguments are sent to the model instance
                      If an argument matches a parameter of the model,
                      that parameter will be held fixed at the specified
                      value.

      NOTE:  If you have data that has already been k-corrected (either
             outside SNooPy or by setting the individual data's K
             attributes, use dokcorr=0 and reset_kcorrs=1.  If you have
             run the self.kcorr() manually and want to keep those
             k-corrections, use dokcorr=0 and reset_kcorrs=0.  Otherwise,
             use the default dokcorr=1 and reset_kcorrs=1.

      Returns:
         None

      Effects:
         If successful, the model instance is updated with the best-fit
         values of the parameters. If self.replot is True, then a plot
         will be generated with the fit. self.ks will be filled in with
         k-corrections.
      '''

      if self.model.external_fitter:
         # We're going to skipp all stuff and just run the fitter and plot
         if bands is None:
            bands = list(self.data.keys())
         self.model.fit(bands, **args)
         if self.replot: self.plot()
         return

      if kcorr is not None:
         warnings.warn("Use of kcorr argument is deprecated. Use dokcorr "
                       " instead", stacklevel=2)
         dokcorr=kcorr

      if bands is None:
         # By default, we fit the bands whose restbands are provided by the model
         bands = [b for b in list(self.data.keys()) \
               if self.restbands[b] in self.model.rbs]

      # Setup initial Robs (in case it is used by the model)
      for band in bands:
         if band not in self.Robs:
            self.Robs[band] = fset[band].R(self.Rv_gal, Ia_w, Ia_f, z=self.z)

      #if self.z <= 0:
      #   raise ValueError("The heliocentric redshift is zero.  Fix this before you fit")

      # Check to make sure we have filters we can fit:
      for filt in bands:
         if self.restbands[filt] not in self.model.rbs:
            raise AttributeError("Error:  filter %s is not supported by this model" % self.restbands[filt] + \
                  ", set self.restbands accordingly")

      if reset_kcorrs:
         self.ks = {}
         self.ks_mask = {}
         self.ks_tck = {}
      if not self.quiet:
         print("Doing Initial Fit to get Tmax...")
      self.model.fit(bands, **args)


      if dokcorr:
         kbands = [band for band in bands if band not in self.ks]
         if len(kbands) > 0:
            if not self.quiet:
               print("Setting up initial k-corrections")
            self.kcorr(kbands, mangle=0, use_stretch=k_stretch)
 
         if not self.quiet:
            if mangle:
               print("Doing first fit...")
            else:
               print("Doing fit...")
         self.model.fit(bands, **args)
 
         if mangle:
            if not self.quiet:
               print("Doing mangled k-corrections")
            self.kcorr(bands, interp=0, use_model=1, use_stretch=k_stretch, **margs)
            if not self.quiet:
               print("Doing final fit...")
            self.model.fit(bands, **args)
      if self.replot:
         self.plot()

   def fitMCMC(self, bands=None, Nwalkers=None, threads=1, Niter=500, 
         burn=200, tracefile=None, verbose=False, plot_triangle=False,
         **args):
      '''Fit the N light curves of filters specified in [bands]
      with the currently set model (see 
      self.choose_model()) using MCMC. Note that this function requires
      the emcee module to do the sampling. You should do an initial fit
      using the regular least-squares fit() function to compute good
      K-corrections and also get a starting point. The parameters that can be 
      varried or held fixed depend on the model being used (try
      help(self.model) for this info).  Parameters can be given priors by
      setting values to them as arguments (e.g., Tmax=0). This can be done in
      one of several ways: 

      * specify a floating point number. The parameter is held fixed at 
        that value.  
      * specify a shorthand string: 
         * U,a,b: Uniform prior with lower/upper limits equal to a/b 
         * G,m,s:  Gaussian prior with mean m and std dev. s.  
         * E,t:    Exponential positive prior: p = exp(-x/tau)/tau x > 0 
      * specify a function that takes a single argument:  the value of the
        parameter and returns thei the prior as a log-probability.

      Args:
         bands (list or None):  a list of observed filters to fit. If None,
                               all filters with valid rest-frame filters are
                               fit.
         Nwalkers (int or None):  Number of emcee walkers to spawn. 
                                   see emcee documentation.
         threads (int):  Number of threads to spawn. Usually not very useful
                         as overhead requires more CPU than computing the model.
         Niter (int): Number of interations to run per walker.
         burn (int): burn-in iterations
         tracefile (str):  Optional name of a file to which the traces of the
                           MCMC will be stored.
         verbose (bool): be verbose?
         plot_triangle (bool): If True, plot a covariance plot. This requires
                               the triangle_plot module (get it from pypi).
         args (dict):  Any extra arguments are sent to the model instance's
                       fit() function. Note that if an argument matches a
                       parameter name, it is treated specially as described
                       above.
      '''
      if snemcee is None:
         print("Sorry, in order to fit with MCMC sampler, you need to install")
         print("the emcee module. Try 'pip install emcee'  or get the source")
         print(" from http://dan.iel.fm/emcee/current/") 
         return None

      if len(self.model._fbands) == 0:
         raise AttributeError("In order to fit with the MCMC sampler, you need to do \
              an initial fit first.")

      if bands is None:
         # By default, we fit the bands whose restbands are provided by the model
         bands = [b for b in list(self.data.keys()) \
               if self.restbands[b] in self.model.rbs]
      if verbose:
         print("Fitting "," ".join(bands))

      Nparam = len(list(self.model.parameters.keys()))
      if Nwalkers is None:
         Nwalkers = Nparam*10
      if verbose:
         print("Setting up %d walkers..." % Nwalkers)
      
      # Set up the inverse covariance matrices and determinants
      #self.invcovar = {}
      #self.detcovar = {}
      self.bcovar = {}
      for band in bands:
         # in fluxes
         thiscov = self.data[band].get_covar(flux=1)
         #thiscov = self.data[band].get_covar(flux=0)
         rband = self.restbands[band]
         if len(thiscov.shape) == 1:
            thiscov = diagflat(thiscov)
         #if self.model.model_in_mags:
         #   modcov = self.model.get_covar(rband, self.data[band].t)
         #   #modcov = modcov*outer(self.data[band].flux, self.data[band].flux)*1.087**2
         #   #thiscov = thiscov + modcov
         #else:
         #   #thiscov = thiscov + self.model.get_covar(rband, self.data[band].t)
         #   #raise RuntimeError, "model must be in mags"
         ##self.invcovar[band] = linalg.inv(thiscov)
         self.bcovar[band] = thiscov

      # This step needs to be done because we are bypassing the usual
      # model setup
      self.model.args = args.copy()

      sampler,vinfo,p0 = snemcee.generateSampler(self, bands, Nwalkers, threads,
            tracefile, **args)
      if verbose:
         print("Doing initial burn-in of %d iterations" % burn)
      if burn > 0:
         pos,prob,state = sampler.run_mcmc(p0, burn)
         if verbose:
            print("Now doing production run of %d iterations" % Niter)
         sampler.reset()
      else:
         pos = p0
      pos,prob,state = sampler.run_mcmc(pos, Niter)

      # The parameters in order of the sampler
      pars = []
      samples = []
      for par in vinfo:
         if type(vinfo[par]) is type({}) and 'index' in vinfo[par] \
               and vinfo[par]['prior_type'] != 'nuissance':
            pars.append(par)
            samples.append(sampler.flatchain[:,vinfo[par]['index']])
      samples = array(samples).T

      # Save the tracefile as requested
      if tracefile is not None:
         d = dict(samples=samples, vinfo=vinfo, pos=pos, prob=prob, state=state,
               pars=pars)
         f = open(tracefile+"_full", 'wb')
         pickle.dump(d, f)
         f.close()
         f = open(tracefile, 'w')
         for i in range(len(pars)):
            f.write('# Col(%d) = %s\n' % (i+1,pars[i]))
         savetxt(f, samples, fmt="%15.10g")
         f.close()
      # Now that we have the samples, we can infer the median and covariance
      meds = median(sampler.flatchain, axis=0)
      covar = cov(sampler.flatchain.T)
      
      self.model.C = {}
      for par in self.model.parameters:
         if not vinfo[par]['fixed']:
            ind = vinfo[par]['index']
            self.model.parameters[par] = meds[ind]
            self.model.errors[par] = sqrt(covar[ind,ind])
            self.model.C[par] = {}
            for par2 in self.model.parameters:
               if not vinfo[par2]['fixed']:
                  ind2 = vinfo[par2]['index']
                  self.model.C[par][par2] = covar[ind,ind2]
      if self.replot:
         self.plot()

      if plot_triangle:
         if triangle is None:
            print("Sorry, but if you want a triangle plot, you have to install")
            print("the corner module (https://github.com/dfm/corner.py)")
         else:
            triangle.corner(samples, labels=pars, truths=meds)

   def systematics(self, **args):
      '''Report any systematic errors that may be present in the 
      fit parameters.

      Args:
         args (dict):  All arguments are sent to the model.systmatics()
                       function.

      Returns:
         dict:  a dictionary of systematic errors keyed by parameter
                name.  It therefore depends on the model being used.  Also
                see the specific model for any extra arguments.  If None
                is returned as a value, no systematic has been estimated
                for it.'''
      return self.model.systematics(**args)

   def plot_filters(self, bands=None, day=0, outfile=None, fill=False):
      '''Plot the filter functions over a typical SN Ia SED.

      Args:
         bands (list or None): filters to plot or, if None, all observed
                               filters.
         day (int): Which epoch (t-tmax) to use for retrieving the Ia SED
         outfile (str): optional filename for graph output.
         fill (bool): If True, use matplotlib's fill_between to fill in
                      filter and SED curves. Useful to gauge overlap.

      Returns:
         matplotlib.figure:  the figure instance for the plot.
      '''
      return plotmod.plot_filters(self, bands, day, outfile=outfile)


   def plot_color(self, f1, f2, epoch=True, deredden=True, interp=False,
         dokcorr=False, kcorr=None, outfile=None, clear=True):
      '''Plot the color curve (color versus time) given by f1 and f2.

      Args:
         f1,f2 (str):  The two filters defining the color (f1-f2)
         epoch (bool):  If True, plot time relative to Tmax
         deredden (bool): If True, remove Milky-Way foreground reddening.
         dokcorr (bool): If True, apply and defined K-corrections to the colors
         outpfile (str): Optional name for graph output to file.

      Returns:
         matplotlib.figure:  the figure istance with the plot.
      '''
      if kcorr is not None:
         warnings.warn("Use of kcorr argument is deprecated. Use dokcorr "
                       " instead", stacklevel=2)
         dokcorr=kcorr

      return plotmod.plot_color(self, f1,f2,epoch, deredden, interp, dokcorr, 
            outfile, clear)

   def compute_w(self, band1, band2, band3, R=None, Rv=3.1, interp=False,
                 use_model=False, model_float=True, dokcorr=False):
      '''Returns the reddening-free magnitude (AKA Wesenheit function)
      in the sense that:
      w = band1 - R(band1,band2,band3)*(band2 - band3)
      for instance compute_w(V,B,V) would give:
      w = V - Rv(B-V)
      
      Args:
         band1,band2,band3 (str): The three filters defining w
         R (float or None):  the R parameter to use (if None,
                            compute assuming reddening due to dust
                            with reddening law R_V = Rv
         Rv (float or None): If specified and R is None, use this
                            Rv to compute R(\lambda) for filters
         interp (bool):  If True, interpolate missing data
         use_model (bool):  If True, use a fit model, rather than an
                         interpolator to do the interpolation
         model_float (bool):  If using model interpolation, allow the
                         model for each filter to float to an indepenent
                         maximum.
         dokcorr (bool):  If True, apply k-corrections to the data.

      Returns:
         (w,ew,flag):  w:  Wesenheit magnitude
                      ew:  error in w
                    flag:  bitwise OR of the following conditions:
                           0:  all magnitudes measured at every epoch
                           1:  some magnitudes interpolated (safely)
                           2:  some magnitudes extrapolated based on model
                           4:  some magnitudes extrapolated (not safe)
                           8:  some k-corrections not valid
      '''
      # First, let's get the proper value of R:
      if R is None:
         R1 = fset[band1].R(Rv, Ia_w, Ia_f)
         R2 = fset[band2].R(Rv, Ia_w, Ia_f)
         R3 = fset[band3].R(Rv, Ia_w, Ia_f)
         R = R1/(R2 - R3)

      if dokcorr:
         if band1 not in self.ks or band2 not in self.ks or \
               band3 not in self.ks:
            raise ValueError("Not all requested bands have k-corrections.  "+\
                              "Try running self.kcorr()")
      # Now, we need the values. First, if not interpolating, just use the 
      # data
      if not interp:
         data = self.get_mag_table([band1,band2,band3])
         gids = less(data[band1],90) & less(data[band2],90) & \
               less(data[band3],90)
         if dokcorr:
            k1 = self.ks[band1][gids]
            k2 = self.ks[band2][gids]
            k3 = self.ks[band3][gids]
         else:
            k1 = k2 = k3 = 0

         mjd = data['MJD']
         m1 = data[band1][gids] - k1; em1 = data['e_'+band1][gids]
         m2 = data[band2][gids] - k2; em2 = data['e_'+band2][gids]
         m3 = data[band3][gids] - k3; em3 = data['e_'+band3][gids]

         w = m1 - R*(m2 - m3)
         ew = sqrt(em1**2 + R**2*(em2**2 + em3**2))
         return mjd,w,ew,logical_not(isnan(w))

      # Here we need interpolation
      mjd,ms,ems,flags = self.interp_table([band1,band2,band3], 
            use_model=use_model, model_float=model_float, dokcorr=dokcorr)

      # Now compute w:
      w = ms[0] - R*(ms[1] - ms[2])
      ew = sqrt(ems[0]**2 + R**2*(ems[1]**2 + ems[2]**2))
      flags = flags[0] | flags[1] | flags[2]
      return(mjd,w,ew,flags)

   def mask_data(self):
      '''Interactively mask out bad data and unmask the data as well.  The only
      two bindings are "A" (click):  mask the data and "u" to unmask the data.
      '''
      return plotmod.mask_data(self)

   def mask_epoch(self, tmin, tmax):
      '''Mask out data in a time range  (relative to B maximum) for all 
      filters.
      
      Args:
         tmin, tmax (float): Time range over which to mask (good) data.
      
      Returns:
         None
         
      Effects:
         the mask attribute for every lc instance in self.data is udpated.
      '''
      for f in self.data:
         self.data[f].mask_epoch(tmin, tmax)

   def mask_emag(self, emax):
      '''Mask out data with (magnitude) error larger than [emax].
      
      Args:
         emax (float):  maximum error allowed. All others masked out.
      
      Returns:
         None
         
      Effects:
         the mask attribute for every lc instance in self.data is udpated.
      '''
      for f in self.data:
         self.data[f].mask_emag(emax)

   def mask_SNR(self, minSNR):
      '''Mask out data with signal-to-noise less than [minSNR].
      
      Args:
         minSNR (float):  minimum signal-to-noise ratio needed for good data.
      
      Returns:
         None
         
      Effects:
         the mask attribute for every lc instance in self.data is udpated.
      '''
      for f in self.data:
         self.data[f].mask_SNR(minSNR)
   
   def plot(self, **kwargs):
      '''Plot out the supernova data in a nice format.  There are many
      options for controlling the output.

      Args:
         xrange (2-tuple): specify the x (time) range to plot (xmin,xmax)
         yrange (2-tuple): specify the y (mag/flux) range to plot (ymin,ymax)
         title (str):  optional title
         single (bool):  If True, plot out as a single (rather than panelled)
                         plot with each filter a separate data set.
         offset (bool):  If True, offset the lightcurves (for single plots) by
                         constant amount such that they don't cross.
         legend (bool): If True and single=True, plot the legend.
         fsize (float):  override the font size used to plot the graphs
                         (default: 12)
         linewidth (int):  override the line width (default: 1)
         symbols (dict):  dictionary of symbols, indexed by filter name.
         colors (dict):  dictionary of colors to use, indexed by band name.
         relative (bool):  If True, plot only relative magnitudes 
                           (normalized to zero). Default: False
         mask (bool):  If True, omit plotting masked data.
         label_bad (bool):  If True, label the masked data with red x's?
         Nxticks (int): maximum number of x-axis tick marks (default: MPL auto)
         JDoffset (bool): If true, compute a JD offset and put it in the x-axis
                          label (useful if x-labels are crowded) default: False
         flux (bool): If True, plot in flux units
         epoch (bool): If True, plot time relative to Tmax
         filtsep (float): If not None, filters whose average wavelength differ
                          by filtsep will be plotted in the same panel.
         combfilt (bool):  If true, filters starting with the same letter will
                          be plotted in the same panel (e.g., g, g2 g_p).
         outfile (str):  optional file name to save the plot
         plotmodel (bool): if True and both a model and spline are present,
                           plot the model instead of spline.
      '''

      return plotmod.plot_sn(self, **kwargs)

   def splot(self):
      from .specobj import InteractiveSpec
      '''Launch the spectrum plotter.'''
      if getattr(self, 'sdata', None) is None:
         raise RuntimeError("Error: this sn object has no spectral data")
      sp = InteractiveSpec(self.sdata)
      return sp

   def plot_kcorrs(self, colors=None, symbols=None, outfile=None):
      '''Plot the derived k-corrections after they have been computed.
      Both mangled and un-mangled k-corrections will be plotted as
      lines and points, respectively.  If mangling was used to
      do the k-corrections, clicking 'm' on a point will bring up
      another plot showing the original and mangled spectrum. 
      
      Args:
         colors (dict or None): specify colors to use by giving a dictionary
                                keyed by filter name
         symbols (dict or None): specfiy symbols to use by givein a dictionary
                                 keyed by filter name
         outfile (str): optional file name for outputting the graph to disk.


      Returns:
         matplotlib.figure:  the figure instance of the plot.
      '''
      return plotmod.plot_kcorrs(self, colors, symbols)

   def bolometric(self, bands, method="direct", DM=None, EBVhost=None,
         Rv=None, redlaw=None, outfile=None, extra_output=False, 
         cosmo='LambdaCDM', **args):
      '''
      Produce a quasi-bolometric flux light-curve based on the input [bands]
      by integrating a template SED from \lambda=lam1 to \lambda=lam2.

      Args:
         bands (list):  list of filters to constrain the bolometric flux
         method (str):  The method to use: 'direct' or 'SED'. See
                        SNooPy documentation for how these methods differ
         lam1,lam2 (float):  limits of the integration. This is igored for
                             the direct method.
         refband (str or None):  For the SED method, the band to use to
                             normalize the SED. If you have more confidence
                             in one particular filter (or want to avoid
                             automatically choosing a bad one), specify this.
                             Ignored by the direct method.
         EBVhost (float or None): The value of E(B-V) to use for the host
                             galaxy. If None, will look for model parameter
                             in the supernova object.
         Rv (float or None):  Value of R_V to use for the host galaxy. If None,
                             will look for parameter in the supernova object.
         redlaw (str): The reddening law to use ('ccm','fm','etc). If None,
                             will look up in sn model.
         extrap_red (str or None):  Specify how we extrapolate to 
                                     \lambda -> lam1 'RJ' specifies 
                                     Rayleigh-Jeans. None turns off 
                                     extrapolation in the red.
         Tmax (float or None): Used by SED method to set the epoch of maximum.
                               If None, will look for it as parameter in model.
         interpolate (str or None):  if not None, missing photometry for each
                              epoch will be interpolated. If 'spline', it will
                              use a spline, if 'model', it will use the
                              currently fit model.
         extrapolate (bool): If interpolating by model, you can choose to
                        extrapolate the model beyond the data. Default: False
         mopts (dict):  You can pass any parameters for mangle_spectrum2()
                        using this option. Only used by SED method.
         SED (str, tuple, function): Specify the SED to use. This can be a
                        string ("H3","H","N", or "91bg"), a (wave,flux)
                        tuple, or a python function f(t) that, given time t,
                        returns a (wave,flux) tuple. For the SED method, this
                        is the SED that will be fluxed, then integrated. For
                        the direct method, it is used to compute effective
                        wavelenth of the filter and convert mag -> flux.
         DM (float or None): The distance modululs. If None, it will be computed
                        using the supernova's z_cmb and the specified cosmology
                        (see cosmo argument below).
         cosmo (string): The cosmology to use. See astropy.cosmology.
         use_stretch(bool): If dm15 or stretch are defined, use this to
                        de-stretch the SED template? Ignored by direct
                        method.
         verbose (bool):  Be verbose?
         outfile (str): If not None, name of output file for bolometric
                        luminosities.
         extra_output(bool): If True, return a dictionary as extra element
                             (see returns)

      Returns:
         3 or 4-tuple:  (epoch, bolo, filters_used, limits)
                 epoch: array of epochs
                 bolo: array of bolometric luminosities in erg/s
                 filters_used: a list of filters used for each epoch
                 limits:  2-type of wavelength limits for each epoch

         If extra_output is True, also return a dict as 4th element
                 containing the following keys:
                 mflux:  the mangled flux (or inferred flux)
                 mfuncs: the mangling functions
                 mwave:  the wavelength vector or effective wavelenghts
                 mags:   the magnitudes to which we mangled the SED
                 masks:  the masks showing where interpolation was done
         
      '''

      if DM is None:
         DM_computed = True
         DM = self.get_distmod(cosmo=cosmo)
      else:
         DM_computed = False

      if method == 'direct':
         fargs = {'sn':self, 'bands':bands, 'DM':DM, 'EBVhost':EBVhost,
               'Rv':Rv, 'redlaw':redlaw, 'extra_output':False}
         for key in ['tmin','tmax','interpolate', 'SED','Tmax']:
            fargs[key] = None
         fargs['extrap_red'] = 'RJ'
         fargs['extrapolate'] = False
         fargs['interp_all'] = False
         fargs['cosmo'] = 'LambdaCDM'
         fargs['verbose'] = False
         fargs['extra_output'] = False
         myupdate(fargs,args)
         res = bolometric.bolometric_direct(**fargs)
         limits = [(l.min(),l.max()) for l in res['lam_effs']]

      else:
         fargs = {'sn':self, 'bands':bands, 'DM':DM, 'EBVhost':EBVhost,
               'Rv':Rv, 'redlaw':redlaw, 'extra_output':False}
         for key in ['lam1','lam2','refband','tmin','tmax',
               'interpolate', 'SED','Tmax']:
            fargs[key] = None
         fargs['extrap_red'] = 'RJ'
         fargs['extrapolate'] = False
         fargs['mopts'] = {}
         fargs['SED'] = 'H3'
         fargs['interp_all'] = False
         fargs['cosmo'] = 'LambdaCDM'
         fargs['use_stretch'] = True
         fargs['extrap_SED'] = True
         fargs['verbose'] = False
         fargs['extra_output'] = False
         myupdate(fargs,args)
         res = bolometric.bolometric_SED(**fargs)
         if len(shape(res['boloflux'])) == 1:
            limits = [(l.min(),l.max()) for l in res['waves']]
         else:
            limits = []
            for i in range(res['waves'].shape[0]):
              limits.append([(l.min(),l.max()) for l in res['waves'][i]])

      if outfile is not None:
         fout = open(outfile, 'w')
         fout.write('# Bolometric luminosity for %s\n' % (self.name))
         fout.write('# Using %s method\n' % (method))
         fout.write('# E(B-V) = %.2f, R_V = %.2f, redlaw=%s\n' % \
               (EBVhost,Rv,redlaw))
         if not DM_computed:
            fout.write('# DM = %.2f\n' % (DM))
         else:
            fout.write(\
                  '# DM = %.2f  (zcmb = %.4f, H0 = 74.3, and %s cosmo)\n' %\
                  (DM, self.zcmb, cosmo))
         fout.write('# Luminosities in units of 10^42 erg/s\n')
         fout.write('# Times are in *observed* frame\n')
         if len(shape(res['boloflux'])) == 2:
            fout.write('#\n#%6s ' % 't')
            for j in range(res['boloflux'].shape[1]):
               fout.write('%8s %8s %6s ' % ('lam1','lam2','L'))
            fout.write('%s\n' % 'filters')
            for i in range(res['epochs'].shape[0]):
               fout.write("%7.2f " % res['epochs'][i])
               for j in range(res['boloflux'].shape[1]):
                  fout.write("%8.1f %8.1f %6.3f " %\
                        (limits[i][j][0], limits[i][j][1],
                         res['boloflux'][i,j]/1e42))
               fout.write("%s\n" % ("".join(res['filters_used'][i])))
         else:
            fout.write('#\n#%6s  %9s  %8s %8s %s\n' % ('t','Lbol','lam1','lam2',
               'filters'))
            for i in range(len(res['epochs'])):
               fout.write('%7.2f  %6.3f  %8.1f %8.1f %s\n' % \
                     (res['epochs'][i], res['boloflux'][i]/1e42, 
                      limits[i][0], limits[i][1],
                      "".join(res['filters_used'][i])))
         fout.close()
      if not extra_output:
         return (res['epochs'],res['boloflux'],res['filters_used'],limits)
      if method == 'direct':
         extra = dict(mflux=res['fluxes'], mwaves=res['lam_effs'],
            mags=res['mags'], masks=res['masks'])
      else:
         extra = dict(mflux=res['fluxes'], mwaves=res['waves'],
            mags=res['mags'], masks=res['masks'],
            mfuncs=res['mfuncs'], pars=res['pars'], refbands=res['refbands'])
      return (res['epochs'],res['boloflux'],res['filters_used'],limits,extra)


   def closest_band(self, band, tempbands=None, lowz=0.15):
      '''Find the rest-frame filter in [tempbands] that is closest to the
      observed filter [band].  If tempbands is None, defaults to 
      self.template_bands.  In the case where the redshift of the
      SN is below [lowz], if [band] is in [tempbands], use [band]
      regardless of whether another band is closer.'''
      if tempbands is None:
         tempbands = self.template_bands
      if self.z < lowz and band in tempbands:  return band

      resps = []
      # normalize responses to the area under the filter response curve
      norm = fset[band].response(fset[band].wave, fset[band].wave*0.0+1.0, z=0,
            zeropad=1, photons=0)
      for temp in tempbands:
         norm2 = fset[temp].response(fset[temp].wave, fset[temp].wave*0.0+1.0, z=0,
                           zeropad=1, photons=0)
         resps.append(fset[band].response(fset[temp].wave, fset[temp].resp,
            z=self.z, zeropad=1, photons=0)*norm/norm2)

      resps = array(resps)
      if max(resps) <= 0:
         # all failed to overlap...  
         dists = absolute(array([fset[temp].ave_wave - fset[band].ave_wave \
               for temp in tempbands]))
         return tempbands[argmin(dists)]
      else:
         return(tempbands[argmax(resps)])

def save(instance, file):
   '''Save a super instance to a file to be loaded back later with load().'''
   f = open(file,'w')
   pickle.dump(instance, f)
   f.close()

def load(fname):
   '''Load a pickle file which has a SNooPy instance stored. If this can't
   be done, return None and issue an error message. This function also tries
   to migrate previous saves.'''
   try:
      f = open(fname, 'rb')
      if six.PY3:
         inst = pickle.load(f, encoding='iso-8859-1')
      else:
         inst = pickle.load(f)
      f.close()
      if getattr(inst, 'k_extrapolate', None) is None:
         inst.k_extrapolate = False
   except:
      print("Warning:  could not load pickle file {}, trying as plain "\
            "text import".format(fname))
      inst = None
   return(inst)

def dump_arrays(file, arrays, formats=None, labels=None, separator=' '):
   f = open(file, 'w')
   if formats is None:
      formats = ["%11.5f"]*len(arrays)
      widths = [11]*len(arrays)
   else:
      widths = [len(form % (0.0)) for form in formats]

   if labels is not None:
      forms = ["%%%ds" % (wid) for wid in widths]
      header = [forms[i] % (labels[i]) for i in range(len(labels))]
      header = separator.join(header)
      header[0] = "#"
      print(header, file=f)

   for i in range(len(arrays[0])):
      line = [formats[j] % (arrays[j][i]) for j in range(len(arrays))]
      print(separator.join(line), file=f)
   f.close()

def import_lc(file):
   '''Import SN data from a datafile in the following format:
   line 1:     name z ra decl
   line 2:     filter {filter name}
   line 3:     Date   magnitude  error
   ...
   line N:     filter {filter name}
   line N+1:   Date   magnitue   error
   ....
   '''
   if type(file) is type(""):
      f = open(file)
   else:
      f = file
   lines = f.readlines()
   fields = lines[0].split()
   if len(fields) != 4:  raise RuntimeError("first line of %s must have 4 " +\
         "fields:  name, redshift, RA, DEC")
   name = fields[0]
   try:
      z,ra,decl = list(map(float, fields[1:]))
   except:
      raise RuntimeError("z, ra and dec must be floats " + \
            " (ra/dec in decimal degrees)")


   s = sn(name, ra=ra, dec=decl, z=z)
   lines = lines[1:]
   this_filter = None
   MJD = {}
   mags = {}
   emags = {}
   sids = {}
   sources = []

   for line in lines:
      if line[0] == "#":  continue
      if line.find('filter') >= 0:
         fields = line.split()
         this_filter = fields[1]
         if len(fields) == 2:
            if this_filter in ['u','g','r','i','B','V','Y','J','H','K',
                               'Ydw','Jrc2']:
               source = 'CSP'
            else:
               source = "unknown"
         else:
            source = " ".join(fields[2:])
         if source not in sources:
            sources.append(source)
         if this_filter not in MJD:
            MJD[this_filter] = []
            mags[this_filter] = []
            emags[this_filter] = []
            sids[this_filter] = []
      elif this_filter is not None:
         try:
            #t,m,em = list(map(float, string.split(string.strip(line))))
            t,m,em = map(float, line.strip().split())
         except:
            raise RuntimeError("Bad format in line:\n %s" % (line))
         MJD[this_filter].append(t)
         mags[this_filter].append(m)
         emags[this_filter].append(em)
         sids[this_filter].append(sources.index(source))
   # Save the sources of photometry
   s.sources = sources
   for f in MJD:
      MJD[f] = array(MJD[f])
      mags[f] = array(mags[f])
      emags[f] = array(emags[f])
      sids[f] = array(sids[f])
      s.data[f] = lc(s, f, MJD[f], mags[f], emags[f], sids=sids[f])
      s.data[f].time_sort()

   s.get_restbands()

   return(s)

def get_sn(str, sql=None, **kw):
   '''Attempt to get a sn object from several possible sources.  First, if str
   corresponds to an existing file name, the function attempts to load the sn
   instance from it as if it were a pickle'd object.  If that fails, it attempts
   to use import_lc() on the file.  If str is not the name of an existing file,
   it is treated as a SN name and is retrieved from the designated sql connection
   object (or default_sql if sql=None), in which case all keyword arguments
   are sent as options to the sql module.'''
   if os.path.isfile(str):
      s = load(str)
      if s is None:
         try:
            s = import_lc(str)
         except RuntimeError:
            raise RuntimeError("Could not load %s into SNPY" % str)
   elif str.find('http://') == 0 or str.find('https://') == 0 or \
         str.find('osc:') == 0:
      from . import get_osc
      s,message = get_osc.get_obj(str, **kw)
      if s is None:
         print(message)
         return None
   else:
      s = sn(str, source=sql, **kw)
   return s

def check_version():
   global __version__
   if six.PY3:
      import urllib.request as urllib
   else:
      import urllib
   #from distutils.version import LooseVersion
   try:
      from packaging.version import parse as LooseVersion
   except:
      from distutils.version import LooseVersion
   try:
      u = urllib.urlopen(
            'ftp://ftp.obs.carnegiescience.edu/pub/cburns/snpy/latest',
            timeout=1)
   except:
      return None,None
   if not u:
      return None,None
   lines = u.readlines()
   ver = lines[0].strip()
   u.close()
   if six.PY3: ver = str(ver, 'utf-8')
   ver = LooseVersion(ver)
   if ver > LooseVersion(__version__):
      return True,lines[1:]
   else:
      return False,None
