#!/usr/bin/env python

import sys,string,os,re
from glob import glob
have_sql = 1
import sqlmod
if 'SQLSERVER' in os.environ:
   sqlmod.default_sql = sqlmod.__dict__['sql_'+os.environ['SQLSERVER']]()
else:
   sqlmod.default_sql = sqlmod.sql_highz()
have_sql = sqlmod.have_sql
   
import types
import plotmod
from lc import lc           # the light-curve class
from distutils.version import StrictVersion
from numpy.oldnumeric import *       # Vectors
import ubertemp             # a template class that contains these two
import kcorr                # Code for generating k-corrections
import utils.IRSA_dust_getval as dust_getval

from utils import stats     # some convenient stats functions
from utils import fit_poly  # polynomial fitter
import scipy                # Scientific python routines
linalg = scipy.linalg       # Several linear algebra routines
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from utils import fit_spline # My Spline fitting routines
from filters import fset    # filter definitions.
from filters import standards as spectra # spectra.
import mangle_spectrum      # SN SED mangling routines
import pickle
import model
from utils.fit1dcurve import list_types

Version = '0.7'     # Let's keep track of this from now on.

# Some useful functions in other modules which the interactive user may want:
getSED = kcorr.get_SED
Robs = kcorr.R_obs
Ia_w,Ia_f = getSED(0, 'H3')

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
      for key in self.parent.data.keys():
         ret += "%s -> %s, " % (key, self.__getitem__(key))
      return ret
   def __repr__(self):
      return self.__str__()

   def keys(self):
      return self.dict.keys()

class sn(object):
   '''This class is the meat of the program.  Create a supernova object by 
   calling the constructor with the name of the superova as the argument.  
   E.g:
   In[1]:  s = sn('SN1999T')
   if the supernova is in the SQL database, its lightcurve data will be loaded 
   into its member data.  Once the object is created, use its member data 
   and functions to do your work.  Of course, you can have multiple 
   supernovae defined at the same time.'''

   def __init__(self, name, source=None, ra=None, dec=None, z=0):
      '''Create the object.  Only required parameter is the [name].  If this 
      is a new object, you can also specify [ra], [dec], and [z].'''
      self.__dict__['data'] = {}        # the photometric data, one for each band.
      self.__dict__['model'] = model.EBV_model(self)
      self.template_bands = ['u','B','V','g','r','i','Y','J','H','K']

      self.Version = Version    # A version-stamp for when we upgrade
      self.name = name
      self.z = z              # Redshift of the SN
      self.ra = ra              # Coordinates
      self.decl = dec
      self.filter_order = None  # The order in which to plot the filters
      self.xrange = None 
      self.yrange = None        # Impose any global plotting ranges?
      #self.Rv_host = 3.1        # Ratio of total to selective absorption in 
      #                          # host and
      self.Rv_gal = 3.1         #  galaxy
      self.EBVgal = 0.0
      self.fit_mag = False     # fit in magnitude space?

      self.restbands = dict_def(self, {})   # the band to which we are fitting for each band
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
               self.read_sql(self.name)
            else:
               print "Warning:  ra and/or decl not specified and no source specified."
               print "   setting ra, decl and z = 0" 
               self.ra = 0;  self.decl = 0;
         else:
            self.ra = ra
            self.decl = dec
      else:
         self.sql = source
         self.read_sql(self.name)

      #self.summary()
      self.getEBVgal()
      self.get_restbands()     # based on z, assign rest-frame BVRI filters to 
                               # data
      self.k_version = 'H3'

   def __getattr__(self, name):
      if 'data' in self.__dict__:
         if name in self.data:
            return(self.data[name])

      if name == 'Tmax':
         if 'model' in self.__dict__:
            if 'Tmax' in self.__dict__['model'].parameters:
               if self.__dict__['model'].parameters['Tmax'] is not None:
                  return self.__dict__['model'].parameters['Tmax']
         for f in self.data:
            if self.restbands[f] == 'B':
               if self.data[self.restbands[f]].Tmax is not None:
                  return self.data[self.restbands[f]].Tmax
         return 0.0

      if name == 'dm15':
         if 'model' in self.__dict__:
            if 'dm15' in self.__dict__['model'].parameters:
               if self.__dict__['model'].parameters['dm15'] is not None:
                  return self.__dict__['model'].parameters['dm15']
         for f in self.data:
            if self.restbands[f] == 'B':
               if self.data[self.restbands[f]].dm15 is not None:
                  return self.data[self.restbands[f]].dm15
         return None
         
      if name == 'parameters':
         if 'model' in self.__dict__:
            return self.__dict__['model'].parameters
         else:
            raise AttributeError, "Error, model not defined, so no paramters"
      if name == 'errors':
         if 'model' in self.__dict__:
            return self.__dict__['model'].errors
         else:
            raise AttributeError, "Error, model not defined, so no errors"
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
      raise AttributeError, "Error:  attribute %s not defined" % (name)

   def __setattr__(self, name, value):
      if 'model' in self.__dict__:
         if name in self.__dict__['model'].parameters:
            self.__dict__['model'].parameters[name] = value
            return
      #if name == 'Rv_host'
      self.__dict__[name] = value

   def choose_model(self, name, stype='dm15'):
      '''A convenience function for selecting a model from the model module.
      [name] is the model to use.  The model will be used when self.fit() is 
      called and will contain all the parameters and errors.'''
      models = []
      for item in model.__dict__:
         obj = model.__dict__[item]
         if type(obj) is types.ClassType:
            if issubclass(obj, model.model):
               models.append(item)
      if name not in models:
         st = "Not a valid model.  Choose one of:  "+str(models)
         raise ValueError, st

      self.model = model.__dict__[name](self, stype=stype)
      self.template_bands = [b for b in self.model.rbs \
            if b not in ['Bs','Vs','Rs','Is']]
     
   def get_mag_table(self, bands=None, dt=0.5, outfile=None):
      '''This routine returns a table of the photometry, where the data from
      different filters are grouped according to day of observation.  The
      desired filters can be specified in [bands], otherwise all filters are
      returned.  The paramter [dt] controls how to group by time:  observations
      separated by less than [dt] in time are grouped.  To have the output
      sent to a file, specify a filename for [outfile].  When data is missing, a
      value of 99.9 is inserted.  The data is returned as a dictionary with the
      following keys: 
      'MJD':  the epoch of observation
      [band] :  The magnitude in filter [band]
      e_[band]:  the error in [band].
      '''

      if bands is None:  bands = self.data.keys()

      ret_data = {}
      # First, we make a list of observation dates from all the bands.
      times = [self.data[band].MJD for band in bands]
      times = sort(concatenate(times))

      # Eliminate repeating days:
      gids = concatenate([[1], greater(absolute(times[0:-1] - times[1:]), dt)])
      times = compress(gids, times)

      ret_data['MJD'] = times
      # Now loop through the bands and see where we need to fill in data
      for band in bands:
         gids = less(absolute(times[:,NewAxis] - \
               self.data[band].MJD[NewAxis,:]), dt)
         temp1 = 0.0*times + 99.9
         temp2 = 0.0*times + 99.9
         for i in range(len(gids)):
            if sum(gids[i]) > 0:
               temp1[i] = sum(self.data[band].mag*gids[i])/sum(gids[i])
               temp2[i] = max(sqrt(sum(power(self.data[band].e_mag,2)*gids[i]))/sum(gids[i]),
                              sqrt(average(power(temp1[i] - \
                                   compress(gids[i], self.data[band].mag),2))))
         ret_data[band] = temp1
         ret_data["e_"+band] = temp2

      if outfile is not None:
         if type(outfile) is types.StringType:
            fp = open(outfile, 'w')
         elif type(outfile) is types.FileType:
            fp = outfile
         else:
            raise TypeError, "outfile must be a file name or file handle"
         JDlen = len(str(int(ret_data['MJD']))) + 3
         title = "MJD" + " "*(JDlen+2)
         for b in bands:  title += "%5s +/-   " % b
         print >> fp, title
         format = "%%%d.2f  " + "%5.2f %4.2f  "*len(bands)
         for i in range(len(ret_data['MJD'])):
            data = []
            for b in bands:  data += [ret_data[b][i], ret_data['e_'+b][i]]
            print >> fp, format % tuple(data)
         fp.close()
         return
      else:
         return(ret_data)

   def lira(self, Bband, Vband, interpolate=0, tmin=30, tmax=90, plot=0):
      '''Use the Lira Law to derive a color excess.  [Bband] and [Vband] 
      should be
      whichever observed bands corresponds to restframe B and V, respectively.
      Use [interpolate]=1 to interpolate missing data.  If [interpolate]=0, then
      no color is computed where data is missing.  The color excess is
      estimated to be the median of the offset between the Lira line and the 
      data.  The uncertainty is 1.49 times the median absolute deviation of
      the offset data from the Lira line.  If you want to restrict the data
      used, use [tmin] and [tmax] to define a window.  If you want a graph use
      [plot]=1.  Returns a 3-tuple:  the E(B-V), error, and the fit splope (which
      can be used as a diagnostic).'''

      # find V-maximum
      t_maxes,maxes,e_maxes,restbands = self.get_rest_max([Vband])
      Tmax = t_maxes[0]

      t,BV,eBV,flag = self.get_color(Bband, Vband, kcorr=1)

      # find all points that have data in both bands
      gids = equal(flag, 0)
      # If we're allowed to interpolate, add flag=1
      if interpolate:
         gids = gids + equal(flag, 1)

      # Now apply a time criterion
      gids = gids*greater_equal(t-Tmax, tmin)*less_equal(t-Tmax, tmax)

      # Now check that we actually HAVE some data left
      if len(nonzero(gids)) == 0:
         raise RuntimeError, "Sorry, no data available between t=%f and t=%f" % (tmin,tmax) 
      
      # extract the data we want and convert to Vmax epochs
      t2 = compress(gids, (t-Tmax)/(1+self.z))
      BV2 = compress(gids, BV)
      eBV2 = compress(gids, eBV)
      
      # Next, solve for a linear fit (as diagnostic)
      w = power(eBV2,-2)
      c,ec = fit_poly.fitpoly(t2, BV2, w=w, k=1, x0=55.0)
      rchisq = sum(power(BV2 - c[0] - c[1]*(t2-55.),2)*w)/(len(BV2) - 2)
      ec = ec*sqrt(rchisq)

      lira_BV = 0.732 - 0.0095*(t2 - 55.0)
      #lira_EBV = stats.median(BV2 - lira_BV)
      #e_lira_EBV = 1.49*stats.median(absolute(BV2 - lira_BV - lira_EBV))
      w = power(eBV2,-2)
      lira_EBV = sum((BV2 - lira_BV)*w)/sum(w)
      e_lira_EBV = power(sum(w), -0.5)

      print "Vmax occurred at %f" % (t_maxes[0])
      print "Slope of (B-V) vs. t-Tvmax was %f(%f)" % (c[1], ec[1])
      print "median E(B-V) = %f    1.49*mad(E(B-V)) = %f" % (lira_EBV, e_lira_EBV)
      if absolute(c[1] + 0.0118) > 3*ec[1]:
         print "WARNING:  fit slope differs from Lira Law by more than three sigma"

      if plot:
         plotmod.plot_lira(t, t2, t_maxes, BV, eBV, BV2, tmin, tmax, c)

      return (lira_EBV, e_lira_EBV, c[1], ec[1])

   def get_rest_max(self, bands, deredden=0):
      return self.get_max(bands, deredden=deredden, restframe=1)

   def get_max(self, bands, restframe=0, deredden=0):
      '''Get the rest-frame maximum magnitue in [bands] based on the currently
      defined model or spline fits.  If you want rest-frame maxima 
      (i.e., have the k-corrections at maximum subtracted), set [restframe]=1.  
      If you want the reddening (galactic and/or host) removed, set 
      [deredden]=1.  If both a model and spline fit are defined for a filter,
      the spline will be taken.  The function returns:
      (Tmax,Mmax,e_Mmax,rband)
      Tmax = array of times of maximum, Mmax = array of maximum magnitudes,
      e_Mamx = error in max magnitudes, rband=rest-band for each filter.'''
      model_bands = [b for b in bands if b in self.model._fbands]
      lc_model_bands = [b for b in bands if self.data[b].Mmax is not None]
      for band in bands:
         if band not in model_bands and band not in lc_model_bands:
            raise ValueError, "Error:  filter %s has not been fit " % band + \
                  "with a light-curve yet, so I cannot compute it's maximum"
      N = len(bands)
      result = (zeros(N, dtype=Float32), zeros(N, dtype=Float32),
            zeros(N, dtype=Float32), [""]*N)
      if len(model_bands) > 0:
         mod_result = self.model.get_max(model_bands, restframe=restframe,
               deredden=deredden)
      for i in range(N):
         b = bands[i]
         if b in model_bands:
            mid = model_bands.index(b)
            for j in range(4): result[j][i] = mod_result[j][mid]
         if b in lc_model_bands:
            result[0][i] = self.data[b].Tmax
            result[1][i] = self.data[b].Mmax
            result[2][i] = self.data[b].e_Mmax
            result[3][i] = b
      return result

   def kcorr(self, bands=None, mbands=None, mangle=1, interp=1, use_model=0, 
         min_filter_sep=400, use_stretch=1, **mopts):
      '''Compute the k-corrections for the filters [bands] (default:  do all
      filters in self.data).  In order to get the best k-corrections possible,
      we warp the SNIa SED (defined by self.k_version) to match the observed
      photometry defined by the filters in [mbands] (default:  same as [bands])
      unless [mangle]=0.  Not all bands will be observed on the same day (or some
      data may be less than reliable), so there are several arguments that
      control how the warping is done.  First, only filters whose effective
      wavelengths are separated by more than [min_filter_sep] will be used (to
      avoid unstable splines).  For days with no data, colors will be
      constructed by GLOEs interpolation, unless [interp]=0, in which case only
      colors that are observed will be used.  If you would rather use the
      light-curve model currently defined, set [use_model]=1 (all points are
      interpolated based on the model).  If [use_stretch]=1, then the value
      of dm15 or st is mapped to a time stretch and applied to the SED templates
      to take into account that faster decliners evolve more quickly.
      Lastly, you can use any options used by mangle_spectrum.mangle_spectrum2
      (see it's docstring).

      Upon successful completion, the following member variables will be
      populated:
        self.ks:      dictionary (indexed by filter) of k-corrections
        self.ks_mask: dictionary indicating valid k-corrections
        self.ks_tck:  dictionary of spline coefficients for the k-corrections
                      (useful for interpolating the k-corrections).
        self.mopts:   If mangling was used, contains the parameters of the mangling
                      function.
      '''
      if use_stretch and self.k_version != '91bg':
         dm15 = getattr(self, 'dm15', None)
         st = getattr(self, 'st', None)
         if dm15 is None and st is None:
            raise AttributeError, "Before you can k-correct with stretch, you"+\
                  " need to solve for dm15 or st, using either a model or LC fit"
         if dm15 is None:
            s = st
         else:
            if dm15 > 1.7:
               print "Warning:  dm15 > 1.7.  Using this stretch on the Hsiao SED"
               print "  is not recommended.  I'm setting stretch to 1.0.  You"
               print "  might consider using the 91bg SED template."
               s = kcorr.dm152s(1.7)
            elif dm15 < 0.7:
               s = kcorr.dm152s(0.7)
            else:
               s = kcorr.dm152s(dm15)
      elif use_stretch and self.k_version == '91bg':
         print "Warning:  you asked for stretching the template SED, but"
         print "you have selected the 91bg template.  Setting stretch to 1.0."
         s = 1.0
      else:
         s = 1.0
      self.ks_s = s

      if bands is None:  bands = self.data.keys()
      if mbands is None:  mbands = [b for b in bands]
      # Check the simple case:
      if not mangle:
         for band in bands:
            x = self.data[band].MJD
            # days since Bmax in the frame of the SN
            days = (x - self.Tmax)/(1+self.z)/s
            days = days.tolist()
            self.ks[band],self.ks_mask[band] = map(array,kcorr.kcorr(days, 
               self.restbands[band], band, self.z, self.EBVgal, 0.0,
               version=self.k_version))
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
      while sometrue(less(dwaves, min_filter_sep)):
         bids = less(dwaves, min_filter_sep)
         mbands = [mbands[i] for i in range(len(bids)) if not bids[i]] + mbands[-1:]
         eff_waves = array([fset[band].eff_wave(Ia_w,Ia_f) for band in mbands])
         dwaves = eff_waves[1:] - eff_waves[0:-1]
      if not self.quiet:  print "Mangling based on filters:", mbands

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
         raise AttributeError, \
               "Error.  self.Tmax must be set in oder to compute K-correctsions"
      t = res['MJD'] - self.Tmax
      if not sometrue(greater_equal(t, -19)*less(t, 70)):
         raise RuntimeError, \
            "Error:  your epochs are all outside -20 < t < 70.  Check self.Tmax"
      kcorrs,mask,Rts,m_opts = kcorr.kcorr_mangle2(t/(1+self.z)/s, bands, 
            mags, masks, restbands, self.z, normfilter=mbands[-1], 
            colorfilts=mbands, version=self.k_version, full_output=1, **mopts)
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


   def get_mangled_SED(self, band, i):
      '''After the mangle_kcorr function has been run, you can use this function to
      retrieve the mangled SED that was used to compute the k-correction for the
      [i]'th day in [band]'s light-curve.  Returns 4 arrays:  wavelength, mangled_flux,
      original flux, and mangling function.'''
      
      if 'ks_mopts' not in self.__dict__:
         raise AttributeError, "Mangling info not found... try running self.kcorr()"
      epoch = self.data[band].t[i]/(1+self.z)/self.ks_s
      wave,flux = kcorr.get_SED(int(epoch), version=self.k_version)
      man_flux = mangle_spectrum.apply_mangle(wave,flux, **self.ks_mopts[band][i])
      return(wave,man_flux,flux,man_flux/flux)
   
   def get_color(self, band1, band2, interp=1, use_model=0, model_float=0, kcorr=0):
      '''return the observed SN color of [band1] - [band2].  If [interp]=1, then
      on days when only one band is measured, the other is interpolated,
      otherwise, only days when both bands are measured will be returned.
      If [use_model]=1, the fit model will be used, otherwise, the light-curve
      will be interpolated using a spline solution (if it exists)
      Set [kcorr]=1 if you want the results k-corrected.
      Returns a 4-tuple:
      (MJD, band1-band2, e_band1-band2, flag).  
       -  Flag is one of:
          0 - both bands measured at given epoch
          1 - only one band measured, other interpolated
          2 - extrapolation (based on template) needed, so reasonably safe
          3 - data interpolated or extrapolated beyond template or spline, not
              safe to use!'''

      if kcorr:
         if band1 not in self.ks_tck or band2 not in self.ks_tck:
            raise RuntimeError, "No k-corrections defined.  Either set " + \
                  "kcorr=0 or run self.kcorr() first"
      # First, get a table of all photometry:
      data = self.get_mag_table([band1,band2])

      if not interp:
         gids = less(data[band1], 90)*less(data[band2], 90)
         mjd = compress(gids, data['MJD'])
         col = compress(gids, data[band1]) - compress(gids, data[band2])
         if kcorr:
            k1 = scipy.interpolate.splev(mjd, self.ks_tck[band1])
            k2 = scipy.interpolate.splev(mjd, self.ks_tck[band2])
            col = col - k1 + k2
         ecol = sqrt(power(data['e_'+band1][gids], 2) + 
                     power(data['e_'+band2][gids], 2))
         return(mjd, col, ecol, floor(0*mjd).astype('l'))

      # Now, if we are interpolating, do so by fitting each band 
      # independently: 
      # Just weighted average between model and data
      interps = []
      masks = []
      doffsets = []
      offsets = []
      if use_model:
         for band in [band1,band2]:
            if float_model:
               temp,etemp,mask = self.model(band, self.data[band].MJD)
               weight = self.data[band].e_flux**2
               weight = power(weight, -1)*mask*self.data[band].mask
               offsets.append(sum((self.data[band].mag - temp)*weight)/sum(weight))
               doffsets.append(stats.median(absolute(self.data[band].mag - \
                                                     temp-offsets[-1])))
            else:
               offsets.append(0)
               doffsets.append(0)
            temp,etemp,mask = self.model(band, data['MJD'])
            interps.append(temp + offsets[-1])
            masks.append(mask)
      else:
         for band in [band1,band2]:
            temp,mask = self.data[band].eval(data['MJD'])
            interps.append(temp)
            masks.append(mask)
            doffsets.append(0)

      # Now we have interpolated values for band1,band2 where needed
      # First, where both bands are measured, flag as 0, otherwise only one 
      #  is measured and we flag as 1
      m1 = less(data[band1], 90);  m2 = less(data[band2], 90)
      flag = where(m1*m2, 0, 1)

      # Get the range where we are doing interpolation
      i1 = max(min(nonzero(m1)), min(nonzero(m2)))
      i2 = min(max(nonzero(m1)), max(nonzero(m2)))
      # from 0 to i1 (non-inclusive) and i2 to end, we have extrapolation, flag
      # as 2
      flag[0:i1] = 2
      flag[i2+1:] = 2
      # flag places where we don't have data and template/spline is not valid
      flag = where(greater(flag,0)*-(masks[0]*masks[1]),3,flag)

      # Now that the flags are set properly, we do the math:
      b1 = where(less(data[band1], 90), data[band1], interps[0])
      e_b1 = where(less(data["e_"+band1], 90), data["e_"+band1], doffsets[0])
      b2 = where(less(data[band2], 90), data[band2], interps[1])
      e_b2 = where(less(data["e_"+band2], 90), data["e_"+band2], doffsets[1])
      colors = b1 - b2
      if kcorr:
         k1 = scipy.interpolate.splev(data['MJD'], self.ks_tck[band1])
         k2 = scipy.interpolate.splev(data['MJD'], self.ks_tck[band2])
         colors = colors - k1 + k2
      e_colors = sqrt(power(e_b1, 2) + power(e_b2, 2))

      return(data['MJD'], colors, e_colors, flag)

   def getEBVgal(self):
      '''Gets the value of E(B-V) due to galactic extinction.  The ra and decl
      member varialbles must be set beforehand.'''
      if self.ra is not None and self.decl is not None:
         self.EBVgal,mask = dust_getval.get_dust_RADEC(self.ra, self.decl)
         self.EBVgal = self.EBVgal[0]
      else:
         print "Error:  need ra and dec to be defined, E(B-V)_gal not computed"

   def summary(self, out=sys.stdout):
      '''Get a quick summary of the data for this SN, along with fitted 
      parameters (if such exist).'''
      print >> out, '-'*80
      print >> out, "SN ",self.name
      if self.z:  print >> out, "z = %.3f         " % (self.z),
      if self.ra:  print >> out, "ra=%9.5f        " % (self.ra),
      if self.decl:  print >> out, "dec=%9.5f" % (self.decl),
      print >> out, ""
      print >> out, "Data in the following bands:",
      for band in self.data:  print >> out, band + ", ",
      print >> out, ""

      print >> out, "Fit results (if any):"
      for band in self.restbands:
         print >> out, "   Observed %s fit to restbad %s" % (band, self.restbands[band])
      for param in self.parameters:
         if self.parameters[param] is not None:
            print >> out, "   %s = %.3f  +/-  %.3f" % (param, self.parameters[param],
                                                       self.errors[param])

   def dump_lc(self, epoch=0, tmin=-10, tmax=70, k_correct=0):
      '''Outputs several files that contain the lc information, both the data,
      uncertainties, and the models themselves.  If epoch=1, remove Tmax from
      the times.  tmin and tmax are the beginning and end epochs for the model.
      If k_correct=1, then apply the k-corrections to the data and model before
      outputting to the file.  Otherwise, the original data will be output.
      Only times for which the model (templates) are valid will be output.  This
      function will create several files with the following template names:
         {SN}_lc_{filter}_data.dat
         {SN}_lc_{filter}_model.dat
         
      which will contain the photometric data and model for each filter (if
      that filter was fit with a model).  In the *_data.dat, there is an
      extra flag column that indicates if the k-corrections are valid (0) or
      invalid (1).'''
      base = self.name + "_lc_"
      if not epoch:
         toff = self.Tmax
      else:
         toff = 0
      for filter in self.data.keys():
         f = open(base+filter+"_data.dat", 'w')
         print >> f, "#  column 1:  time"
         print >> f, "#  column 2:  oberved magnitude"
         print >> f, "#  column 3:  error in observed magnitude"
         print >> f, "#  column 4:  Flag:  0=OK  1=Invalid K-correction"
         for i in range(len(self.data[filter].mag)):
            if k_correct:
               flag = 0
               if filter not in self.ks:
                  flag = 1
                  ks = 0
               else:
                  flag = (not self.ks_mask[filter][i])
                  ks = self.ks[filter][i]
               print >> f, "%.2f  %.3f  %.3f  %d" % \
                     (self.data[filter].t[i]+toff, 
                     self.data[filter].mag[i] - ks, 
                     self.data[filter].e_mag[i], flag)
            else:
               print >> f, "%.2f  %.3f  %.3f  %d" % \
                     (self.data[filter].t[i]+toff, 
                     self.data[filter].mag[i], self.data[filter].e_mag[i], 0)
         f.close()
         if filter in self.model._fbands:
            ts = arange(tmin, tmax+1, 1.0)
            ms,e_ms,mask = self.model(filter, ts+self.Tmax)
            if k_correct and filter in self.ks_tck:
               ks = scipy.interpolate.splev(ts + self.Tmax, self.ks_tck[filter])
               # mask out valid k-corrections
               mids = argmin(absolute(ts[:,NewAxis]-self.data[filter].MJD[NewAxis,:]+\
                     self.Tmax))
               ks_mask = self.ks_mask[filter][mids]*greater_equal(ts, -19)*less_equal(ts, 70)
               mask = mask*ks_mask
               ms = ms - ks
            ms = ms[mask]
            ts = ts[mask]
            f = open(base+filter+"_model.dat", 'w')
            print >>f, "# column 1: time"
            print >>f, "# column 2:  model magnitude"
            for i in range(len(ts)):
               print >> f, "%.1f, %.3f" % (ts[i]+toff, ms[i])
            f.close()
         if self.data[filter].tck is not None:
            f = open(base+filter+"_spline.dat", 'w')
            ts = arange(self.data[filter].tck[0][0], 
                  self.data[filter].tck[0][-1]+1, 1.0)
            m,mask = self.data[filter].eval(ts, t_tol=-1)
            print >> f, "# column 1:  time"
            print >> f, "# column 2:  splined magnitude"
            for i in range(len(ts)):
               if not mask[i]:  continue
               print >> f, "%.1f  %.3f" % (ts[i]+toff-self.Tmax, m[i])
            f.close()

   def update_sql(self, attributes=None, dokcorr=1):
      '''Updates the current information in the SQL database, creating a new SN
      if needed.   If attributes are specified (as a list of strings), then only
      these attributes are updated.'''
      if have_sql:
         N = self.sql.connect(self.name)
         if N == 0:
            self.sql.create_SN(self.ra, self.decl, self.z)
            data = {}
            for f in self.data:
               if dokcorr:
                  data[f] = [self.data[f].JD, self.data[f].magnitude, self.data[f].e_mag,
                          self.data[f].K]
               else:
                  data[f] = [self.data[f].JD, self.data[f].magnitude, self.data[f].e_mag,
                          None]
            self.sql.create_SN_photometry(data)
         elif dokcorr:
            for f in self.data:
               self.sql.update_photometry(f, self.data[f].MJD, "K", self.data[f].K)
         attr_list = ['z','ra','decl'] + self.parameters.keys()
         for attr in attr_list:
            try:
               self.sql.set_SN_parameter(attr, self.__getattr__[attr])
            except:
               pass
         self.sql.close()

   def read_sql(self, name):
      '''Get the data from the SQL server for supernova 'name'.'''
      if have_sql:
         N = self.sql.connect(name)
         if N == 0:
            print "%s not found in database, starting from scratch..." % (name)
            self.sql.close()
            return(-1)
         try:
            self.z = self.sql.get_SN_parameter('z')
            self.ra = self.sql.get_SN_parameter('ra')
            self.decl = self.sql.get_SN_parameter('decl')
            for param in self.parameters:
               try:
                  self.parameters[param] = self.sql.get_SN_parameter(param)
                  self.errors[param] = self.sql.get_SN_parameter('e_'+param)
               except:
                  pass
            data = self.sql.get_SN_photometry()
            for filter in data:
               d = data[filter]
               #if d[3] is not None:
               #   self.data[filter] = lc(self, filter, d[0], d[1], d[2], K=d[3])
               #else:
               self.data[filter] = lc(self, filter, d[0], d[1], d[2])
         finally:
            self.sql.close()

   def get_restbands(self):
      '''Automatically populates the restbands member data with one of 'B','V',
      'R','I', whichever's effective wavelength is closest to the observed 
      bands.'''
      for band in self.data:
         self.restbands[band] = self.closest_band(band)

   def lc_offsets(self, min_off=0.5):
      '''Find offsets such that the lcs, when plotted, won't overlap.  Specify
      [filters] the order you want them, otherwise, they are chosen in
      order of increasing effective wavelength.'''

      if self.filter_order is None:
         bands = self.data.keys()
         eff_wavs = []
         for filter in bands:
            eff_wavs.append(fset[filter].ave_wave)
         eff_wavs = asarray(eff_wavs)
         ids = argsort(eff_wavs)
         self.filter_order = [bands[i] for i in ids]

      offs = [0]
      filter = self.filter_order[0]

      f = interp1d(self.data[filter].MJD, self.data[filter].mag, bounds_error=False,
            fill_value=self.data[filter].mag.min())
      for filter in self.filter_order[1:]:
         off = offs[-1] + min_off
         while not alltrue(greater(self.data[filter].mag+off - f(self.data[filter].MJD),
            0.5)):
            off += 0.5
         offs.append(off)
         if self.data[filter].MJD.shape[0] > 1:
            f = interp1d(self.data[filter].MJD, self.data[filter].mag+off, bounds_error=False,
                           fill_value=self.data[filter].mag.max())
         else:
            f = lambda x:  self.data[filter].mag+off
      return offs



   def save(self, filename):
      '''Save this SN instance to a pickle file, which can be loaded again
      using the get_sn() function.'''
      f = open(filename, 'w')
      pickle.dump(self, f)
      f.close()
   
   def fit(self, bands=None, mangle=1, kcorr=1, reset_kcorrs=1, k_stretch=True, 
         margs={}, **args):
      '''Fit the N light curves with the currently set model (see 
      self.choose_model()).  The parameters that can be varried or held 
      fixed depending on the model being used (try help(self.model)
      for this info).  If one of these parameters is specified with a 
      value as an argument, it is held fixed.  Otherwise it is varied.  
      If you set a parameter to None, it will be automoatically chosen by 
      self.model.guess().

      There are several optional arguments that change how the fit deals
      with k-corrections:
         - kcorr=1:   If non-zero, compute k-corrections and save them to 
                      self.ks before doing the fit.  These k-corrections 
                      will be used as part of the model and will replace 
                      any previously computed k-corrections.
         - mangle=1:  if non-zero, compupte the k-corrections based on the
                      Hsiao SED mangled to fit the model of the photometry.
         - reset_kcorrs=1:  before the first fit, zero any previous 
                      k-corrections.  Setting  this to false and kcorr to 
                      false will use the k-corrections previously computed 
                      (i.e., those saved in s.ks).
      NOTE:  If you have data that has already been k-corrected (either
             outside SNooPy or by setting the individual data's K
             attributes, use kcorr=0 and reset_kcorrs=1.  If you have
             run the self.kcorr() manually and want to keep those
             k-corrections, use kcorr=0 and reset_kcorrs=0.  Otherwise,
             use the default kcorr=1 and reset_kcorrs=1.'''

      if bands is None:
         # By default, we fit the bands whose restbands are provided by the model
         bands = [b for b in self.data.keys() if self.restbands[b] in self.model.rbs]
      # Setup initial Robs (in case it is used by the model)
      for band in bands:
         if band not in self.Robs:
            self.Robs[band] = fset[band].R(self.Rv_gal, Ia_w, Ia_f, z=self.z)

      if self.z <= 0:
         raise ValueError, "The heliocentric redshift is zero.  Fix this before you fit"

      # Check to make sure we have filters we can fit:
      for filter in bands:
         if self.restbands[filter] not in self.model.rbs:
            raise AttributeError, \
                  "Error:  filter %s is not supported by this model" % filter+\
                  ", set self.restbands accordingly"

      if reset_kcorrs:
         self.ks = {}
         self.ks_mask = {}
         self.ks_tck = {}
      if not self.quiet:
         print "Doing Initial Fit to get Tmax..."
      self.model.fit(bands, **args)


      if kcorr:
         kbands = [band for band in bands if band not in self.ks]
         if len(kbands) > 0:
            if not self.quiet:
               print "Setting up initial k-corrections"
            self.kcorr(kbands, mangle=0)
 
         if not self.quiet:
            if mangle:
               print "Doing first fit..."
            else:
               print "Doing fit..."
         self.model.fit(bands, **args)
 
         if mangle:
            if not self.quiet:
               print "Doing mangled k-corrections"
            self.kcorr(bands, interp=0, use_model=1, use_stretch=k_stretch, **margs)
            if not self.quiet:
               print "Doing final fit..."
            self.model.fit(bands, **args)
      if self.replot:
         self.plot()

   def systematics(self, **args):
      '''Returns a dictionary of systematics errors keyed by paramter
      name.  It therefore depends on the model being used.  Also
      see the specific model for any extra arguments.  If None
      is returned as a value, no systematic has been estimated
      for it.'''
      return self.model.systematics(**args)

   def plot_filters(self, bands=None, day=0, **args):
      return plotmod.plot_filters(self, bands, day, **args)

   def compute_w(self, band1, band2, band3, R=None):
      '''Returns the reddeining-free magnitude in the sense that:
      w = band1 - R(band1,band2,band3)*(band2 - band3)
      for for instance compute_w(V,B,V) would give:
      w = V - Rv(B-V)'''
      # First, let's get the proper value of R:
      if R is None:
         R1 = fset[band1].R(self.Rvhost, Ia_w, Ia_f)
         R2 = fset[band2].R(self.Rvhost, Ia_w, Ia_f)
         R3 = fset[band3].R(self.Rvhost, Ia_w, Ia_f)
         R = R1/(R2 - R3)

      # Now, we're probably going to have to interpolate band2 and band3 to get
      # the colors at times of band1, so let's spline it.
      t = self.data[band1].MJD - self.Tmax
      m = self.data[band1].mag
      e_m = self.data[band1].e_mag
      t2 = self.data[band2].MJD - self.Tmax
      t3 = self.data[band3].MJD - self.Tmax
      m2 = self.data[band2].mag
      e_m2 = self.data[band2].e_mag
      m3 = self.data[band3].mag
      e_m3 = self.data[band3].e_mag

      ev2 = fit_spline.interp_spline(t2, m2, e_m2, t, k=3)
      ev3 = fit_spline.interp_spline(t3, m3, e_m3, t, k=3)

      # Now compute w:
      w = m - R*(ev2 - ev3)
      return(w)

   def mask_data(self):
      '''Interactively mask out bad data and unmask the data as well.  The only
      two bindings are "A" (click):  mask the data and "u" to unmask the data.
      '''
      return plotmod.mask_data(self)
   
   def plot(self, xrange=None, yrange=None,  
         title=None, interactive=0, single=0, dm=1, fsize=None, linewidth=None,
         symbols=None, colors=None, relative=0, legend=1, mask=1, label_bad=0,
         flux=0, epoch=1, msize=6, outfile=None):
      '''Plot out the supernova data in a nice format.  There are several 
      options:
         - xrange,yrange:  specify the ranges to plot as lists [xmin,xmax], 
           [ymin,ymax]
         - title:  optional title
         - interactive:  allows for an interactive plot (PGPLOT only).
         - single:  plot out as a single (rather than panelled) plot?
         - dm:  offset in magnitudes between the lightcurves (for single plots)
         - fsize:  override the font size used to plot the graphs
         - linewidth:  override the line width
         - symbols:  dictionary of symbols, indexed by band name.
         - colors:  dictionary of colors to use, indexed by band name.
         - relative:  plot only relative magnitudes (normalized to zero)?
         - legend:  do we plot the legend?
         - mask:  Omit plotting masked out data?
         - label_bad:  label the masked data with red x's?
         - flux:  (boolean) plot in flux units?
         - epoch:  (boolean) plot time relative to Tmax?
         - outfile:  if supplied, save the plot to [outfile]
      '''

      return plotmod.plot_sn(self, xrange, yrange,
         title, interactive, single, dm, fsize, linewidth,
         symbols, colors, relative, legend, mask, label_bad,
         flux, epoch, msize, outfile)

   def plot_kcorrs(self, colors=None, symbols=None, outfile=None):
      '''Plot the derived k-corrections after they have been computed.
      Both mangled and un-mangled k-corrections will be plotted as
      lines and points, respectively.  If mangling was used to
      do the k-corrections, clicking 'm' on a point will bring up
      another plot showing the original and mangled spectrum. You can
      specify colors and symbols with [colors] and [symbols].  If
      [outfile] is proviced, output the graph to [outfile]'''
      return plotmod.plot_kcorrs(self, colors, symbols)

   def bolometric(self, bands, lam1=None, lam2=None, refband=None, 
         normband=None, remangle=0, extrap_red='RJ', extrap_blue=None, 
         outfile=None, verbose=0, **mopts):
      '''Produce a quasi-bolometric flux light-curve based on the input [bands]
      by integrating a template SED from \lambda=lam1 to \lambda=lam2.
      The bands are used to mangle a SNIa SED template (determined by
      self.k_version).  If mangled k-corrections have already been
      performed (self.kcorr(mangle=1)), then those mangling functions
      will be used here, unless remangle=1.   refband is used to determine
      the cadence (so should be the band with the *least* coverage.
      Photometry will be interpolated to this cadence.  The bolometric
      flux will be normalized to match the data in normband (default is
      to use the first filter in bands).  To extrapolate the red end of
      the spectrum using Rayleigh-Jeans, use extrap_red='RJ', otherwise
      set it to None.  To specify a function of wavelength for the blue
      end, use [extrap_blue].  Specifying [outfile] will save the graph.'''
      for b in bands:  
         if b not in self.data:
            raise AttributeError, "band %s not defined in data set" % (b)

      if lam1 is None:
         lam1 = array([fset[b].waverange()[0] for b in bands]).min()
      if lam2 is None:
         lam2 = array([fset[b].waverange()[1] for b in bands]).max()
      if verbose:
         print "Integrating from %.1f to %.1f" % (lam1, lam2)
      if refband is None:
         # take the band with the fewest data points
         nps = array([len(self.data[b].MJD) for b in bands])
         id = argmin(nps)
         refband = bands[id]
      if verbose:  print "Using %s as the reference band" % refband

      if refband not in bands:
         raise ValueError, "refband must be one of bands"

      if normband is None:
         normband = bands[0]
      if verbose:  print "Normalizing to flux in %s band" % normband

      if 'ks_mopt' not in self.__dict__:
         ks_mopt = {}
      else:
         ks_mopt = self.ks_mopt
      if not alltrue(array([b in ks_mopt for b in bands])):
         remangle = 1

      if remangle:
         self.kcorr(bands, mangle=1, **mopts)

      MJD = []
      bol = []
      bol_mask = []
      # Now we compute the bolometric lightcurve (in magnitudes)
      #  by doing a cross-band k-correction to the box filter
      for i in range(len(self.data[refband].t)):
         MJD.append(self.data[refband].MJD[i])
         if not self.ks_mask[refband][i]:
            bol.append(0.0)
            bol_mask.append(False)
            continue
         wave,mflux,flux,ratio = self.get_mangled_SED(refband, i)
         if (wave[0] > lam1 and not extrap_blue) or \
               (wave[-1] < lam2 and not extrap_red):
            raise ValueError, "The template SED does not cover the requested range."
         l1 = lam1;   l2 = lam2

         # normalize to observed flux
         mobs,mask = self.data[normband].eval(MJD[-1], t_tol=-1)
         if not mask:
            bol.append(0.0)
            bol_mask.append(False)
         fobs = power(10, -0.4*(mobs - fset[normband].zp))
         fint = fset[normband].response(wave*(1+self.z), mflux)
         mflux = mflux*fobs/fint

         # Make sure we are integrating within the SED
         if wave[0] > lam1:  
            l1 = wave[0]
         if wave[-1] < lam2:
            l2 = wave[-1]
         i_min = argmin(absolute(wave - l1))
         i_max = argmin(absolute(wave - l2))

         # integrate the array from l1 to l2, adding (or subtracting) any
         #   bits by interpolation.
         flx = trapz(mflux[i_min:i_max+1], x=wave[i_min:i_max+1])
         if wave[i_min] > l1:
            flx += trapz(mflux[i_min-1:i_min+1], x=wave[i_min-1:i_min+1])*\
                  (wave[i_min]-l1)/(wave[i_min]-wave[i_min-1])
         else:
            flx -= trapz(mflux[i_min:i_min+2], x=wave[i_min:i_min+2])*\
                  (l1 - wave[i_min])/(wave[i_min+1]-wave[i_min])
         if wave[i_max] > l2:
            flx -= trapz(mflux[i_max-1:i_max+1], x=wave[i_max-1:i_max+1])*\
                  (wave[i_max]-l2)/(wave[i_max]-wave[i_max-1])
         else:
            flx += trapz(mflux[i_max:i_max+2], x=wave[i_max:i_max+2])*\
                  (l2 - wave[i_max])/(wave[i_max+1]-wave[i_max])

         # Now we need to handle extrapolation
         if lam2 > wave[-1] and extrap_red is not None:
            if extra_red == 'RJ':
               # Do a Rayleigh-Jeans extrapolation (~ 1/lam^4)
               # normalize to wave[-1]
               flx += mflux[-1]/3*lam2
            else:
               raise ValueError, "Unrecognized red extrapolation method"
         if lam1 < wave[0] and extrap_blue is not None:
            flx += extrap_blue

         bol_mask.append(True)
         bol.append(flx)
      return(self.dadta[refband].t, array(bol),array(bol_mask))

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

def load(file):
   try:
      f = open(file, 'r')
      inst = pickle.load(f)
      f.close()
   except:
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
      header = string.join(header, separator)
      header[0] = "#"
      print >>f, header

   for i in range(len(arrays[0])):
      line = [formats[j] % (arrays[j][i]) for j in range(len(arrays))]
      print >>f, string.join(line, separator)
   f.close()

def fix_arrays(node):
   '''A recursive function that seeks out Numeric arrays and replaces them
   with numpy arrays.'''
   from Numeric import ArrayType
   if type(node) is ArrayType:
      return array(node)
   elif type(node) is types.InstanceType:
      for key in node.__dict__:
         if key != 'parent':
            node.__dict__[key] = fix_arrays(node.__dict__[key])
      return node
   elif type(node) is types.DictType:
      for key in node.keys():
         if key != 'parent':
            node[key] = fix_arrays(node[key])
      return node
   elif type(node) is types.ListType:
      return [fix_arrays(item) for item in node]
   elif type(node) is types.TupleType:
      return tuple([fix_arrays(item) for item in node])
   else:
      return node

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
   f = open(file)
   lines = f.readlines()
   fields = lines[0].split()
   if len(fields) != 4:  raise RuntimeError, "first line of %s must have 4 " +\
         "fields:  name, redshift, RA, DEC"
   name = fields[0]
   try:
      z,ra,decl = map(float, fields[1:])
   except:
      raise RuntimeError, "z, ra and dec must be floats " + \
            " (ra/dec in decimal degrees)"


   s = sn(name, ra=ra, dec=decl, z=z)
   lines = lines[1:]
   this_filter = None
   MJD = {}
   mags = {}
   emags = {}

   for line in lines:
      if line[0] == "#":  continue
      if line.find('filter') >= 0:
         this_filter = line.split()[1]
         MJD[this_filter] = []
         mags[this_filter] = []
         emags[this_filter] = []
      elif this_filter is not None:
         try:
            t,m,em = map(float, string.split(string.strip(line)))
         except:
            raise RuntimeError, "Bad format in line:\n %s" % (line)
         MJD[this_filter].append(t)
         mags[this_filter].append(m)
         emags[this_filter].append(em)

   for f in MJD:
      MJD[f] = array(MJD[f])
      mags[f] = array(mags[f])
      emags[f] = array(emags[f])
      s.data[f] = lc(s, f, MJD[f], mags[f], emags[f])
      s.data[f].time_sort()

   s.get_restbands()

   return(s)

def get_sn(str, sql=None):
   '''Attempt to get a sn object from several possible sources.  First, if str
   corresponds to an existing file name, the function attempts to load the sn
   instance from it as if it were a pickle'd object.  If that fails, it attempts
   to use import_lc() on the file.  If str is not the name of an existing file,
   it is treated as a SN name and is retrieved from the designated sql connection
   object (or default_sql if sql=None).'''
   if os.path.isfile(str):
      try:
         f = open(str, 'r')
         s = pickle.load(f)
         return s
      except:
         try:
            s = import_lc(str)
         except RuntimeError:
            raise RuntimeError, "Could now load %s into SNPY" % str
   else:
      s = sn(str, source=sql)
   return s
