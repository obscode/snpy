#!/usr/bin/env python
'''
A module for computing k-corrections. This code has been assembled from many
sources, mostly from IRAF scripts written by Mark Phillips and IDL code 
written by Mark Sullivan.

Most of the heavy lifting w.r.t. integration of filters over the SED has
been moved to the :mod:`snpy.filter` module.
'''

import os,sys,string,re
import numpy as num
import scipy.interpolate
from utils import deredden
try:
   from astropy.io import fits as pyfits
except ImportError:
   try:
      import pyfits
   except ImportError:
      sys.stderr.write('Error:  You need pyfits to run snpy.  You can get it\n')
      sys.stderr.write('        from:  http://www.stsci.edu/resources/'+\
                       'software_hardware/pyfits/\n')
      raise ImportError
import filters
from mangle_spectrum import mangle_spectrum2

base = os.path.dirname(globals()['__file__'])
spec_base = os.path.join(base,'typeIa')

debug=0
h = 6.626068e-27
c = 2.997925e18
ch = c * h

# This converts dm15 into a stretch that can be used to
#    stretch the Hsiao template.  This is a first correction
#    to "warp" the SED template
def dm152s(dm15):
   return 2.13 - 2.44*dm15 + 2.07*dm15**2 - 0.7*dm15**3

# Load all the SED templates:
# Hsiao's uberspectrum:
f = pyfits.open(os.path.join(spec_base, 'Hsiao_SED_V2.fits'))
h_sed = f[0].data
head = f[0].header
h_wav = head['CRVAL1'] + (num.arange(head['NAXIS1'],dtype=num.float32) - \
      head['CRPIX1'] + 1)*head['CDELT1']
f.close()
# Hsiao's new OPT+NIR uberspectrum
f = pyfits.open(os.path.join(spec_base, 'Hsiao_SED_V3.fits'))
h3_sed = f[0].data
head = f[0].header
h3_wav = head['CRVAL1']+(num.arange(head['NAXIS1'],dtype=num.float32) - \
      head['CRPIX1'] + 1)*head['CDELT1']
f.close()
# Nugent's uberspectrum:
f = pyfits.open(os.path.join(spec_base, "Nugent_SED.fits"))
n_sed = f[0].data
head = f[0].header
n_wav = head['CRVAL1']+(num.arange(head['NAXIS1'],dtype=num.float32) - \
      head['CRPIX1'] + 1)*head['CDELT1']
# Nugent's 91bg-like SED:
f = pyfits.open(os.path.join(spec_base, "Nugent_91bg_SED.fits"))
n91_sed = f[0].data
head = f[0].header
n91_wav = head['CRVAL1']+(num.arange(head['NAXIS1'],dtype=num.float32) - \
      head['CRPIX1'] + 1)*head['CDELT1']

def linterp(spec1, spec2, day1, day2, day):
   if day1 == day2:
      return spec1
   if day > day2 or day < day1:
      raise ValueError, "day must be in interval (day1,day2)"
   if abs(day1-day) < 1e-9:
      return spec1
   if abs(day2 - day) < 1e-9:
      return spec2
   spec = spec1 + (spec2 - spec1)/(day2 - day1)*(day - day1)
   return spec

SED_lims = {
      'H':(-19,70),
      'H3':(-19,70),
      'N':(-19,70),
      '91bg':(-13,100)}

def get_SED(day, version='H3', interpolate=True, extrapolate=False):
   '''Retrieve the SED for a SN for a particular epoch.
   
   Args:
      day (int or float): The integer day w.r.t. time of B-maximum
      version (str): The version of SED sequence to use:

         * 'H': Old Hsiao Ia SED (Hsiao, private communication)
         * 'H3': Hsiao+2007 Ia SED
         * 'N': Nugent+2002 Ia SED
         * '91bg': a SN1991bg Ia SED (Peter Nugent)
      interpolate(bool): If and day is not an integer, interpolate
                         the spectrum linearly. Otherwise, choose
                         nearest spectrum.
      extrapolate(bool): If True and the date is outside the range
                         of defined SED, simply take the first/last
                         SED to extend before/after range.

   Returns:
      2-tuple: (wave,flux):

      * wave (array):  Wavelength in Angstroms
      * flux (array):  arbitrarily normalized flux
   '''

   if type(day) is type(1.0) and not interpolate:
      day = round(day)
   
   # Check limits
   if day < SED_lims[version][0]:
      if extrapolate:
         day = SED_lims[version][0]
      else:
         return (None,None)
   if day > SED_lims[version][1]:
      if extrapolate:
         day = SED_lims[version][1]
      else:
         return (None,None)

   day1 = int(num.floor(day))
   day2 = int(num.ceil(day))
   if version == 'H':
      return (h_wav, 
            linterp(h_sed[day1+20,:],h_sed[day2+20,:],day1,day2,day))
   elif version == 'H3':
      return (h3_wav, 
            linterp(h3_sed[day1+20,:],h3_sed[day2+20,:],day1,day2,day))
   elif version == 'N':
      return (n_wav, 
            linterp(n_sed[day1+19,:],n_sed[day2+19,:],day1,day2,day))
   elif version == '91bg':
      return (n91_wav, 
            linterp(n91_sed[day1+13, :],n91_sed[day2+13,:],day1,day2,day))
   else:
      raise AttributeError, "version %s not recognized" % version

def redden(wave, flux, ebv_gal, ebv_host, z, R_gal=3.1, R_host=3.1,
      redlaw='ccm', strict_ccm=False):
   '''Artificially redden the spectral template to simulate dust reddening, a la
   Cardelli et al.
   
   Args:
      wave (float array): Input wavelength in Angstroms
      flux (float array): arbitrarily scaled SED flux
      ebv_gal (float): color excess to be applied in rest-frame (due to MW)
      ebv_host (floag): color excess to be applied at host redshift
      z (float): redshift of the host extinction
      R_gal (float): Ratio of total to selective absoption in V for restframe
                     component of extinction.
      R_host (float): Ratio of total to selective absorption in V for host
                      frame extinction.
      redlaw (str): Form of the dust extinction curve. Possible values are
                    'ccm', 'f99', or 'fm07'. See :mod:`snpy.utils.deredden`.
   
   Returns:
      float array: reddened flux.
   '''

   #First we redden due to galactic extinction:
   newflux = 1.0*flux
   # ebv_host is in the frame of the SN
   if ebv_host != 0:
      newflux,a,b = deredden.unred(wave, flux, -ebv_host, R_host,redlaw=redlaw,
            strict_ccm=strict_ccm)
   # ebv_gal is in the frame of the observer
   if ebv_gal != 0:
      newflux,a,b = deredden.unred(wave, newflux, -ebv_gal, R_gal, z, 
            redlaw=redlaw, strict_ccm=strict_ccm)

   return(newflux)

def K(wave, spec, f1, f2, z, photons=1):
   '''compute single K-correction based on a single spectrum and set of 
   observed and rest-frame filters.
   
   Args:
      wave (float array): input wavelength in Angstroms
      flux (float array): arbitrarily scaled flux
      f1 (filter instance): Rest-frame filter.
      f2 (filter instance): Observed filter. This could be the same as f1 or
                            a redder filter for cross-band K-correction
      z (float): redshift
      photons (bool): If True, fluxes are computed in units of photons rather
                      than energy (see Nugent+2002)
   
   Returns:
      2-tuple: (K,flag)

      * K: K-correction
      * flag: 1 -> success, 0->failed
   '''

   # The zero-points
   zpt1 = f1.zp
   zpt2 = f2.zp

   # compute the response through each filter
   f1flux_0 = f1.response(wave, spec, photons=photons)
   f2flux_z = f2.response(wave, spec, z=z, photons=photons)

   if f1flux_0 < 0 or f2flux_z <= 0:
      # Something clearly went wrong
      return (0.0, 0)
   else:
      # Finally calculate the cross-band K Correction
      # Since we blueshift the spectrum (instead of redshift the filter)
      # the sign of the 2.5 is now positive
      kf1f2 =  2.5*num.log10(1+z) + 2.5*num.log10(f1flux_0/f2flux_z) - \
               zpt1 + zpt2
      return (kf1f2, 1)

def kcorr(days, filter1, filter2, z, ebv_gal=0, ebv_host=0, R_gal=3.1, 
      R_host=3.1, version="H3", photons=1):
   '''Find the cross-band k-correction for a series of type Ia SED from
   SNooPy's catalog. These can be thought of as "empirical" K-corrections.
   
   Args:
      days (float array): epochs (t-T(Bmax)) to compute
      filter1 (str):  rest-frame filter
      filter2 (str):  observed filter. This can be the same as filter1,
                      or another, redder, filter for cross-band K-corrections
      z (float): redshift
      ebv_gal (float): restframe (foreground) color excess to be applied to
                       SED before computing K-corrections
      ebv_host (float): host-galaxy color excess to be applied to SED before
                       computing K-correction
      R_gal (float): Ratio of selective to total absorption at V for restframe
                     extinction.
      R_host (float): Ratio of selective to total absorption at V for host
                     extinction.
      version (str): Which SED sequence to use. See :func:`.get_SED`
      photons (bool): If True, compute fluxes in units of photons rather
                      than energy. Default is true and should be used unless
                      filter definition is in energy units.

   Returns
      2-tuple:  (K,mask)

      * K (float array): K-corrections
      * mask (bool array): True where K-corrections are valid.
                     
   '''

   if filter1 not in filters.fset:
      raise AttributeError, "filter %s not defined in filters module" % filter1
   if filter2 not in filters.fset:
      raise AttributeError, "filter %s not defined in filters module" % filter2

   kcorrs = []
   mask = []      # Masks the good values (1) and bad (not defined) values (0)
   # Loop through the list of days
   for day in days:
      day = int(day)
      spec_wav,spec_f = get_SED(day, version)
      if spec_wav is None:
         # print "Warning:  no spectra for day %d, setting Kxy=0" % day
         kcorrs.append(0.0)
         mask.append(0)
         continue

      # Do the reddening, if required
      if ebv_gal > 0 or ebv_host > 0:
         sp_f = redden(spec_wav, spec_f, ebv_gal, ebv_host, z, R_gal, R_host)
      else:
         sp_f = spec_f

      f1 = filters.fset[filter1]
      f2 = filters.fset[filter2]
      k,f = K(spec_wav, sp_f, f1, f2, z, photons=photons)
      kcorrs.append(k)
      mask.append(f)
   return(kcorrs,mask)

def kcorr_mangle2(waves, spectra, filts, mags, m_mask, restfilts, z, 
      colorfilts=None, full_output=0, **mopts): 
   '''Compute (cross-)band K-corrections with "mangling" using provided
   spectral SEDs. The SEDs are first multiplied by a smooth spline such that
   the synthetic colors match the observed colors.

   Args:
      waves (list of float arrays):  Input wavelengths in Angstroms
      spectra (list of float arrays):  Input fluxes in arbitrary units
      filts (list of str): list of observed filters
      mags (2d float array): Observed magnitude array indexed by
                             [spectrum index,filter index]
      m_mask (2d bool array): mask array indicating valid magnitudes. Indexed
                              by [spectrum index,filter index]
      restfilts (list of str): Rest-frame filters corresponing to filts.
      z (float):  redshift
      colorfilts (list of str): (optional) Sub set of filters to use in 
                              mangling colors (filters that have very similar
                              effective wavelengths can make for unstable
                              splines).
      full_output (bool):  If True, output more information than just the
                          K-corrections and mask.
      mopts (dict): All additional arguments to function are sent to 
                    :func:`snpy.mangle_spectrum.mangle_spectrum2`.
   
   Returns:
      tuple:

         * if not full_output: 2-tuple (K,mask):
            * K (flaot array):  K-corrections for filts
            * mask (bool array): mask of valid K-corrections
         * if full_output: 5-tuple (K,mask,anchors,factors,funcs)
            * anchors (float array): wavelengths of anchor points
            * factors (float array): factors corresponding to anchors
            * funcs (float array): mangling function evaluated at anchors
   '''
   if colorfilts is None:  colorfilts = filts
   for filter1 in filts + restfilts + colorfilts:
      if filter1 not in filters.fset:
         raise AttributeError, "filter %s not defined in filters module" % filter1

   if len(num.shape(waves)) < 2:
      scalar = 1
      waves = num.array([waves])
      spectra = num.array([spectra])
      mags = num.array([mags])
      m_mask = num.array([m_mask])
   else:
      scalar = 0


   kcorrs = []
   mask = []      # Masks the good values (1) and bad (not defined) values (0)
   waves_a = []
   manf_a = []
   factors_a = []

   for j in range(len(spectra)):
      kcorrs.append([])
      mask.append([])
      spec_wav,spec_f = waves[j],spectra[j]

      # Now determine which colors to use:
      fs = [colorfilts[i] for i in range(len(colorfilts)) if m_mask[j,i]]
      if len(fs) <= 1:
         # only one filter, so no color information, leave the SED alone:
         man_waves,man_spec_f,factors = spec_wav,spec_f,spec_wav*0.0+1.0
      else:
         #cs = num.compress(m_mask[j],mags[j])[0:-1] - \
         #      num.compress(m_mask[j],mags[j])[1:]
         ms = num.compress(m_mask[j],mags[j])
 
         # Now we mangle the spectrum:
         man_spec_f,man_waves,factors = mangle_spectrum2(spec_wav*(1+z),spec_f,fs, 
               ms, **mopts)
      if full_output:
         waves_a.append(man_waves)
         manf_a.append(man_spec_f)
         factors_a.append(factors)
 
      for i in range(len(filts)):
         f1 = filters.fset[restfilts[i]]
         f2 = filters.fset[filts[i]]

         k,f = K(spec_wav, man_spec_f[0], f1, f2, z)
         if f == 1:
            kcorrs[-1].append(k)
            mask[-1].append(len(fs))
         else:
            kcorrs[-1].append(0)
            mask[-1].append(0)
   kcorrs = num.array(kcorrs)
   mask = num.array(mask)

   if full_output:
      if scalar:
         return(kcorrs[0], mask[0], waves_a[0], factors_a[0], manf_a[0])
      else:
         return(kcorrs, mask, waves_a, factors_a, manf_a)
   else:
      if scalar:
         return(kcorrs[0],mask[0])
      else:
         return(kcorrs,mask)

def kcorr_mangle(days, filts, mags, m_mask, restfilts, z, version='H', 
      colorfilts=None, full_output=0, mepoch=False, **mopts):
   '''Compute (cross-)band K-corrections with "mangling" using built-in library
   of spectral SEDs. The SEDs are first multiplied by a smooth spline such that
   the synthetic colors match the observed colors.

   Args:
      days (float array): epochs (t-Tmax(B)) at which to compute K-corrections
      filts (list of str): list of observed filters
      mags (2d float array): Observed magnitude array indexed by
                             [spectrum index,filter index]
      m_mask (2d bool array): mask array indicating valid magnitudes. Indexed
                              by [spectrum index,filter index]
      restfilts (list of str): Rest-frame filters corresponing to filts.
      z (float):  redshift
      version (str): Specify which spectral sequence to use. See
                     :func:`.get_SED`.
      colorfilts (list of str): (optional) Sub set of filters to use in 
                              mangling colors (filters that have very similar
                              effective wavelengths can make for unstable
                              splines).
      full_output (bool):  If True, output more information than just the
                          K-corrections and mask.
      mepoch (bool): If True, a single mangling function is solved for
                     all epochs. EXPERIMENTAL.
      mopts (dict): All additional arguments to function are sent to 
                    :func:`snpy.mangle_spectrum.mangle_spectrum2`.
   
   Returns:
      tuple:

         * if not full_output: 2-tuple (K,mask):
            * K (flaot array):  K-corrections for filts
            * mask (bool array): mask of valid K-corrections
         * if full_output: 5-tuple (K,mask,anchors,factors,funcs)
            * anchors (float array): wavelengths of anchor points
            * factors (float array): factors corresponding to anchors
            * funcs (float array): mangling function evaluated at anchors
   '''

   if 'method' in mopts:
      method = mopts['method']
   else:
      method = 'tspline'

   if colorfilts is None:  colorfilts = filts
   for filter1 in filts + restfilts + colorfilts:
      if filter1 not in filters.fset:
         raise AttributeError, "filter %s not defined in filters module" % filter1

   kcorrs = []
   mask = []      # Masks the good values (1) and bad (not defined) values (0)
   m_opts = []
   Rts = []
   if debug:  mopts['verbose'] = 1

   if mepoch:
      # Doing multi epoch simultaneously with one mangling function...
      spec_wavs = []
      spec_fs = []
      sids = []
      for j in range(len(days)):
         day = int(days[j])
         s,f = get_SED(day, version)
         if s is None:
            spec_wavs.append(num.arange(980.,24981.0,10.))
            spec_fs.append(num.zeros((2401,), dtype=num.float64))
            sids.append(False)
         else:
            spec_wavs.append(s)
            spec_fs.append(f)
            sids.append(True)
      spec_wavs = num.array(spec_wavs)
      spec_fs = num.array(spec_fs)
      sids = num.array(sids)
      fs = [colorfilts[i] for i in range(len(colorfilts)) \
            if num.sometrue(m_mask[:,i])]
      if len(fs) <=1 :
         waves,man_spec_fs,factors = spec_wavs,spec_fs,spec_wavs*0.0+1.0
      else:
         #cs = mags[:,:-1] - mags[:,1:]
         #gids = m_mask[:,:-1]*m_mask[:,1:]
         #gids = gids*sids[:,num.newaxis]
         #cs[-gids] = 99.9   # flag invalid value
         gids = m_mask*sids[:,num.newaxis]
         ms = where(gids, ms, 99.9)
         man_spec_fs,waves,factors = mangle_spectrum2(spec_wavs*(1+z),
               spec_fs, fs, ms, **mopts)


      for j in range(len(days)):
         kcorrs.append([])
         mask.append([])
         if not sids[j]:
            kcorrs[-1] = num.zeros((len(filts),), dtype=num.float32)
            mask[-1] = num.zeros((len(filts),), dtype=num.int8)
            Rts.append(kcorrs[-1] - 1.0)
            m_opts.append(None)
            continue
         if full_output:
            args = {'sw':waves, 'sf':factors}
            for key in mopts:
               args[key] = mopts[key]
            m_opts.append(args)
 
         for i in range(len(filts)):
             f1 = filters.fset[restfilts[i]]
             f2 = filters.fset[filts[i]]
             k,f = K(spec_wavs[j], man_spec_fs[j], f1, f2, z)
             kcorrs[-1].append(k)
             mask[-1].append(0)
         Rts.append(R_obs_spectrum(filts, spec_wavs[j], man_spec_fs[j], z, 
            0.01, 0.0))
   else:
      for j in range(len(days)):
         kcorrs.append([])
         mask.append([])
         day = int(days[j])
         spec_wav,spec_f = get_SED(day, version)
         if spec_wav is None:
            # print "Warning:  no spectra for day %d, setting Kxy=0" % day
            kcorrs[-1] = num.zeros((len(filts),), dtype=num.float32)
            mask[-1] = num.zeros((len(filts),), dtype=num.int8)
            Rts.append(kcorrs[-1] - 1.0)
            m_opts.append(None)
            continue
 
         # Now determine which colors to use:
         fs = [colorfilts[i] for i in range(len(colorfilts)) if m_mask[j,i]]
         if len(fs) <= 1:
            # only one filter, so no color information, leave the SED alone:
            waves,man_spec_f,factors = spec_wav,spec_f,spec_wav*0.0+1.0
            man_spec_f = [man_spec_f]
         else:
            #cs = num.compress(m_mask[j],mags[j])[0:-1] - \
            #     num.compress(m_mask[j],mags[j])[1:]
            ms = num.compress(m_mask[j], mags[j])
            if debug:
               print "filters and colors for day %f:" % (days[j])
               print fs
               print ms[:-1]-ms[1:]
  
            # Now we mangle the spectrum.  Note, we are redshifting the spectrum
            # here, so do NOT set z in mangle_spectrum2.
            man_spec_f,waves,factors = mangle_spectrum2(spec_wav*(1+z),spec_f,
                  fs, ms, **mopts)
 
            if debug:  print "factors = ",factors
            if debug:
               # check the colors
               for i in range(len(fs)-1):
 
                  print "input color:  %s-%s = %f" % (fs[i],fs[i+1], ms[i]-ms[i+1]),
                  f1 = filters.fset[fs[i]]
                  f2 = filters.fset[fs[i+1]]
                  col = f1.synth_mag(spec_wav*(1+z), man_spec_f[0]) - \
                        f2.synth_mag(spec_wav*(1+z), man_spec_f[0])
                  print "  output color:  %f" % (col)
  
         if full_output:
            args = {'sw':waves, 'sf':factors}
            for key in mopts:
               args[key] = mopts[key]
            m_opts.append(args)
  
         for i in range(len(filts)):
            f1 = filters.fset[restfilts[i]]
            f2 = filters.fset[filts[i]]
            k,f = K(spec_wav,man_spec_f[0], f1, f2, z)
            kcorrs[-1].append(k)
            mask[-1].append(f)
         Rts.append(R_obs_spectrum(filts, spec_wav, man_spec_f[0], z, 0.01, 
            0.0))
   Rts = num.array(Rts)
   gids = num.greater(Rts, 0)
   Rtave = num.array([num.average(num.compress(gids[:,k], Rts[:,k])) \
         for k in range(len(gids[0]))])
   Rts = num.array([num.where(gids[i], Rts[i], Rtave) for i in range(len(gids))])

   kcorrs = num.array(kcorrs)
   mask = num.array(mask)

   if not full_output:
      return(kcorrs,mask)
   else:
      return(kcorrs, mask, Rts, m_opts)

def R_obs_abc(filter1, filter2, filter3, z, days, EBVhost, EBVgal, 
      Rv_host=3.1, Rv_gal=3.1, version='H', redlaw='f99', strict_ccm=False):
   '''Compute the observed value of the selective-to-total extinction, R,
   by applying an extinction curve to a set of library Ia spectral SEDs and
   computing synthetic photometry:

   .. math::

      A(\lambda_1) = R\ E(\lambda_2-\lambda_3)

   where

   .. math::

      E(\lambda_2-\lambda_3) = (m_{\lambda_2} - m_{\lambda_3}) - 
                                 (m_{\lambda_2} - m_{\lambda_3})_o 

   ie, the color excess for filters :math:`\lambda_2` and :math:`\lambda_3`.
   
   Args:
       filter1,filter2,filter3 (str): the 3 filters defining R
       z (float): redshift of host compoenent of extinction
       days (int array): epochs at which to comptue R (t-Tmax(B))
       EBVhost (float): host component of extinction to apply
       EBVgal (float): Milky-way (foreground) component of extinction to apply
       Rv_host (float): R_V for host component
       Rv_gal (float): R_V for MW component
       version (str): Version of Ia SED library. See :func:`.get_SED`.
       redlaw (str): Which reddening law to use. See :mod:`snpy.utils.deredden`
       strict_ccm (bool): If True and using CCM reddening law, do not apply
                          the corrections due to O'Donnel.'''
   try:
      N = len(days)
      outarr = 1
   except:
      days = [days]
      outarr = 0

   Rs = []
   for day in days:
      spec_wav,spec_f = get_SED(int(day), version)
      if spec_wav is None:
         Rs.append(99.9)
      # Redden the spectrum based on Cardelli et al. and assumed EBVgal + EBVhost
      red_f = redden(spec_wav, spec_f, EBVgal, EBVhost, z, Rv_gal, Rv_host,
            redlaw=redlaw, strict_ccm=strict_ccm)
      for filter in [filter1,filter2,filter3]:
         if filter not in filters.fset:
            raise AttributeError, "filter %s not defined in filters module" % filter
 
      # Now, we get the response across the filters:
      resp = {}
      resp_red = {}
      for filter in [filter1,filter2,filter3]:
         if filter not in resp:
            resp[filter] = filters.fset[filter].response(spec_wav, spec_f, z=z)
            resp_red[filter] = filters.fset[filter].response(spec_wav, red_f, z=z)
 
      A1 = -2.5*num.log10(resp_red[filter1]/resp[filter1])
      A2 = -2.5*num.log10(resp_red[filter2]/resp[filter2])
      A3 = -2.5*num.log10(resp_red[filter3]/resp[filter3])
      Rs.append(A1/(A2 - A3))
   if outarr:
      return(num.array(Rs))
   else:
      return(Rs[0])

def A_obs(filter, z, days, EBVhost, EBVgal, Rv_host=3.1, Rv_gal=3.1, 
      version='H3'):
   try:
      N = len(days)
      outarr = 1
   except:
      days = [days]
      outarr = 0

   As = []
   for day in days:
      if day < -19:  day=-18
      if day > 70:  day = 69
      spec_wav,spec_f = get_SED(int(day), version)
      if spec_wav is None:
         Rs.append(99.9)
         continue
      # rEDDEN the spectrum based on Cardelli et al. and assumed EBVgal+EBVhost
      red_f = redden(spec_wav, spec_f, EBVgal, EBVhost, z, Rv_gal, Rv_host)

      if filter not in filters.fset:
         raise AttributeError, "filter %s not defined in filters module" % filter
 
      # Now, we get the response across the filters:
      resp = filters.fset[filter].response(spec_wav, spec_f, z=z, photons=1)
      resp_red = filters.fset[filter].response(spec_wav, red_f, z=z, photons=1)
 
      As.append(-2.5*num.log10(resp_red/resp))
   if outarr:
      return(num.array(As))
   else:
      return(As[0])

def R_obs(filter, z, days, EBVhost, EBVgal, Rv_host=3.1, Rv_gal=3.1, 
      version='H', redlaw='f99', strict_ccm=False):
   '''Compute the 'true' value of R based on a fiducial value of Rv for both Galactic and
   host extinction and the SED of a supernova.  The filter is such that:

      A(filter) = R(filter)*E(B-V).
   '''
   try:
      N = len(days)
      outarr = 1
   except:
      days = [days]
      outarr = 0

   Rs = []
   for day in days:
      if day < -19:  day=-18
      if day > 70:  day = 69
      spec_wav,spec_f = get_SED(int(day), version)
      if spec_wav is None:
         Rs.append(99.9)
         continue
      # rEDDEN the spectrum based on Cardelli et al. and assumed EBVgal + EBVhost
      red_f = redden(spec_wav, spec_f, EBVgal, EBVhost, z, Rv_gal, Rv_host,
            redlaw=redlaw, strict_ccm=strict_ccm)

      if filter not in filters.fset:
         raise AttributeError, "filter %s not defined in filters module" % filter
 
      # Now, we get the response across the filters:
      resp = filters.fset[filter].response(spec_wav, spec_f, z=z, photons=1)
      resp_red = filters.fset[filter].response(spec_wav, red_f, z=z, photons=1)
 
      A_obs = -2.5*num.log10(resp_red/resp)
      Rs.append(A_obs/(EBVhost + EBVgal))
   if outarr:
      return(num.array(Rs))
   else:
      return(Rs[0])

def R_obs_spectrum(filts, wave, flux, z, EBVgal, EBVhost, Rv_gal=3.1, 
      Rv_host=3.1, redlaw='f99', strict_ccm=False):
   '''Compute the 'true' value of R based on a fiducial value of Rv for both Galactic and
   host extinction and the SED given by wave,flux for each filter in filters.  The filter 
   is such that: A(filter) = R(filter)*E(B-V).'''

   Rs = []
   if len(num.shape(filts)) == 0:
      outarr = 0
   else:
      outarr = 1

   # redden the spectrum based on EBVhost and EBVgal:
   red_f = redden(wave, flux, EBVgal, EBVhost, z=z, R_gal=Rv_gal, 
         R_host=Rv_host, redlaw=redlaw, strict_ccm=strict_ccm)

   # Compute absoption in each filter:
   for filter in filts:
      # Redden the spectrum based on Cardelli et al. and assumed EBVgal + EBVhost
      if filter not in filters.fset:
         raise AttributeError, "filter %s not defined in filters module" % filter
 
      # Now, we get the response across the filters:
      resp = filters.fset[filter].response(wave, flux, z=z, photons=1)
      resp_red = filters.fset[filter].response(wave, red_f, z=z, photons=1)
 
      A_obs = -2.5*num.log10(resp_red/resp)
      #Rs.append(A_obs/EBV_obs)
      Rs.append(A_obs/(EBVgal + EBVhost))
   if outarr:
      return(num.array(Rs))
   else:
      return(Rs[0])

