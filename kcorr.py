#!/usr/bin/env python
'''

A module for computing k-corrections.  This module contains the following 
functions:

   get_SED(day, version):    Retrieve the SED of a SNIa.  The version can
                           refer to : 'H':  Hsiao (published)
                                      'H3': Hsiao (unpublished with NIR)
                                      'N':  Nugent
                                      '91bg': Nugent's SN1991bg template
   redden:  Use CCM+ODonnel to redden an SED.  Optionally redden by both a
            local(z=0) reddening law and high-z reddening law
   K:  compute the k-correction for a single spectrum and filter combo
   kcorr: compute the k-correction (without mangling) for a list of
          epochs
   kcorr_mangle:  compute the k-corrections (with mangling) for a list of
          epochs
   kcorr_mangle2:  comptute the k-corrections (with mangling) for a list of
          spectra.

Original comments from Mark's IRAF code:
 T_NEWXKCORR -- Calculate cross-band K Correction from a 1-d spectrum.
 This version converts the spectrum to photons before multiplying by the filter trans.
  
   Program written by MMP
     Re-written in python by CRB
  Added reddening correction for galactic and host extinctions'''

import os,sys,string,re
#import pygplot
import numpy as num
import scipy.interpolate
from utils import deredden
try:
   import FITS
   have_FITS=1
except ImportError:
   try:
      import pyfits
      have_FITS=0
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
if have_FITS:
   # Hsiao's uberspectrum:
   f = FITS.FITS(os.path.join(spec_base, 'Hsiao_SED_V2.fits'))
   h_sed = f.data()
   h_wav = f['CRVAL1'] + (num.arange(f['NAXIS1'], dtype=num.float32) - \
         f['CRPIX1'] + 1)*f['CDELT1']
   f.close()
   # Hsiao's new OPT+NIR uberspectrum
   f = FITS.FITS(os.path.join(spec_base, 'Hsiao_SED_V3.fits'))
   h3_sed = f.data()
   h3_wav = f['CRVAL1'] + (num.arange(f['NAXIS1'], dtype=num.float32) - \
         f['CRPIX1'] + 1)*f['CDELT1']
   f.close()
   # Nugent's uberspectrum:
   f = FITS.FITS(os.path.join(spec_base, "Nugent_SED.fits"))
   n_sed = f.data()
   n_wav = f['CRVAL1'] + (num.arange(f['NAXIS1'], dtype=num.float32) - \
         f['CRPIX1'] + 1)*f['CDELT1']
   # Nugent's 91bg-like SED:
   f = FITS.FITS(os.path.join(spec_base, "Nugent_91bg_SED.fits"))
   n91_sed = f.data()
   n91_wav = f['CRVAL1'] + (num.arange(f['NAXIS1'], dtype=num.float32) - \
         f['CRPIX1'] + 1)*f['CDELT1']
else:
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

def get_SED(day, version='H3'):
   '''Retrieve the SED for a SN on day 'day' (where day=0 corresponds to Bmax). 
   If version = 'H', use Eric Hsiao's version.  Otherwise, use Peter Nugent's.'''
   if version in ['H','H3','N']:
      if day < -19 or day > 70:  return (None,None)
      if version == 'H':
         return (h_wav, h_sed[day+20,:])
      elif version == 'H3':
         return (h3_wav, h3_sed[day+20,:])
      elif version == 'N':
         return (n_wav, n_sed[day+19,:])
   elif version == '91bg':
      if day < -13 or day > 100:  return(None,None)
      return (n91_wav, n91_sed[day+13, :])
   else:
      raise AttributeError, "version %s not recognized" % version

def redden(wave, flux, ebv_gal, ebv_host, z, R_gal=3.1, R_host=3.1,
      redlaw='ccm', strict_ccm=False):
   '''Artificially redden the spectral template to simulate dust reddening, a la
   Cardelli et al.  Unless over-ridden, the standard reddening coefficients is 
   assumed (Rv=3.1) for both the Milky Way and host.  You can choose which
   redening law (CCM or Fitzpatrick) by specifying redlaw='ccm' or
   redlaw='fm', respectively.'''
   
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
   '''compute single K-correction based on spectrum defined by [wave] and 
   [spec].  [f1] is the filter instance of the rest-frame filter.  [f2] is a 
   filter instance of the observer-frame filter.  So f1 should be either be 
   the same as f2 (for low-z) or bluer than f2.'''

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
   '''Find the cross-band k-correction for rest band [filter1] to observed
   band [filter2] at redshift [z].  In other words, [filter1] should be bluer 
   than [filter2] if [z] is a at red-shift.'''

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
   '''Find the (cross-band) k-correction for each filter in [restfilts] to
   observed corresponding filter in [fitls] at redshift [z], using the provided
   [colorfilts] (or [filts] if [colorfilts] is not defined) and magnitudes
   [mags] to mangle the SED provided in [waves],[spectra].  The SED must be in
   the REST frame.  The array [mags] and [mask] should have dimensions
   (len([spectra]),len( [colorfilts])) so that mags[i,j] coorespond to the
   magnitude for spectra[i] in filter filts[j].  The mask is used to determine
   which magnitudes are good and which are bad (for whatever reason).  This
   version is good if you have your own spectra to mangle instead of the Nugent
   or Hsiao templates.'''

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
   '''Find the (cross-band) k-correction for each filter in [restfilts] to
   observed corresponding filter in [filts] at redshift [z], using the provided
   [colorfilts] (or [filts] if [colorfilts] is not defined) and magnitudes to
   mangle the SED.  The array [mags] and [mask] should have dimensions
   (len(days),len( colorfilts)) so that mags[i,j] correspond to the magnitude
   on day days[i] in filter filts[j].  The mask is used to determine which
   magnitudes are good and which are bad (for whatever reason).  Use this
   version if your data needs to be masked in any way.  If [full_output]=True,
   returm a list of [waves,factors] of

   Output: 
      kcorrs,mask      array of k-corrections and associated masks

      if full_output=1: 
      kcorrs,mask,Rts,m_opts   Rts are the observed reddening
         coefficients computed from the mangled SED.  m_opts is a list of
         dictionaries that can be used to re compute the mangled spectrum.  For
         example: apply_mangle(wave,spec,**m_opts[2]) will give back the mangled
         spectrum for day days[2].  '''

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
      Rv_host=3.1, Rv_gal=3.1, version='H'):
   '''Compute the observed value of R based on a fiducial value of Rv for 
   both Galactic and host extinction and the SED of a supernova.  filter[1-3]
   are used in the sense that absorption in band filter1 is:
      A(filter1) = R(filter1,filter2,filter3)*E(filter2-filter3)
   where
      E(filter2-filter3) = (filter2-filter3) - (filter2-filter3)_0,
   ie, the color excess for filters filter2 and filter3.'''
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
      red_f = redden(spec_wav, spec_f, EBVgal, EBVhost, z, Rv_gal, Rv_host)
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

def R_obs(filter, z, days, EBVhost, EBVgal, Rv_host=3.1, Rv_gal=3.1, version='H'):
   '''Compute the 'true' value of R based on a fiducial value of Rv for both Galactic and
   host extinction and the SED of a supernova.  The filter is such that:
      A(filter) = R(filter)*E(B-V).'''
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
      red_f = redden(spec_wav, spec_f, EBVgal, EBVhost, z, Rv_gal, Rv_host)

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

def R_obs_spectrum(filts, wave, flux, z, EBVgal, EBVhost, Rv_gal=3.1, Rv_host=3.1):
   '''Compute the 'true' value of R based on a fiducial value of Rv for both Galactic and
   host extinction and the SED given by wave,flux for each filter in filters.  The filter 
   is such that: A(filter) = R(filter)*E(B-V).'''

   Rs = []
   if len(num.shape(filts)) == 0:
      outarr = 0
   else:
      outarr = 1

   # redden the spectrum based on EBVhost and EBVgal:
   red_f = redden(wave, flux, EBVgal, EBVhost, z=z, R_gal=Rv_gal, R_host=Rv_host)

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

