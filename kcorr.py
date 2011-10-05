#!/usr/bin/env python

# T_NEWXKCORR -- Calculate cross-band K Correction from a 1-d spectrum.
# This version converts the spectrum to photons before multiplying by the filter trans.
#  
#   Program written by MMP
#     Re-written in python by CRB
#  Added reddening correction for galactic and host extinctions

import os,sys,string,re
#import pygplot
import numpy.oldnumeric as num
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
   h_wav = f['CRVAL1'] + (num.arange(f['NAXIS1'], typecode=num.Float32) - \
         f['CRPIX1'] - 1)*f['CDELT1']
   f.close()
   # Hsiao's new OPT+NIR uberspectrum
   f = FITS.FITS(os.path.join(spec_base, 'Hsiao_SED_V3.fits'))
   h3_sed = f.data()
   h3_wav = f['CRVAL1'] + (num.arange(f['NAXIS1'], typecode=num.Float32) - \
         f['CRPIX1'] - 1)*f['CDELT1']
   f.close()
   # Nugent's uberspectrum:
   f = FITS.FITS(os.path.join(spec_base, "Nugent_SED.fits"))
   n_sed = f.data()
   n_wav = f['CRVAL1'] + (num.arange(f['NAXIS1'], typecode=num.Float32) - \
         f['CRPIX1'] - 1)*f['CDELT1']
   # Nugent's 91bg-like SED:
   f = FITS.FITS(os.path.join(spec_base, "Nugent_91bg_SED.fits"))
   n91_sed = f.data()
   n91_wav = f['CRVAL1'] + (num.arange(f['NAXIS1'], typecode=num.Float32) - \
         f['CRPIX1'] - 1)*f['CDELT1']
else:
   # Hsiao's uberspectrum:
   f = pyfits.open(os.path.join(spec_base, 'Hsiao_SED_V2.fits'))
   h_sed = f[0].data
   head = f[0].header
   h_wav = head['CRVAL1'] + (num.arange(head['NAXIS1'],typecode=num.Float32) - \
         head['CRPIX1'] - 1)*head['CDELT1']
   f.close()
   # Hsiao's new OPT+NIR uberspectrum
   f = pyfits.open(os.path.join(spec_base, 'Hsiao_SED_V3.fits'))
   h3_sed = f[0].data
   head = f[0].header
   h3_wav = head['CRVAL1']+(num.arange(head['NAXIS1'],typecode=num.Float32) - \
         head['CRPIX1'] - 1)*head['CDELT1']
   f.close()
   # Nugent's uberspectrum:
   f = pyfits.open(os.path.join(spec_base, "Nugent_SED.fits"))
   n_sed = f[0].data
   head = f[0].header
   n_wav = head['CRVAL1']+(num.arange(head['NAXIS1'],typecode=num.Float32) - \
         head['CRPIX1'] - 1)*head['CDELT1']
   # Nugent's 91bg-like SED:
   f = pyfits.open(os.path.join(spec_base, "Nugent_91bg_SED.fits"))
   n91_sed = f[0].data
   head = f[0].header
   n91_wav = head['CRVAL1']+(num.arange(head['NAXIS1'],typecode=num.Float32) - \
         head['CRPIX1'] - 1)*head['CDELT1']

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

def redden(wave, flux, ebv_gal, ebv_host, z, R_gal=3.1, R_host=3.1):
   '''Artificially redden the spectral template to simulate dust reddening, a la
   Cardelli et al.  Unless over-ridden, the standard reddening coefficients is 
   assumed (Rv=3.1) for both the Milky Way and host.'''
   
   #First we redden due to galactic extinction:
   newflux = 1.0*flux
   # ebv_host is in the frame of the SN
   if ebv_host != 0:
      newflux,a,b = deredden.unred(wave, flux, -ebv_host, R_host)
   # ebv_gal is in the frame of the observer
   if ebv_gal != 0:
      newflux,a,b = deredden.unred(wave, newflux, -ebv_gal, R_gal, z)

   return(newflux)

def kcorr(days, filter1, filter2, z, ebv_gal=0, ebv_host=0, R_gal=3.1, R_host=3.1,
      version="H3"):
   '''Find the cross-band k-correction for rest band filter1 to observed
   band filter2 at redshift z.  In other words, filter1 should be bluer than filter2
   if z is a red-shift.'''

   if filter1 not in filters.fset:
      raise AttributeError, "filter %s not defined in filters module" % filter1
   if filter2 not in filters.fset:
      raise AttributeError, "filter %s not defined in filters module" % filter2
   f1 = filters.fset[filter1]
   f2 = filters.fset[filter2]
   zpt1 = f1.zp
   zpt2 = f2.zp

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

      # Now compute the fluxes using Simpson's composite rule:
      #  But now we keep the filter response unchanged, so instead
      #  we blue-shift the SED
      f1flux_0 = f1.response(spec_wav, sp_f)
      f2flux_z = f2.response(spec_wav, sp_f, z=z)
      if f1flux_0 < 0 or f2flux_z <= 0:
         kcorrs.append(0.0)
         mask.append(0)
      else:
         # Finally calculate the cross-band K Correction
         # Since we blueshift the spectrum (instead of redshift the filter)
         # the sign of the 2.5 is now positive
         kf1f2 =  2.5*num.log10(1+z) + 2.5*num.log10(f1flux_0/f2flux_z) - \
                  zpt1 + zpt2
         kcorrs.append(kf1f2)
         mask.append(1)
   return(kcorrs,mask)

def kcorr2(wave, spectrum, filter1, filter2, z):
   '''Find the cross-band k-correction for rest band filter1 to observed band
    filter2 at redshift z for an in put spectrum given by wave,spectrum. 
    The SED should be in the REST frame.'''

   if filter1 not in filters.fset:
      raise AttributeError, "filter %s not defined in filters module" % filter1
   if filter2 not in filters.fset:
      raise AttributeError, "filter %s not defined in filters module" % filter2
   f1 = filters.fset[filter1]
   f2 = filters.fset[filter2]
   zpt1 = f1.zp
   zpt2 = f2.zp

   kcorrs = []
   mask = []      # Masks the good values (1) and bad (not defined) values (0)
   spec_wav,spec_f = wave,spectrum

   # Now compute the fluxes using Simpson's composite rule:
   f1flux_0 = f1.response(spec_wav, spec_f, z=0)
   f2flux_z = f2.response(spec_wav, spec_f, z=z)
   if f1flux_0 < 0 or f2flux_z <= 0:
      kcorr = 0
      mask = 0
   else:
      # Finally calculate the cross-band K Correction
      kf1f2 =  2.5*num.log10(1+z) + 2.5*num.log10(f1flux_0/f2flux_z) - \
               zpt1 + zpt2
      kcorr = kf1f2
      mask = 1
   return(kcorr,mask)

def kcorr_mangle3(waves, spectra, filts, mags, m_mask, restfilts, z, colorfilts=None,
      full_output=0, **mopts): 
   '''Find the (cross-band) k-correction for each filter in restfilts to
   observed corresponding filter in fitls at redshift z, using the provided
   colorfilts (or filts if colorfilts is not defined) and magnitudes to mangle
   the SED provided in waves,spectra.  The SED must be in the REST frame.
   The array mags and mask should have
   dimensions (len(spectra),len( colorfilts)) so that mags[i,j] coorespond to
   the magnitude for spectra[i] in filter filts[j].  The mask is used to
   determine which magnitudes are good and which are bad (for whatever reason).
   This version is good if you have your own spectra to mangle instead of the
   Nugent or Hsiao templates.'''

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
         waves,man_spec_f,factors = spec_wav,spec_f,spec_wav*0.0+1.0
      else:
         cs = num.compress(m_mask[j],mags[j])[0:-1] - \
               num.compress(m_mask[j],mags[j])[1:]
 
         # Now we mangle the spectrum:
         man_spec_f,waves,factors = mangle_spectrum2(spec_wav*(1+z),spec_f,fs, 
               cs, **mopts)
      if full_output:
         waves_a.append(waves)
         manf_a.append(man_spec_f)
         factors_a.append(factors)
 
      # Let's plot these guys out
      #if debug:
      #   p = pygplot.Plot(device='/XWIN', xrange=[3500,9500])
      #   p.line(spec_wav, spec_f/max(spec_f), color='blue')
      #   p.line(spec_wav, man_spec_f/max(spec_f), color='green')
      #   p.point(waves/(1+z), factors, symbol=3, size=2,color='orange')
      #   tck = scipy.interpolate.splrep(waves/(1+z), factors, k=3, s=0)
      #   xx = num.arange(waves[0]/(1+z), waves[-1]/(1+z), 10)
      #   p.line(xx, scipy.interpolate.splev(xx, tck), color='orange')
 
      for i in range(len(filts)):
         f1 = filters.fset[restfilts[i]]
         zpt1 = f1.zp
         f2 = filters.fset[filts[i]]
         zpt2 = f2.zp
         # Now compute the fluxes using Simpson's composite rule:
         f1flux_0 = f1.response(spec_wav, man_spec_f, z=0)
         f2flux_z = f2.response(spec_wav, man_spec_f, z=z)
         #if debug:
         #   p.line(f1.wave, f1.resp)
         #   p.line(f2.wave, f2.resp)
         #   p.line(f2.wave/(1.0+z), f2.resp, color='red')
         #   p.plot()
         #   p.close()
         #   del p.lines[-3:]
         if f1flux_0 < 0 or f2flux_z < 0:
            kcorrs[-1].append(0.0)
            mask[-1].append(0)
         else:
            # Finally calculate the cross-band K Correction
            kf1f2 = 2.5*num.log10(1+z) + 2.5*num.log10(f1flux_0/f2flux_z) - \
                     zpt1 + zpt2
            kcorrs[-1].append(kf1f2)
            mask[-1].append(len(fs))
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

def kcorr_mangle2(days, filts, mags, m_mask, restfilts, z, version='H', colorfilts=None,
      full_output=0, normfilter=None, **mopts):
   '''Find the (cross-band) k-correction for each filter in [restfilts] to
   observed corresponding filter in [filts] at redshift [z], using the provided
   [colorfilts] (or [filts] if [colorfilts] is not defined) and magnitudes to mangle
   the SED.  The array [mags] and [mask] should have dimensions (len(days),len(
   colorfilts)) so that mags[i,j] correspond to the magnitude on day days[i] in
   filter filts[j].  The mask is used to determine which magnitudes are good
   and which are bad (for whatever reason).  Use this version if your data
   needs to be masked in any way.  If [full_output]=True, returm a list of
   [waves,factors] of

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

   for j in range(len(days)):
      kcorrs.append([])
      mask.append([])
      day = int(days[j])
      spec_wav,spec_f = get_SED(day, version)
      if spec_wav is None:
         # print "Warning:  no spectra for day %d, setting Kxy=0" % day
         kcorrs[-1] = num.zeros((len(filts),), typecode=num.Float32)
         mask[-1] = num.zeros((len(filts),))
         Rts.append(kcorrs[-1] - 1.0)
         m_opts.append(None)
         continue

      # Now determine which colors to use:
      fs = [colorfilts[i] for i in range(len(colorfilts)) if m_mask[j,i]]
      if len(fs) <= 1:
         # only one filter, so no color information, leave the SED alone:
         waves,man_spec_f,factors = spec_wav,spec_f,spec_wav*0.0+1.0
      else:
         cs = num.compress(m_mask[j],mags[j])[0:-1] - num.compress(m_mask[j],mags[j])[1:]
         if debug:
            print "filters and colors for day %f:" % (days[j])
            print fs
            print cs
 
         # Now we mangle the spectrum:
         man_spec_f,waves,factors = mangle_spectrum2(spec_wav*(1+z),spec_f,fs, 
               cs, **mopts)
 
         if debug:  print "factors = ",factors
         if debug:
            # check the colors
            for i in range(len(fs)-1):

               print "input color:  %s-%s = %f" % (fs[i],fs[i+1], cs[i]),
               f1 = filters.fset[fs[i]]
               f2 = filters.fset[fs[i+1]]
               col = f1.synth_mag(spec_wav*(1+z), man_spec_f) - \
                     f2.synth_mag(spec_wav*(1+z), man_spec_f)
               print "  output color:  %f" % (col)
 
      if full_output:
         args = {'sw':waves, 'sf':factors}
         for key in mopts:
            args[key] = mopts[key]
         m_opts.append(args)
      # Let's plot these guys out
      #if debug:
      #   p = pygplot.Plot(device='/XWIN', yrange=[0,2])
      #   p.line(spec_wav, spec_f/max(spec_f), color='blue')
      #   p.line(spec_wav, man_spec_f/max(spec_f), color='green')
      #   p.point(waves/(1+z), factors, color='orange', symbol=4)
      #   p.line(spec_wav, man_spec_f/spec_f, color='orange')
 
      for i in range(len(filts)):
         f1 = filters.fset[restfilts[i]]
         zpt1 = f1.zp
         f2 = filters.fset[filts[i]]
         zpt2 = f2.zp
         # Now compute the fluxes using Simpson's composite rule:
         f1flux_0 = f1.response(spec_wav, man_spec_f, z=0)
         f2flux_z = f2.response(spec_wav, man_spec_f, z=z)
         #if debug:
         #   p.line(f1.wave, f1.resp)
         #   p.line(f2.wave, f2.resp)
         #   p.line(f2.wave/(1.0+z), f2.resp, color='red')
         #   p.plot()
         #   p.close()
         #   del p.lines[-3:]
         if f1flux_0 < 0 or f2flux_z < 0:
            kcorrs[-1].append(0.0)
            mask[-1].append(0)
         else:
            # Finally calculate the cross-band K Correction
            kf1f2 = 2.5*num.log10(1+z) + 2.5*num.log10(f1flux_0/f2flux_z) - \
                     zpt1 + zpt2
            kcorrs[-1].append(kf1f2)
            mask[-1].append(len(fs))
      Rts.append(R_obs_spectrum(filts, spec_wav, man_spec_f, z, 0.01, 0.0))
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

def kcorr_mangle(days, filts, colors,  restfilts, z, version='H', **mopts):
   '''Find the cross-band k-correction for each filter in restfilts to observed
   corresponding filter in fitls at redshift z, using the provided filters and 
   colors to mangle the SED.  len(filts) must be one more than len(colors) and
   be equal to len(restfilts).  colors must be 2D with len(colors) = len(days)
   and len(colors[i])=len(filts)-1.  The argument days must be the epoch in the
   frame of the supernova.'''

   for filter1 in filts + restfilts:
      if filter1 not in filters.fset:
         raise AttributeError, "filter %s not defined in filters module" % \
               filter1

   kcorrs = []
   mask = []      # Masks the good values (1) and bad (not defined) values (0)

   j = 0
   for day in days:
      kcorrs.append([])
      mask.append([])
      day = int(day)
      spec_wav,spec_f = get_SED(day, version)
      if spec_wav is None:
         # print "Warning:  no spectra for day %d, setting Kxy=0" % day
         kcorrs[-1] = num.zeros((len(filts),), typecode=num.Float32)
         mask[-1] = num.zeros((len(filts),))
         continue
 
      # Now we mangle the spectrum:
      man_spec_f,waves,factors = mangle_spectrum2(spec_wav*(1+z),spec_f,filts,
           colors[j], **mopts)

      if debug:
         # check the colors
         for i in range(len(filts)-1):
            print "input color:  %s-%s = %f" % (filts[i],filts[i+1],colors[i]),
            f1 = filters.fset[filts[i]]
            f2 = filters.fset[filts[i+1]]
            col = f1.synth_mag(spec_wav*(1+z), man_spec_f) - \
                  f2.synth_mag(spec_wav*(1+z), man_spec_f)
            print "  output color:  %f" % (col)
 
      # Let's plot these guys out
      #if debug:
      #   p = pygplot.Plot(device='/XWIN')
      #   p.line(spec_wav, spec_f/max(spec_f), color='blue')
      #   p.line(spec_wav, man_spec_f/max(spec_f), color='green')
      #   #p.line(waves/(1+z), factors, color='orange')
      #   p.point(waves/(1+z), factors, symbol=3, size=2,color='orange')
      #   tck = scipy.interpolate.splrep(waves/(1+z), factors, k=3, s=0)
      #   xx = num.arange(waves[0]/(1+z), waves[-1]/(1+z), 10)
      #   p.line(xx, scipy.interpolate.splev(xx, tck), color='orange')
 
      for i in range(len(filts)):
         f1 = filters.fset[restfilts[i]]
         zpt1 = f1.zp
         f2 = filters.fset[filts[i]]
         zpt2 = f2.zp
         # Now compute the fluxes using Simpson's composite rule:
         f1flux_0 = f1.response(spec_wav, man_spec_f, z=0)
         f2flux_z = f2.response(spec_wav, man_spec_f, z=z)
         #if debug:
         #   p.line(f1.wave, f1.resp)
         #   p.line(f2.wave, f2.resp)
         #   p.line(f2.wave/(1.0+z), f2.resp, color='red')
         #   p.plot()
         #   p.close()
         #   del p.lines[-3:]
         if f1flux_0 < 0 or f2flux_z < 0:
            kcorrs[-1].append(0.0)
            mask[-1].append(0)
         else:
            # Finally calculate the cross-band K Correction
            kf1f2 = 2.5*num.log10(1+z) + 2.5*num.log10(f1flux_0/f2flux_z) - \
                     zpt1 + zpt2
            kcorrs[-1].append(kf1f2)
            mask[-1].append(1)
      j = j + 1
   kcorrs = num.array(kcorrs)
   mask = num.array(mask)

   return(kcorrs,mask)

def R_obs_abc(filter1, filter2, filter3, z, days, EBVhost, EBVgal, Rv_host=3.1, Rv_gal=3.1, 
      version='H'):
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

