'''Module to compute the bolometric luminosity of a SNIa.

There are serveral steps:
   1) Use color_model (or other) to measure E(B-V) and R_V. This is done
      exterior to the bolometric function (and is given as input)
   2) mangle the Hsiao spectra to match the observed (unreddened) colors,
      which are computed either:
       a) by using a fit template model (from color_model, say)
       b) by using a pre-computed interpolation for each filter
       c) only using observed colors
   3) Normalize the spectra to match one filter's absolute flux
   4) de-reddeng flux
   5) Integrate spectra over (l1,l2)
   6) Correct for distance.
'''

from snpy.filters import fset, ch, standards
from .kcorr import get_SED
from . import kcorr
from .filters import standards
Vega = standards['VegaB']
from .utils import deredden
from .mangle_spectrum import mangle_spectrum2
from scipy.integrate import trapz
from scipy.interpolate import splrep,splev
from numpy import *
import sys
import types

def log(msg):
   sys.stderr.write(msg+"\n")

def bolometric_SED(sn, bands=None, lam1=None, lam2=None, refband=None,
              tmin=None, tmax=None,
              EBVhost=None, Rv=None, redlaw=None, extrap_red='RJ',
              Tmax=None, interpolate=None, interp_all=False, extrapolate=False,
              mopts={}, SED='H3', DM=None, cosmo='LambdaCDM', use_stretch=True, 
              extrap_SED=True, extra_output=False, verbose=False):

   w,f = get_SED(0, version='H3')
   if verbose: log("Starting bolometric calculation for %s\n" % sn.name)

   if EBVhost is None:
      EBVhost = getattr(sn, 'EBVhost', None)
      if EBVhost is None:
         raise AttributeError("Error:  you must either specify E(B-V) or "\
                               "fit a model that has EBVhost as a parameter")
   if Rv is None:
      Rv = getattr(sn, 'Rv', None)
      if Rv is None:
         raise AttributeError("Error:  you must either specify Rv or "\
                               "fit a model that has Rv as a parameter")

   if redlaw is None:
      redlaw = getattr(sn, 'redlaw', None)
      if redlaw is None:
         raise AttributeError("Error:  you must either specify redlaw or "\
                               "fit a model that has redlaw is defined")

   if bands is None:
      bands = getattr(sn.model, '_fbands', None)
      if bands is None:
         bands = list(sn.data.keys())
   for b in bands:
      if b not in sn.data:
         raise AttributeError("band %s not defined in data set" % (b))
   # Bands must be increasing in wavelength
   eff_waves = array([fset[b].eff_wave(w,f) for b in bands])
   sids = argsort(eff_waves)
   bands = [bands[i] for i in sids]
   pars0 = {}
   for b in bands:  pars0[b] = 1.0

   # Get itegration limits in SN restframe. So if we get it from the observer
   #  frame filter limits, we need to divide by (1+z)
   if lam1 is None:
      lam1 = array([fset[b].waverange()[0] for b in bands]).min()/(1+sn.z)
   if lam2 is None:
      lam2 = array([fset[b].waverange()[1] for b in bands]).max()/(1+sn.z)

   # Allow a list of lam1,lam2 which will allow to sub-integrate the SED
   # if you want to know what fraction comes out in a particular wavelength
   # regime
   scalar = False
   if len(shape(lam1)) == 0:
      scalar = True
      lam1 = [lam1]
   if len(shape(lam2)) == 0:
      lam2 = [lam2]

   if refband is not None:
      if refband not in bands:
         raise ValueError("refband %s is not one of your observed filters" % \
               refband)


   # We need a time of maximum to set the scale of the Hsiao templates
   if Tmax is None:
      Tmax = getattr(sn, 'Tmax', None)
      if Tmax is None:
         raise ValueError("You must supply a time of B maximum or fit a "\
               "model that has Tmax as a parameter")

   # Now check that we can interpolate if needed
   if interpolate is not None:
      if interpolate == 'spline':
         if verbose: log("   Using spline interpolation")
         for b in bands:
            if getattr(sn.data[b], 'interp', None) is None:
               raise ValueError("You asked for spline interpolation, but "\
                     "filter %s has no interpolator defined" % b)
      else:
         if verbose: log("   Using model interpolation")
         for b in bands:
            if b not in sn.model._fbands:
               raise ValueError("You asked for model interpolation, but "\
                     "filter %s was not fit with the model" % b)
   else:
      if verbose: log("   Not using interpolation")

   if type(SED) is type(""):
      # Assume it is a spectrum by name
      if SED in ['H3','H','N','91bg']:
         fSED = lambda x: get_SED(x, version=SED, extrapolate=True)
         if verbose: log("   Using SED template '%s'" % SED)
      elif SED in standards:
         fSED = lambda x: (standards[SED].wave,standards[SED].flux)
         if verbose: log("   Using standards['%s'] to compute effective wavelengths" % (SED))
      else:
         raise KeyError("SED '%s' not found in standards database" % SED)
   elif type(SED) in [list,tuple]:
      if len(SED) != 2:
         raise ValueError("SED must be tuple or list of length 2")
      fSED = lambda x: SED
   elif type(SED) is types.FunctionType:
      try:
         w,f = SED(0)
      except:
         raise ValueError("If SED is a function, it must take single" \
               " argument (epoch) and return (wave,flux) tuple")
      fSED = SED
   else:
      raise ValueError("Unrecognized type (%s) for SED" % (type(SED)))


   s = 1.0
   if use_stretch:
      s = getattr(sn, 'st', None)
      dm15 = getattr(sn, 'dm15', None)
      if s is None:
         if dm15 is None:
            raise ValueError("If you want to apply a stretch to the SED's, "\
                  "you must fit a model that uses dm15 or st as a parameter")
         if dm15 > 1.7:
            if verbose: log("Warning:  dm15 > 1.7. Hsiao template is not "
                  "compatible with fast decliners. Proceed a your own risk")
            s = kcorr.dm152s(1.7)
         elif dm15 < 0.7:
            if verbose: log("Warning:  dm15 < 0.7. This is a very slow "
                  "decliner. Proceed at your own risk")
            s = kcorr.dm152s(0.7)
         else:
            s = kcorr.dm152s(dm15)
      if verbose: log("   Using a stretch of %f" % s)

   # Now, we build up the times at which we will be integrating and the 
   # assocuated fluxes. We do things in the frame of the SN and de-stretch
   # the times.
   res = sn.get_mag_table(bands)
   ts = (res['MJD'] - Tmax)/(1+sn.z)
   # Restrict to valid interval of Eric's templates if extrap_SED is False
   if not extrap_SED:
      gids = greater_equal(ts/s, -19)*less_equal(ts/s,70)
   else:
      gids = ~isnan(ts/s)
   ts = ts[gids]
   mags = []
   masks = []
   for b in bands:
      mags.append(res[b][gids])
      masks.append(less(res[b][gids],90))

   if interpolate == 'spline':
      # we fill in (where we can) missing data using splines
      for i,b in enumerate(bands):
         # interpolation is in the absolute time, so use MJD
         mag,mask = sn.data[b].interp(res['MJD'][gids])
         mags[i] = where(masks[i], mags[i], mag)
         if extrapolate:
            masks[i] = -isnan(mags[i])
         else:
            masks[i] = masks[i] + mask

   elif interpolate == 'model':
      # we fill in (where we can) missing data using the model
      for i,b in enumerate(bands):
         mag,emag,mask = sn.model(b, res['MJD'][gids], extrap=extrapolate)
         mags[i] = where(masks[i], mags[i], mag)
         masks[i] = masks[i] + mask

   mags = transpose(array(mags))
   masks = transpose(array(masks))
   if verbose:  log("   Working on data matrix with size (%d,%d)" % mags.shape)
   # see if any days have zero data
   gids = greater(sum(masks, axis=1), 0)
   mags = mags[gids,:]
   masks = masks[gids,:]
   ts = ts[gids]

   # restrict on wanted interval if necessary
   if tmin is not None and tmin > ts.min():
      gids = greater(ts, tmin)
   if tmax is not None and tmax < ts.max():
      gids = less(ts, tmax)
   ts = ts[gids]

   # Now mangle the spectra, deredden and integrate-em
   filters_used = []
   boloflux = []
   epochs = []
   mfuncs = []
   fluxes = []
   waves = []
   parss = []
   refbands = []

   mids = []
   for i in range(len(ts)):
      t = ts[i]
      wave,flux = fSED(t/s)
      # Check limits of integration (in rest frame of SN)
      if min(lam1) < wave.min() or (max(lam2) > wave.max() and not extrap_red):
         raise RuntimeError("Error: your limits of integration (%.3f,%.3f) "\
                             "are outside the limits of the SED (%.3f,%.3f)" %\
                             (min(lam1),max(lam2),wave.min(), wave.max()))

      # integration limits
      i1 = [searchsorted(wave, lam) for lam in lam1]
      i2 = [searchsorted(wave, lam) for lam in lam2]
      bs = [bands[j] for j in range(masks.shape[1]) if masks[i,j]]

      if refband is None:
         idx = bands.index(bs[0])
         filt = fset[bs[0]]
         if verbose:  log("   Using %s as reference filter" % bs[0])
      else:
         if refband not in bs:
            if verbose: log("Warning: refband %s has no observation or "
                  "interpolation for epoch %f" % ts[i])
            mids.append(False)
            continue 
         idx = bands.index(refband)
         filt = fset[refband]
      mids.append(True)
      refbands.append(filt.name)

      if len(bs) == 1:
         # No mangling possible
         mflux = flux
         mfunc = (flux*0+1)
         #mfuncs.append(mfunc/mfunc.max())
      else:
         init = [pars0.get(b, 1.0) for b in bs]
         mflux,ave_wave,pars = mangle_spectrum2(wave*(1+sn.z), flux, bs, 
               mags[i,masks[i]], normfilter=refband, init=init, **mopts)
         mfunc = (mflux[0]/flux)
         #mfuncs.append(mfunc/mfunc.max())
         mflux = mflux[0]
         for k,b in enumerate(bs):
            pars0[b] = pars[k]
         parss.append(pars)

      # Scale to match photometry. We need to be careful here. The magnitude
      # measures the response of the filter to the *redshifted* spectrum
      mag = mags[i, idx]
      mflux = mflux*power(10, 
            -0.4*(mag - filt.zp))/filt.response(wave, mflux/(1+sn.z), z=sn.z)
      # Note:  the quantity power()/filt.response() is actually
      # dimensionless. Therefore mflux is in erg/s/cm^2/AA
      # and *not* in photons

      # Next, de-redden MW extinction and host extinction
      mflux,a,b = deredden.unred(wave*(1+sn.z),mflux,sn.EBVgal,R_V=3.1,redlaw=redlaw)
      mflux,a,b = deredden.unred(wave,mflux,EBVhost, R_V=Rv, redlaw=redlaw)

      # Finally!  integrate!
      fbol = []
      ws = []
      fs = []
      mfs = []
      for j in range(len(i1)):
         ws.append(wave[i1[j]:i2[j]])
         fs.append(mflux[i1[j]:i2[j]])
         mfs.append(mfunc[i1[j]:i2[j]])
         fbol.append(trapz(mflux[i1[j]:i2[j]], x=wave[i1[j]:i2[j]]))
         if lam2[j] > wave.max():
            # add Rayleigh-Jeans extrapolation (~ 1/lam^4)
            fbol[-1] += mflux[-1]*wave[-1]/3*(1 - power(wave[-1]/lam2[j],3))

      filters_used.append(bs)
      if scalar:
         boloflux.append(fbol[0])
         waves.append(ws[0])
         fluxes.append(fs[0])
         mfuncs.append(mfs[0])
      else:
         boloflux.append(fbol)
         waves.append(ws)
         fluxes.append(fs)
         mfuncs.append(mfs)
      epochs.append(ts[i])

   mids = array(mids)
   mags = mags[mids,:]
   masks = masks[mids, :]
   boloflux = array(boloflux)
   epochs = array(epochs)
   waves = array(waves)
   fluxes = array(fluxes)
   mfuncs = array(mfuncs)

   # lastly, inverse-square law
   if DM is None:
      DM = sn.get_distmod(cosmo=cosmo)
   dlum = power(10, 0.2*(DM+5))*3.086e18
   boloflux = boloflux*4*pi*dlum**2

   return(dict(epochs=array(epochs), 
               boloflux=array(boloflux), 
               filters_used=filters_used,waves=waves,fluxes=fluxes,
               mfuncs=mfuncs, mags=mags, masks=masks, pars=parss,
               refbands=refbands))

def bolometric_direct(sn, bands=None, tmin=None, tmax=None,
              EBVhost=None, Rv=None, redlaw=None, extrap_red='RJ',
              interpolate=None, interp_all=False, extrapolate=False, 
              SED=None, Tmax=None, DM=None, cosmo='LambdaCDM', verbose=False,
              extra_output=False):

   if verbose: log("Starting bolometric calculation for %s\n" % sn.name)

   if EBVhost is None:
      EBVhost = getattr(sn, 'EBVhost', None)
      if EBVhost is None:
         raise AttributeError("Error:  you must either specify E(B-V) or "\
                               "fit a model that has EBVhost as a parameter")
   if Rv is None:
      Rv = getattr(sn, 'Rv', None)
      if Rv is None:
         raise AttributeError("Error:  you must either specify Rv or "\
                               "fit a model that has Rv as a parameter")

   if redlaw is None:
      redlaw = getattr(sn, 'redlaw', None)
      if redlaw is None:
         raise AttributeError("Error:  you must either specify redlaw or "\
                               "fit a model that has redlaw is defined")

   # We need a time of maximum to set the scale of the Hsiao templates
   if Tmax is None:
      Tmax = getattr(sn, 'Tmax', None)

   if SED is None:
      if verbose: log("   Using Vega SED to compute effective wavelengths")
      fSED = lambda x: (Vega.wave,Vega.flux)
   elif type(SED) is type(""):
      if SED in ['H3','H','N','91bg']:
         fSED = lambda x: get_SED(x, version=SED, extrapolate=True)
         if verbose: log("   Using SED template '%s'" % SED)
      # Assume it is a spectrum by name
      elif SED in standards:
         fSED = lambda x: (standards[SED].wave,standards[SED].flux)
         if verbose: log("   Using standards['%s'] to compute effective wavelengths" % (SED))
      else:
         raise KeyError("SED '%s' not found in standards database" % SED)
   elif type(SED) in [list,tuple]:
      if len(SED) != 2:
         raise ValueError("SED must be tuple or list of length 2")
      fSED = lambda x: SED
   elif type(SED) is types.FunctionType:
      try:
         w,f = SED(0)
      except:
         raise ValueError("If SED is a function, it must take single" \
               " argument (epoch) and return (wave,flux) tuple")
      fSED = SED
   else:
      raise ValueError("Unrecognized type (%s) for SED" % (type(SED)))

   if bands is None:
      bands = getattr(sn.model, '_fbands', None)
      if bands is None:
         bands = list(sn.data.keys())
   for b in bands:
      if b not in sn.data:
         raise AttributeError("band %s not defined in data set" % (b))

   # Bands must be increasing in wavelength
   eff_waves = array([fset[b].eff_wave(Vega) for b in bands])
   sids = argsort(eff_waves)
   bands = [bands[i] for i in sids]

   # Now check that we can interpolate if needed
   if interpolate is not None:
      if interpolate == 'spline':
         if verbose: log("   Using spline interpolation")
         for b in bands:
            if getattr(sn.data[b], 'interp', None) is None:
               raise ValueError("You asked for spline interpolation, but "\
                     "filter %f has not interpolator defined" % b)
      else:
         if verbose: log("   Using model interpolation")
         for b in bands:
            if b not in sn.model._fbands:
               raise ValueError("You asked for model interpolation, but "\
                     "filter %f was not fit with the model" % b)
   else:
      if verbose: log("   Not using interpolation")

   # Now, we build up the times at which we will be integrating and the 
   # assocuated fluxes. We do things in the frame of the SN
   res = sn.get_mag_table(bands)
   if Tmax is None:
      ts = (res['MJD'] - res['MJD'][0])/(1+sn.z)
   else:
      ts = (res['MJD'] - Tmax)/(1+sn.z)
   mags = []
   masks = []
   for b in bands:
      mags.append(res[b])
      masks.append(less(res[b],90))

   if interpolate == 'spline':
      # we fill in (where we can) missing data using splines
      for i,b in enumerate(bands):
         # interpolation is in the absolute time, so use MJD
         mag,mask = sn.data[b].interp(res['MJD'])
         if interp_all:
            mags[i] = mag
         else:
            mags[i] = where(masks[i], mags[i], mag)
         if extrapolate:
            masks[i] = -isnan(mags[i])
         else:
            masks[i] = masks[i] + mask

   elif interpolate == 'model':
      # we fill in (where we can) missing data using the model
      for i,b in enumerate(bands):
         mag,emag,mask = sn.model(b, res['MJD'], extrap=extrapolate)
         mags[i] = where(masks[i], mags[i], mag)
         masks[i] = masks[i] + mask

   mags = transpose(array(mags))
   masks = transpose(array(masks))

   # Now we de-redden, host in the rest frame of the SN, MW in observed
   # frame
   if verbose:  log("   Working on data matrix with size (%d,%d)" % mags.shape)
   # see if any days have zero data
   gids = greater(sum(masks, axis=1), 0)
   mags = mags[gids,:]
   masks = masks[gids,:]
   ts = ts[gids]

   # Restrict on wanted interval if necessary
   if tmin is not None and tmin > ts.min():
      gids = greater(ts, tmin)
   if tmax is not None and tmax < ts.max():
      gids = less(ts, tmax)
   ts = ts[gids]

   # Now mangle the spectra, deredden and integrate-em
   filters_used = []
   boloflux = []
   epochs = []
   fluxes = []
   lam_effs = []
   for i in range(len(ts)):
      wave,flux = fSED(ts[i])

      # Deredden
      A_lamb_host = array([fset[f].R(Rv=Rv, wave=wave, flux=flux, z=sn.z) \
            for f in bands])*EBVhost
      A_lamb_MW = array([fset[f].R(Rv=3.1, wave=wave, flux=flux) \
            for f in bands])*sn.EBVgal
      mags[i] = mags[i] - A_lamb_host - A_lamb_MW
      bs = [bands[j] for j in range(masks.shape[1]) if masks[i,j]]
      zp = array([fset[f].zp for f in bands])

      # Now assume f_\lambda = A*f(\lambda), where f() is dimensionless.
      # A = 10^(-0.4*(m-zpZ))/filter.reponse(f)
      flam = power(10, -0.4*(mags[i]-zp)) # In photons/s/cm^2
      flam = flam/array([fset[f].response(wave,flux, z=sn.z) for f in bands])

      # need to evaluate SED at \lambda_{eff}
      lam_eff = array([fset[f].eff_wave(wave,flux, z=sn.z) for f in bands])

      # Now, we do a weighted average of flam over the response of the filter
      flam = flam*\
            array([fset[f].response(wave,flux,z=sn.z,photons=0)/\
                   fset[f].response(wave,wave*0+1, z=sn.z, photons=0) \
                   for f in bands])

      lam_effs.append(lam_eff[masks[i]])
      fluxes.append(flam[masks[i]])

      # Finally!  integrate!
      fbol = trapz(flam[masks[i]], x=lam_eff[masks[i]])
      if not len(bs) > 1: fbol = flam[masks[i]][0]

      filters_used.append(bs)
      boloflux.append(fbol)
      epochs.append(ts[i])

   boloflux = array(boloflux)
   epochs = array(epochs)

   # lastly, inverse-square law
   if DM is None:
      DM = sn.get_distmod(cosmo=cosmo)
   dlum = power(10, 0.2*(DM+5))*3.086e18
   boloflux = boloflux*4*pi*dlum**2

   return(dict(epochs=array(epochs), 
               boloflux=array(boloflux), 
               filters_used=filters_used,
               fluxes=fluxes, lam_effs=lam_effs,
               mags=mags, masks=masks))
