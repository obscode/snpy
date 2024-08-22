#!/usr/bin/env python
'''
A module for computing k-corrections. This code has been assembled from many
sources, mostly from IRAF scripts written by Mark Phillips and IDL code 
written by Mark Sullivan.

Most of the heavy lifting w.r.t. integration of filters over the SED has
been moved to the :mod:`snpy.filter` module.
'''

from __future__ import print_function
import os,sys,string,re
import numpy as num
import scipy.interpolate
import h5py
import onnxruntime as onnx_rt
from scipy.integrate import simpson as simps
from .utils import deredden
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
from . import filters
from .mangle_spectrum import mangle_spectrum2, default_method

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
# read in building blocks for Lu2023 NIR template:
NIR_temp_path = os.path.join(spec_base, 'NIR_template_buildingblocks.h5')
NIR_temp_dict = {}
# Open the HDF5 file in read mode and load the onnx model too 
with h5py.File(NIR_temp_path, "r") as file:
    # Iterate over the dataset keys and values
    for key, value in file.items():
        # check if they are a group, sorry for so many layers of nested groups ><  
        if isinstance(value, h5py.Group):
            NIR_temp_dict[key] = {}
            for subkey, subvalue in value.items():
                if isinstance(subvalue, h5py.Group):
                    NIR_temp_dict[key][subkey] = {}
                    for subsubkey, subsubvalue in subvalue.items():
                        if isinstance(subsubvalue, h5py.Group):
                            NIR_temp_dict[key][subkey][subsubkey] = {}
                            for subsubsubkey, subsubsubvalue in subsubvalue.items():
                                if subsubsubkey == "GPR_onx_string": # it's byte datatype
                                    NIR_temp_dict[key][subkey][subsubkey][subsubsubkey] = num.bytes_(subsubsubvalue[()])
                                    # load the model here is better than in the function - since this is the slow part 
                                    NIR_temp_dict[key][subkey][subsubkey]["GPR_onx"] = onnx_rt.InferenceSession(
                                                                                          subsubsubvalue[()], providers=["CPUExecutionProvider"]
                                                                                    )
                                else:
                                    NIR_temp_dict[key][subkey][subsubkey][subsubsubkey] = subsubsubvalue[()]
                        else:
                            NIR_temp_dict[key][subkey][subsubkey] = subsubvalue[()]
                else:
                    NIR_temp_dict[key][subkey] = subvalue[()]
        else:
            NIR_temp_dict[key] = value


def linterp(spec1, spec2, day1, day2, day):
   if day1 == day2:
      return spec1
   if day > day2 or day < day1:
      raise ValueError("day must be in interval (day1,day2)")
   if abs(day1-day) < 1e-9:
      return spec1
   if abs(day2 - day) < 1e-9:
      return spec2
   spec = spec1 + (spec2 - spec1)/(day2 - day1)*(day - day1)
   return spec

SED_lims = {
      'H':(-19,70),
      'H3':(-19,70),
      'H3+L':(-19,70),
      'N':(-19,70),
      '91bg':(-13,100)}

sBV_lims = {
      'H3+L':(0.4,1.3)}


''' Functions needed for the new CSP2  NIR templates '''
def GPR_predict_PC(Wbin,epoch,sBV):
    '''Function of predicting Principle components from the fitted 
    Gaussian Process Regression (GPR) model for the NIR templates

    Args:
        Wbin (str): 'W1' to 'W7', correspinding to z,Y,J,telluric1,H,telluric2,K band
        epoch (float): restframe days since B-band maximum of the spectrum
        sBV (float): color stretch of the SN

    Returns:
        2-tuple: (pred_PCs,pred_PC_sigmas)

        * pred_PCs: (float array) Gaussian Process predicted pricinple compoennts
        * pred_PC_sigmas: (float array) uncertainty ofGaussian Process
                          predicted pricinple compoennts (CURRENLTY NOT IMPLEMENTED)
    '''
    #TODO? Try to include the uncertainty prediction in saved onnx model 
    #(had issue with white kernel diagonal matrix inversion)

    ## initialize the list to collect the predicted PCs
    pred_PCs = []
    for i in range(len(NIR_temp_dict["GPR"][Wbin].keys())):
        PC = 'PC'+str(i+1)
        ## locate the original PC range and min to get the normalization factor for later
        yrange = NIR_temp_dict["GPR"][Wbin][PC]["yrange"]
        ymin = NIR_temp_dict["GPR"][Wbin][PC]['ymin']
        ## locate the fitted GPR onx_string and load the onnx model
        GPR_onx = NIR_temp_dict["GPR"][Wbin][PC]["GPR_onx"] 
        y_gp = GPR_onx.run(None, {"X": num.array([[num.float64(epoch),num.float64(sBV)]])})[0]
        ## unnormalize the predicted y using the yrange and min above
        y_gp_unnorm = y_gp * yrange + ymin #y_sigma_unnorm = y_sigma*yrange
        pred_PCs.append(y_gp_unnorm.flatten()[0])
    return pred_PCs,None

def PCA_inverse_transform(Wbin, PC_values,PC_masks):
   '''Function of inverse transforming the principle components to the flux space 

   Args:
      Wbin (str): 'W1' to 'W7', correspinding to z,Y,J,telluric1,H,telluric2,K band
      PC_values (array): the predicted principle component projection values

   Returns:
      PC_inversed (array): the inverse transformed flux from the principle components
   '''
   PC_components = NIR_temp_dict["PCA_components"][Wbin]["components"]
   PC_masked_values = PC_values*PC_masks
   PC_inversed = num.sum(PC_masked_values.reshape(-1,1) * PC_components,axis=0) \
                 + NIR_temp_dict["PCA_components"][Wbin]["pca_mean"]
   return PC_inversed


def merge_spec(wave_b,flux_b,wave_r,flux_r,interp_option=0,normalize='blue_side',plot_ax=None):
    '''Function of merging spectra, *_b as bluer wavelength side, *_r as
       redder wavelength side

    Args:
        wave_b (float array): the bluer wavelength side of the spectra wavelength
                              need to be merged
        flux_b (float array): the bluer wavelength side of the spectra flux
                              need to be merged
        wave_r (float array): the redder wavelength side of the spectra wavelength
                              need to be flux_merged
        flux_r (float array): the redder wavelength side of the spectra flux
                              need to be flux_merged

        interp_option (int): the option of how the wavelength sampling point
                             in the overlap region is dealed with, the options can be:
            * 0: take both wavelength sampling points from blue and red side
            * 1: only take the wavelength sampling points from blue side
            * 2: only take the wavelength sampling points from red side

        normalize (str): how the flux is normalized for two side of the spectra,
                         the options can be:
            * 'blue_side': scale the red side flux by matching the overlapped
                           region flux to the overlapped region flux of blue side
            * 'red_side': scale the blue side flux by matching the overlapped
                           region flux to the overlapped region flux of red side
            * else: no scaling, just merge smoothly by weights

        plot_ax: default no plot outputs (None), can also given ax to visualize
                 the merge process


    Returns:
        2-tuple: (wave_merged,flux_merged)

        * wave_merged: (float array) merged wavelength
        * flux_merged: (float array) merged flux
    '''
    ## converted from Eric's IDL code EYH_MERGE_SPEC
    weight_lo=0.
    w1 = num.where(wave_b>=min(wave_r))[0]
    w2 = num.where(wave_r<=max(wave_b))[0]
    if (len(w1)<2) or (len(w2)<2):
        print('WARNING! Not enough overlap')
        return
    else:
        ### decide what wavelegnth to use in the overlap region
        if interp_option ==0: ## combine the overlap wavelengths
            wave_overlap = num.concatenate([wave_b[w1], wave_r[w2]],axis=0)
            wave_overlap = num.sort(wave_overlap)
            wave_overlap = num.unique(wave_overlap)
        elif interp_option ==1:
            wave_overlap = wave_b[w1]
        elif interp_option ==2:
            wave_overlap = wave_r[w2]
        ### decide how to normalize
        if normalize == 'blue_side':
            ### match the overlappen region flux based on the blue side
            norm_b = 1
            norm_r = simps(flux_b[w1],x=wave_b[w1])/\
                     simps(flux_r[w2],x=wave_r[w2])
        elif normalize == 'red_side':
            norm_b = simps(flux_r[w2],x=wave_r[w2])/\
                  simps(flux_b[w1],x=wave_b[w1])
            norm_r = 1
        else:
            norm_b,norm_r = 1,1 ## let the spectra merge by weight, smooothly merging together

        ### inteplate the flux in overlaped region
        x = [min(wave_overlap), max(wave_overlap)]
        f1 = scipy.interpolate.interp1d(x, [1.,weight_lo])
        f2 = scipy.interpolate.interp1d(x, [weight_lo,1.])
        weight1 = f1(wave_overlap)
        weight2 = f2(wave_overlap)
        f1_f = scipy.interpolate.interp1d(wave_b[w1],norm_b*flux_b[w1],fill_value="extrapolate")
        f2_f = scipy.interpolate.interp1d(wave_r[w2],norm_r*flux_r[w2],fill_value="extrapolate")
        flux_ol = 1./(weight1+weight2)*(f1_f(wave_overlap)*weight1 + f2_f(wave_overlap)*weight2)

        ### now combine the flux
        w1_rest = num.where(wave_b<min(wave_r))[0]
        w2_rest = num.where(wave_r>max(wave_b))[0]
        wave_out = num.concatenate([wave_b[w1_rest],wave_overlap,wave_r[w2_rest]],axis=0)
        flux_out = num.concatenate([norm_b*flux_b[w1_rest],flux_ol,norm_r*flux_r[w2_rest]],axis=0)

        if plot_ax is not None:
            ax.plot(wave_b,norm_b*flux_b,'b',alpha=0.4)
            ax.plot(wave_r,norm_r*flux_r,'r',alpha=0.4)
            ax.plot(wave_out,flux_out,'k',alpha=0.2)
        return wave_out,flux_out

def get_single_NIR_template(epoch,sBV):
    '''function of getting one template spectra based on given epoch and sBV

    Args:
        epoch (float): restframe days since B-band maximum of the spectrum
        sBV (float): color stretch of the SN

    Returns:
        tuple:

        2 tuple: (wavelength, flux)
            * wavelength: (float array) NIR template wavelength in unit of Angstroms
            * flux: (float array)  NIR template flux in arbitrary unit
    '''
    ## define some parameters
    GPR_score_threshold=0.2   ## lower limit of GPR score when selecting PC
    ## check data type and of value is in range
    if isinstance(epoch,list) or isinstance(epoch,num.ndarray) or isinstance(sBV,list) or isinstance(sBV,num.ndarray):
        raise ValueError('Input single value for epoch and sBV please, one at a time :)')
    else:
        if (epoch < SED_lims['H3+L'][0]) or (epoch > SED_lims['H3+L'][1]):
            raise ValueError('Epoch not in supported range, please input epoch \
            between %d to %d'%(SED_lims['H3+L'][0],SED_lims['H3+L'][1]))
        if (sBV < sBV_lims['H3+L'][0]) or (sBV > sBV_lims['H3+L'][1]):
            raise ValueError('sBV not in supported range, please input sBV \
            between %.1f to %.1f'%(sBV_lims['H3+L'][0],sBV_lims['H3+L'][1]))
    x1,x2 = epoch,sBV
    PC_names = ['PC'+str(i+1) for i in range(len(NIR_temp_dict["GPR"]["W1"].keys()))]
    for i in range(len(NIR_temp_dict["PCA_components"].keys())-1):
        ## read the pca transformation and wave grid in the adjacent blocks
        if i==0: # start with the bluest block as the first "merged product"
            Wbin_blue = 'W1'
            wave_merged = NIR_temp_dict["PCA_components"][Wbin_blue]["wavelength"] 
            # get the predicted PC values from GPR
            PC_values_blue,_ = GPR_predict_PC(Wbin_blue,x1,x2)
            # mask the PCs as 0 if the GPR for this PC is under threshold
            if GPR_score_threshold is not None:
               gp_scores = [NIR_temp_dict["GPR"][Wbin_blue][PC_key]["gp_score"] for PC_key in PC_names]
               PC_masks_blue = num.array([1. if (score >= GPR_score_threshold) else 0. for score in gp_scores])
            else:
               PC_masks_blue = num.ones(len(gp_scores))
            flux_merged = PCA_inverse_transform(Wbin_blue, PC_values_blue,PC_masks_blue)
        # attach the next wavelength block to previous merged product
        Wbin_red = 'W'+str(i+2) 
        wave_red = NIR_temp_dict["PCA_components"][Wbin_red]["wavelength"] 
        PC_values_red,_ = GPR_predict_PC(Wbin_red,x1,x2)
        if GPR_score_threshold is not None:
           gp_scores = [NIR_temp_dict["GPR"][Wbin_red][PC_key]["gp_score"] for PC_key in PC_names]
           PC_masks_red = num.array([1. if (score >= GPR_score_threshold) else 0. for score in gp_scores])
        else:
           PC_masks_red = num.ones(len(NIR_temp_dict["GPR"][Wbin_red].keys()))
        flux_red = PCA_inverse_transform(Wbin_red, PC_values_red,PC_masks_red)
        ## merge the flux from the W bins
        wave_merged,flux_merged = merge_spec(wave_merged,flux_merged,wave_red,flux_red)

    return wave_merged*10000,flux_merged

def get_SED_H3_plus_L(epoch,sBV,extrapolate=False):
    '''Function to substitude the NIR part of  the Hsiao template (H3) with
       the new CSP2 NIR template (L)

    args:
        epoch (float): restframe days since B-band maximum of the spectrum
        sBV (float): color stretch of the SN

    Returns:
        2-tuple: (wave,flux)

        * wave (array):  Wavelength in Angstroms
        * flux (array):  arbitrarily normalized flux
    '''
    # get the NIR template
    wave_NIR,flux_NIR = get_single_NIR_template(epoch,sBV)
    # get the H3 templates
    stretched_epoch = epoch/sBV
    # Check limits
    if stretched_epoch <= SED_lims['H3'][0]:
       if extrapolate:
          day = SED_lims['H3'][0]
       else:
          return (None,None)
    elif stretched_epoch > SED_lims['H3'][1]:
       if extrapolate:
          day = SED_lims['H3'][1]
       else:
          return (None,None)
    else:
        day = stretched_epoch
    wave_hsiao,flux_hsiao = get_SED(day, version='H3')

    ## substitude the NIR region of H3
    w_optical = num.where(wave_hsiao<=wave_NIR[0]+300) # set overlap region to be 200A
    wave_optical,flux_optical = wave_hsiao[w_optical],flux_hsiao[w_optical]
    # merge the optical with NIR
    wave_merged,flux_merged = merge_spec(wave_optical,flux_optical,wave_NIR,flux_NIR)
    # merge with the H3 NIR tail since L NIR template only goes to 2.33um
    w_NIRtail = num.where(wave_hsiao>=wave_NIR[-1]-300)
    wave_tail,flux_tail = wave_hsiao[w_NIRtail],flux_hsiao[w_NIRtail]
    wave_merged,flux_merged = merge_spec(wave_merged,flux_merged,wave_tail,flux_tail)

    return wave_merged,flux_merged
''' End of Functions needed for the new CSP2  NIR templates '''


def get_SED(day, version='H3',sBV=1.0, interpolate=True, extrapolate=False):
   '''Retrieve the SED for a SN for a particular epoch.
   
   Args:
      day (int or float): The integer day w.r.t. time of B-maximum
      sBV (float): The color stretch of the SN, only appliable to version='H3+L'
      version (str): The version of SED sequence to use:

         * 'H': Old Hsiao Ia SED (Hsiao, private communication)
         * 'H3': Hsiao+2007 Ia SED
         * 'H3+L': Hsiao+2007 Ia SED + Lu+2022 NIR SED
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
   epoch = day

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
   elif version =='H3+L':
       return (get_SED_H3_plus_L(epoch,sBV,extrapolate=extrapolate))
   else:
      raise AttributeError("version %s not recognized" % version)

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

def S(wave, spec, f1, f2, z):
   '''compute single S-correction based on a single spectrum and set of 
   observed filters. This is like a K-correction, except we are transforming
   from observer frame filters to observer frame filters, so the only 
   redshifting involves is in the underlying SED.  The output S-correction
   is in the sense that m(f2) = m(f1) + S

   Args:
      wave (float array): input wavelength in Angstroms
      flux (float array): arbitrarily scaled flux
      f1 (filter instance): Observed source filter.
      f2 (filter instance): Observed destination filter.
      z (float): redshift
      photons (bool): If True, fluxes are computed in units of photons rather
                      than energy (see Nugent+2002)
   
   Returns:
      2-tuple: (S,flag)

      * S: S-correction
      * flag: 1 -> success, 0->failed
   '''

   # compute the magnitude through each filter
   mag1 = f1.synth_mag(wave, spec, z=z)
   mag2 = f2.synth_mag(wave, spec, z=z)

   if num.isnan(mag1) or num.isnan(mag2):
      # Something clearly went wrong
      return (0.0, 0)
   else:
      S = mag2 - mag1 
      return (S, 1)

def kcorr(days, filter1, filter2, z, sBV=1.0, ebv_gal=0, ebv_host=0, R_gal=3.1, 
      R_host=3.1, version="H3", photons=1, Scorr=False, extrapolate=False):
   '''Find the cross-band k-correction for a series of type Ia SED from
   SNooPy's catalog. These can be thought of as "empirical" K-corrections.
   
   Args:
      days (float array): epochs (t-T(Bmax)) to compute
      filter1 (str):  rest-frame filter
      filter2 (str):  observed filter. This can be the same as filter1,
                      or another, redder, filter for cross-band K-corrections
      z (float): redshift
      sBV (float): The color stretch of the SN, only appliable to version='H3+L'
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
      Scorr (bool):  If True, return an S-correction rather than a K-correction.

   Returns
      2-tuple:  (K,mask)

      * K (float array): K-corrections
      * mask (bool array): True where K-corrections are valid.
                     
   '''

   if filter1 not in filters.fset:
      raise AttributeError("filter %s not defined in filters module" % filter1)
   if filter2 not in filters.fset:
      raise AttributeError("filter %s not defined in filters module" % filter2)

   kcorrs = []
   mask = []      # Masks the good values (1) and bad (not defined) values (0)
   # Loop through the list of days
   for day in days:
      if version != 'H3+L':
          day = int(day)
      spec_wav,spec_f = get_SED(day, version, extrapolate=extrapolate)
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
      if not Scorr:
         k,f = K(spec_wav, sp_f, f1, f2, z)
      else:
         k,f = S(spec_wav, sp_f, f1, f2, z)
      kcorrs.append(k)
      mask.append(f)
   return(kcorrs,mask)

def kcorr_mangle2(waves, spectra, filts, mags, m_mask, restfilts, z, 
      colorfilts=None, full_output=0, Scorr=False, **mopts): 
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
      Scorr (bool):  If True, compute S-corrections rather than K-corrections
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
         raise AttributeError("filter %s not defined in filters module" % filter1)

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
         man_spec_f,state,pars = mangle_spectrum2(spec_wav*(1+z),spec_f,fs, 
               ms, **mopts)
      if full_output:
         waves_a.append(state['ave_waves'])
         manf_a.append(man_spec_f)
         factors_a.append(pars)
 
      for i in range(len(filts)):
         f1 = filters.fset[restfilts[i]]
         f2 = filters.fset[filts[i]]

         if Scorr:
            k,f = S(spec_wav, man_spec_f[0], f1, f2, z)
         else:
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

def kcorr_mangle(days, filts, mags, m_mask, restfilts, z,sBV = 1.0, version='H3', 
      colorfilts=None, full_output=0, mepoch=False, Scorr=False, 
      extrapolate=False, **mopts):
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
      sBV (float): The color stretch of the SN, only appliable to version='H3+L'
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
      Scorr (bool): If True, compute S-corrections rather than K-corrections
      extrapolate (bool): If True, extrapolate beyond limits of SED sequence
                          (see :func:`.get_SED`.
      mopts (dict): All additional arguments to function are sent to 
                    :func:`snpy.mangle_spectrum.mangle_spectrum2`.
   
   Returns:
      tuple:

         * if not full_output: 2-tuple (K,mask):
            * K (flaot array):  K-corrections for filts
            * mask (bool array): mask of valid K-corrections
         * if full_output: 4-tuple (K,mask,Rts,mopts)
            * Rts (float array): The total-to-selective absorption ratio
                                 based on filter functions and mangled SED
            * mopts (dict): mangling function state dictionary
   '''

   if 'method' in mopts:
      method = mopts['method']
   else:
      method = default_method

   if colorfilts is None:  colorfilts = filts
   for filter1 in filts + restfilts + colorfilts:
      if filter1 not in filters.fset:
         raise AttributeError("filter %s not defined in filters module" % filter1)

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
         if version != 'H3+L':
             day = int(days[j])
         else:
             day = days[i]
         s,f = get_SED(day, version, extrapolate=extrapolate)
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
            if num.any(m_mask[:,i])]
      if len(fs) <=1 :
         waves,man_spec_fs,factors = spec_wavs,spec_fs,spec_wavs*0.0+1.0
         state,pars = None,None
      else:
         #cs = mags[:,:-1] - mags[:,1:]
         #gids = m_mask[:,:-1]*m_mask[:,1:]
         #gids = gids*sids[:,num.newaxis]
         #cs[-gids] = 99.9   # flag invalid value
         gids = m_mask*sids[:,num.newaxis]
         ms = where(gids, ms, 99.9)
         man_spec_fs,state,pars = mangle_spectrum2(spec_wavs*(1+z),
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
            args = {}
            if state is not None:
               args['state'] = state
               args['pars'] = pars
            for key in mopts:
               args[key] = mopts[key]
            m_opts.append(args)
 
         for i in range(len(filts)):
             f1 = filters.fset[restfilts[i]]
             f2 = filters.fset[filts[i]]
             if Scorr:
                k,f = S(spec_wavs[j], man_spec_fs[j], f1, f2, z)
             else:
                k,f = K(spec_wavs[j], man_spec_fs[j], f1, f2, z)
             kcorrs[-1].append(k)
             mask[-1].append(0)
         Rts.append(R_obs_spectrum(filts, spec_wavs[j], man_spec_fs[j], z, 
            0.01, 0.0))
   else:
      for j in range(len(days)):
         kcorrs.append([])
         mask.append([])
         if version != 'H3+L':
             day = int(days[j])
         else:
             day = days[j]
         spec_wav,spec_f = get_SED(day, version, extrapolate=extrapolate)
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
            state,pars = None,None
         else:
            #cs = num.compress(m_mask[j],mags[j])[0:-1] - \
            #     num.compress(m_mask[j],mags[j])[1:]
            ms = num.compress(m_mask[j], mags[j])
            if debug:
               print("filters and colors for day %f:" % (days[j]))
               print(fs)
               print(ms[:-1]-ms[1:])
  
            # Now we mangle the spectrum.  Note, we are redshifting the spectrum
            # here, so do NOT set z in mangle_spectrum2.
            man_spec_f,state,pars = mangle_spectrum2(spec_wav*(1+z),spec_f,
                  fs, ms, **mopts)
 
            if debug:  print("factors = ",factors)
            if debug:
               # check the colors
               for i in range(len(fs)-1):
 
                  print("input color:  %s-%s = %f" % (fs[i],fs[i+1], ms[i]-ms[i+1]), end=' ')
                  f1 = filters.fset[fs[i]]
                  f2 = filters.fset[fs[i+1]]
                  col = f1.synth_mag(spec_wav*(1+z), man_spec_f[0]) - \
                        f2.synth_mag(spec_wav*(1+z), man_spec_f[0])
                  print("  output color:  %f" % (col))
  
         if full_output:
            args = {}
            if state is not None:
               args['state'] = state
               args['pars'] = pars
            for key in mopts:
               args[key] = mopts[key]
            m_opts.append(args)
  
         for i in range(len(filts)):
            f1 = filters.fset[restfilts[i]]
            f2 = filters.fset[filts[i]]
            if Scorr:
               k,f = S(spec_wav,man_spec_f[0], f1, f2, z)
            else:
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
      Rv_host=3.1, Rv_gal=3.1, version='H', redlaw='f99', strict_ccm=False,
      extrapolate=False):
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
       extrapolate (bool): Extrpolate beyond SED limits. See :func:`.get_SED`.
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
      spec_wav,spec_f = get_SED(int(day), version, extrapolate=extrapolate)
      if spec_wav is None:
         Rs.append(99.9)
      # Redden the spectrum based on Cardelli et al. and assumed EBVgal + EBVhost
      red_f = redden(spec_wav, spec_f, EBVgal, EBVhost, z, Rv_gal, Rv_host,
            redlaw=redlaw, strict_ccm=strict_ccm)
      for filter in [filter1,filter2,filter3]:
         if filter not in filters.fset:
            raise AttributeError("filter %s not defined in filters module" % filter)
 
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
      spec_wav,spec_f = get_SED(int(day), version, extrapolate=extrapolate)
      if spec_wav is None:
         Rs.append(99.9)
         continue
      # rEDDEN the spectrum based on Cardelli et al. and assumed EBVgal+EBVhost
      red_f = redden(spec_wav, spec_f, EBVgal, EBVhost, z, Rv_gal, Rv_host)

      if filter not in filters.fset:
         raise AttributeError("filter %s not defined in filters module" % filter)
 
      # Now, we get the response across the filters:
      resp = filters.fset[filter].response(spec_wav, spec_f, z=z, photons=1)
      resp_red = filters.fset[filter].response(spec_wav, red_f, z=z, photons=1)
 
      As.append(-2.5*num.log10(resp_red/resp))
   if outarr:
      return(num.array(As))
   else:
      return(As[0])

def R_obs(filter, z, days, EBVhost, EBVgal, Rv_host=3.1, Rv_gal=3.1, 
      version='H', redlaw='f99', strict_ccm=False, extrapolate=False):
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
      spec_wav,spec_f = get_SED(int(day), version, extrapolate=extrapolate)
      if spec_wav is None:
         Rs.append(99.9)
         continue
      # rEDDEN the spectrum based on Cardelli et al. and assumed EBVgal + EBVhost
      red_f = redden(spec_wav, spec_f, EBVgal, EBVhost, z, Rv_gal, Rv_host,
            redlaw=redlaw, strict_ccm=strict_ccm)

      if filter not in filters.fset:
         raise AttributeError("filter %s not defined in filters module" % filter)
 
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
         raise AttributeError("filter %s not defined in filters module" % filter)
 
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

