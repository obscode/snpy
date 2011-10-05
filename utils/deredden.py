## Automatically adapted for numpy.oldnumeric Feb 04, 2009 by ipython

import numpy.oldnumeric as num

def ccm(wave, strict_ccm=0):
   '''Returns the Cardelli, Clayton, and Mathis (CCM) reddening curve.'''
   x = 10000./ wave                ; #Convert to inverse microns 
   a = 0.0*x
   b = 0.0*x
  
   # **************************************************************
   good = num.greater(x, 0.3) * num.less(x, 1.1)   # Infrared
   if len(num.nonzero(good)) > 0: 
      a = num.where(good, 0.574 * num.power(x, 1.61), a)
      b = num.where(good, -0.527 * num.power(x, 1.61), b)
  
   #****************************************************************
   good = num.greater_equal(x,1.1) * num.less(x, 3.3)        #Optical/NIR
   if len(num.nonzero(good)) > 0:      #Use new constants from O'Donnell (1994)
      y = x - 1.82
      if strict_ccm:
         c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085,     #Original
                       0.01979, -0.77530,  0.32999 ]            #coefficients
         c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434,     #from CCM89
                      -0.62251,  5.30260, -2.09002 ]
      else:
         # New coefficents from O'Donnell (1994)
         c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137,
                      -1.718,   -0.827,    1.647, -0.505 ]
         c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985,  
                      11.102,    5.491,  -10.805,  3.347 ]
  
      poly = 0.0
      for i in range(len(c1)): poly = poly + c1[i]*num.power(y, i)
      a = num.where(good, poly, a)
      poly = 0.0
      for i in range(len(c2)): poly = poly + c2[i]*num.power(y, i)
      b = num.where(good, poly, b)

   #******************************************************************
   good = num.greater_equal(x, 3.3) * num.less(x,8)          #Mid-UV
   if len(num.nonzero(good)) > 0:
      y = 1.0*x
      good1 = num.greater(y, 5.9)
      F_a = y*0.0    ; F_b = y*0.0
      if len(num.nonzero(good1)) > 0:
         y1 = y - 5.9
         F_a = -0.04473 * num.power(y1,2) - 0.009779 * num.power(y1,3)
         F_b = 0.2130 * num.power(y1,2)  +  0.1207 * num.power(y1,3)
         F_a = num.where(good1, F_a, 0.0)
         F_b = num.where(good1, F_b, 0.0)
      
      a = num.where(good, 
           1.752 - 0.316*y - (0.104 / ( num.power(y-4.67,2) + 0.341 )) + F_a, a)
      b = num.where(good,
           -3.090 + 1.825*y + (1.206 / ( num.power(y-4.67,2) + 0.263 )) + F_b, b)
  
   #   *******************************
  
   good = num.greater_equal(x, 8) * num.less_equal(x, 11)  #Far-UV
   if len(num.nonzero(good)) > 0:
      y = x - 8.0
      c1 = [ -1.073, -0.628,  0.137, -0.070 ]
      c2 = [ 13.670,  4.257, -0.420,  0.374 ]
      poly = 0.0*y
      for i in range(len(c1)):  poly = poly + c1[i]*num.power(y,i)
      a = num.where(good, poly, a)
      poly = 0.0*y
      for i in range(len(c2)):  poly = poly + c2[i]*num.power(y,i)
      b = num.where(good, poly, b)

   return(a,b)


def unred(wave, flux, ebv, R_V = 3.1, z=0, strict_ccm=0):
   '''
      pro ccm_UNRED, wave, flux, ebv, funred, R_V = r_v
   
   NAME:
       CCM_UNRED
   PURPOSE:
       Deredden a flux vector using the CCM 1989 parameterization 
   EXPLANATION:
       The reddening curve is that of Cardelli, Clayton, and Mathis (1989 ApJ.
       345, 245), including the update for the near-UV given by O'Donnell 
       (1994, ApJ, 422, 158).   Parameterization is valid from the IR to the 
       far-UV (3.5 microns to 0.1 microns).    
   
       Users might wish to consider using the alternate procedure FM_UNRED
       which uses the extinction curve of Fitzpatrick (1999).
   CALLING SEQUENCE:
       CCM_UNRED, wave, flux, ebv, funred, [ R_V = ]      
               or 
       CCM_UNRED, wave, flux, ebv, [ R_V = ]      
   INPUT:
       WAVE - wavelength vector (Angstroms)
       FLUX - calibrated flux vector, same number of elements as WAVE
               If only 3 parameters are supplied, then this vector will
               updated on output to contain the dereddened flux.
       EBV  - color excess E(B-V), scalar.  If a negative EBV is supplied,
               then fluxes will be reddened rather than deredenned.
   
   OUTPUT:
       FUNRED - unreddened flux vector, same units and number of elements
               as FLUX
   
   OPTIONAL INPUT KEYWORD
       R_V - scalar specifying the ratio of total selective extinction
               R(V) = A(V) / E(B - V).    If not specified, then R_V = 3.1
               Extreme values of R(V) range from 2.75 to 5.3
   
   EXAMPLE:
       Determine how a flat spectrum (in wavelength) between 1200 A and 3200 A
       is altered by a reddening of E(B-V) = 0.1.   Assume an "average"
       reddening for the diffuse interstellar medium (R(V) = 3.1)
   
         IDL> w = 1200 + findgen(40)*50      ;Create a wavelength vector
         IDL> f = w*0 + 1                    ;Create a "flat" flux vector
         IDL> ccm_unred, w, f, -0.1, fnew  ;Redden (negative E(B-V)) flux vector
         IDL> plot,w,fnew                   
   
   NOTES:
       (1) The CCM curve shows good agreement with the Savage & Mathis (1979)
               ultraviolet curve shortward of 1400 A, but is probably
               preferable between 1200 and 1400 A.
       (2)  Many sightlines with peculiar ultraviolet interstellar extinction 
               can be represented with a CCM curve, if the proper value of 
               R(V) is supplied.
       (3)  Curve is extrapolated between 912 and 1000 A as suggested by
               Longo et al. (1989, ApJ, 339,474)
       (4) Use the 4 parameter calling sequence if you wish to save the 
                 original flux vector.
       (5) Valencic et al. (2004, ApJ, 616, 912) revise the ultraviolet CCM
               curve (3.3 -- 8.0 um-1).    But since their revised curve does
               not connect smoothly with longer and shorter wavelengths, it is
               not included here.
   
   REVISION HISTORY:
         Written   W. Landsman        Hughes/STX   January, 1992
         Extrapolate curve for wavelengths between 900 and 1000 A   Dec. 1993
         Use updated coefficients for near-UV from O'Donnell   Feb 1994
         Allow 3 parameter calling sequence      April 1998
         Converted to IDLV5.0                    April 1998

   Added a redshift argument so that we can (de)redden in a different rest frame.  It is used
   in the sense that the wavelengths are redshifted by z before the calculation is
   performed.  We do this because the SED templates are defined in the observed frame.
   '''

   if z > 0:
      wave = wave*(1+z)
  
   a,b = ccm(wave, strict_ccm)
   #; Now apply extinction correction to input flux vector
  
   #A_V = R_V * ebv
   A_lambda = ebv * (a*R_V + b)
   unred_flux = flux * num.power(10.0, 0.4*A_lambda)
   return(unred_flux, a , b)


def R_z(wave, z, R_V = 3.1, strict_ccm=0):
   '''Compute host reddening law for effective observed wavelength wave, at redshift
   z, assuming R_V.'''
   wave = wave/(1+z)
   a,b = ccm(wave, strict_ccm)
   R = R_V*a + b
   return(R, a , b)
