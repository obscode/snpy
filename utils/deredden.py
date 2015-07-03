from scipy.interpolate import interp1d

import numpy as num
from scipy.interpolate import UnivariateSpline

def poly(x, c):
   ret = 0
   for i in range(len(c)):
      ret += c[i]*num.power(x,i)
   return ret

def ccm(wave, strict_ccm=0):
   '''Returns the Cardelli, Clayton, and Mathis (CCM) reddening curve.'''
   x = 10000./ wave                ; #Convert to inverse microns 
   a = 0.0*x
   b = 0.0*x
  
   # **************************************************************
   good = num.greater(x, 0.3) * num.less(x, 1.1)   # Infrared
   if num.sometrue(good):
      a = num.where(good, 0.574 * num.power(x, 1.61), a)
      b = num.where(good, -0.527 * num.power(x, 1.61), b)
  
   #****************************************************************
   good = num.greater_equal(x,1.1) * num.less(x, 3.3)        #Optical/NIR
   if num.sometrue(good):      #Use new constants from O'Donnell (1994)
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
   if num.sometrue(good):
      y = 1.0*x
      good1 = num.greater(y, 5.9)
      F_a = y*0.0    ; F_b = y*0.0
      if num.sometrue(good1) > 0:
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
   if num.sometrue(good):
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

def fm(wave, R_V=3.1, avglmc=False, lmc2=False):
   # NAME:
   #    FM_UNRED
   # PURPOSE:
   #    Deredden a flux vector using the Fitzpatrick (1999) parameterization
   # EXPLANATION:
   #    The R-dependent Galactic extinction curve is that of Fitzpatrick & Massa
   #    (Fitzpatrick, 1999, PASP, 111, 63; astro-ph/9809387 ).    
   #    Parameterization is valid from the IR to the far-UV (3.5 microns to 0.1
   #    microns).    UV extinction curve is extrapolated down to 912 Angstroms.
   #
   # CALLING SEQUENCE:
   #    FM_UNRED, wave, flux, ebv, [ funred, R_V = , /LMC2, /AVGLMC, ExtCurve= 
   #                       gamma =, x0=, c1=, c2=, c3=, c4= ]
   # INPUT:
   #    WAVE - wavelength vector (Angstroms)
   #    FLUX - calibrated flux vector, same number of elements as WAVE
   #             If only 3 parameters are supplied, then this vector will
   #             updated on output to contain the dereddened flux.
   #    EBV  - color excess E(B-V), scalar.  If a negative EBV is supplied,
   #             then fluxes will be reddened rather than dereddened.
   #
   # OUTPUT:
   #    FUNRED - unreddened flux vector, same units and number of elements
   #             as FLUX
   #
   # OPTIONAL INPUT KEYWORDS
   #    R_V - scalar specifying the ratio of total to selective extinction
   #             R(V) = A(V) / E(B - V).    If not specified, then R = 3.1
   #             Extreme values of R(V) range from 2.3 to 5.3
   #
   #    /AVGLMC - if set, then the default fit parameters c1,c2,c3,c4,gamma,x0 
   #           are set to the average values determined for reddening in the 
   #           general Large Magellanic Cloud (LMC) field by Misselt et al. 
   #          (1999, ApJ, 515, 128)
   #    /LMC2 - if set, then the fit parameters are set to the values determined
   #             for the LMC2 field (including 30 Dor) by Misselt et al.
   #             Note that neither /AVGLMC or /LMC2 will alter the default value
   #             of R_V which is poorly known for the LMC. 
   #             
   #   The following five input keyword parameters allow the user to customize
   #   the adopted extinction curve.    For example, see Clayton et al. (2003,
   #   ApJ, 588, 871) for examples of these parameters in different interstellar
   #   environments.
   #
   #      x0 - Centroid of 2200 A bump in microns (default = 4.596)
   #      gamma - Width of 2200 A bump in microns (default  =0.99)
   #      c3 - Strength of the 2200 A bump (default = 3.23)
   #      c4 - FUV curvature (default = 0.41)
   #      c2 - Slope of the linear UV extinction component 
   #           (default = -0.824 + 4.717/R)
   #      c1 - Intercept of the linear UV extinction component 
   #           (default = 2.030 - 3.007*c2
   #            
   # OPTIONAL OUTPUT KEYWORD:
   #    ExtCurve - Returns the E(wave-V)/E(B-V) extinction curve, interpolated
   #                 onto the input wavelength vector
   #
   # EXAMPLE:
   #    Determine how a flat spectrum (in wavelength) between 1200 A and 3200 A
   #    is altered by a reddening of E(B-V) = 0.1.   Assume an "average"
   #    reddening for the diffuse interstellar medium (R(V) = 3.1)
   #
   #    IDL> w = 1200 + findgen(40)*50      ;Create a wavelength vector
   #    IDL> f = w*0 + 1                    ;Create a "flat" flux vector
   #    IDL> fm_unred, w, f, -0.1, fnew  ;Redden (negative E(B-V)) flux vector
   #    IDL> plot,w,fnew                   
   #
   # NOTES:
   #    (1) The following comparisons between the FM curve and that of Cardelli,
   #        Clayton, & Mathis (1989), (see ccm_unred.pro):
   #        (a) - In the UV, the FM and CCM curves are similar for R < 4.0, but
   #              diverge for larger R
   #        (b) - In the optical region, the FM more closely matches the
   #              monochromatic extinction, especially near the R band.
   #    (2)  Many sightlines with peculiar ultraviolet interstellar extinction 
   #            can be represented with the FM curve, if the proper value of 
   #            R(V) is supplied.
   #    (3) Use the 4 parameter calling sequence if you wish to save the 
   #            original flux vector.
   # PROCEDURE CALLS:
   #       CSPLINE(), POLY()
   # REVISION HISTORY:
   #       Written   W. Landsman        Raytheon  STX   October, 1998
   #       Based on FMRCurve by E. Fitzpatrick (Villanova)
   #       Added /LMC2 and /AVGLMC keywords,  W. Landsman   August 2000
   #       Added ExtCurve keyword, J. Wm. Parker   August 2000
   #       Assume since V5.4 use COMPLEMENT to WHERE  W. Landsman April 2006

   #if not keyword_set(R_V) then R_V = 3.1
   
   x = 10000./ wave.astype(num.float64)  # Convert to inverse microns 
   curve = x*0.
   
   # Set default values of c1,c2,c3,c4,gamma and x0 parameters
   
   if lmc2:
      x0,gamma,c4,c3,c1,c1 = (4.626, 1.05, 0.42, 1.92, 1.31, -2.16)
   elif avglmc:
      x0,gamma,c4,c3,c2,c1 = (4.596, 0.91, 0.64, 2.73, 1.11, -1.28)
   else:
      x0,gamma,c3,c4 = (4.596,0.99,3.23,0.41)
      c2 = -0.824 + 4.717/R_V 
      c1 = 2.030 - 3.007*c2
   
   # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and 
   # R-dependent coefficients
    
   xcutuv = 10000.0/2700.0
   xspluv = 10000.0/num.array([2700.0,2600.0])
   iuv = num.greater_equal(x, xcutuv); N_UV = sum(iuv)
   iopir = -iuv;  Nopir = sum(iopir)
   if N_UV > 0: 
      xuv = num.concatenate([xspluv,x[iuv]])
   else:  
      xuv = xspluv
   
   yuv = c1  + c2*xuv
   xuv2 = num.power(xuv,2)
   yuv = yuv + c3*xuv2/(num.power(xuv2-x0**2,2) + xuv2*gamma**2)
   mxuv = num.where(xuv < 5.9, 5.9, xuv)
   yuv = yuv + c4*(0.5392*num.power(mxuv-5.9,2)+0.05644*num.power(mxuv-5.9,3))
   yuv = yuv + R_V
   yspluv  = yuv[0:2]                  # save spline points
   
   if (N_UV > 0): 
      curve[iuv] = yuv[2:]      # remove spline points
    
   # Compute optical portion of A(lambda)/E(B-V) curve
   # using cubic spline anchored in UV, optical, and IR
   
   xsplopir = num.concatenate([[0],
       10000.0/num.array([26500.0,12200.0,6000.0,5470.0,4670.0,4110.0])])
   ysplir   = num.array([0.0,0.26469,0.82925])*R_V/3.1 
   ysplop   = num.array([poly(R_V, [-4.22809e-01, 1.00270, 2.13572e-04] ),
               poly(R_V, [-5.13540e-02, 1.00216, -7.35778e-05] ),
               poly(R_V, [ 7.00127e-01, 1.00184, -3.32598e-05] ),
               poly(R_V, [ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04 
                        -4.45636e-05] ) ])
     
   ysplopir = num.concatenate([ysplir,ysplop])
   
   if Nopir > 0:
      spl = UnivariateSpline(num.concatenate([xsplopir,xspluv]),
                     num.concatenate([ysplopir,yspluv]), s=0)
      curve[iopir] = spl(x[iopir])
   
   # Now apply extinction correction to input flux vector
   
   return curve
   #if N_params() EQ 3 then flux = flux * 10.^(0.4*curve) else $
   #   funred = flux * 10.^(0.4*curve)       ;Derive unreddened flux
   #
   #ExtCurve = Curve - R_V

def fm07_full(wave, x0, gamma, c1,c2,c3,c4,c5,o1,o2,o3,R,iscale,ipower=1.84):
   '''The fully-parametrized Fitzpatrick & Massa reddening curve. Wave
   is expected to be in Angstroms. Returns E(\lambda-V)/E(B-V), which
   is R_\lambda - Rv'''
   x = 10000./wave.astype(num.float64) # Convert ot inverse microns
   curve = x*0.
   xcutuv = 10000.0/2700.0
   iuv = num.greater_equal(x, xcutuv); N_UV = sum(iuv)
   iopir = -iuv;  Nopir = sum(iopir)

   # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and 
   # R-dependent coefficients
   xspluv = 10000.0/num.array([2700.0,2600.0])
   x1,x2,x3 = (10000./3300, 10000./4000, 10000./5430)
   if N_UV > 0: 
      xuv = num.concatenate([xspluv,x[iuv]])
   else:  
      xuv = xspluv
   
   yuv = c1  + c2*xuv
   xuv2 = num.power(xuv,2)
   yuv = yuv + c3*xuv2/(num.power(xuv2-x0**2,2) + xuv2*gamma**2)
   yuv = yuv + num.greater(xuv,c5)*c4*num.power(xuv-c5,2)
   yspluv  = yuv[0:2]                  # save spline points
   
   if (N_UV > 0): 
      curve[iuv] = yuv[2:]      # remove spline points
    
   # Compute optical portion of A(lambda)/E(B-V) curve
   # using cubic spline anchored in UV, optical, and IR
   
   xsplir = num.array([0,0.25,0.50,0.75,1.0])
   ysplir = iscale*num.power(xsplir, ipower) - R
   xsplop = num.array([x3,x2,x1])
   ysplop = num.array([o3,o2,o1])
   xspl = num.concatenate([xsplir,xsplop,xspluv])
   yspl = num.concatenate([ysplir,ysplop,yspluv])

   if Nopir > 0:
      spl = UnivariateSpline(xspl, yspl, k=3, s=0)
      curve[iopir] = spl(x[iopir])
   return curve


def fm07(wave, R_V=3.1):
   '''Average reddening law from Fitzpatrick & Massa (2007). This includes
   a variation with R(V) and uses the correlation between R(V) and
   k_NIR to make a one-parameter curve (keeping all other curve parameters
   fixed. Caveats abound.  Read the paper!'''
   # Set default values of fixed constants
   x0,gamma,c1,c2,c3,c4,c5 = (4.592,0.922,-0.175, 0.807, 2.991, 0.319, 6.097)
   O1,O2,O3 = (2.055,1.322,0.000)
   ipower = 1.84
   iscale = -0.83 + 0.63*R_V
   ElEv = fm07_full(wave,x0,gamma,c1,c2,c3,c4,c5,O1,O2,O3,R_V,iscale,ipower)

   return ElEv+R_V

def nataf(wave, R_V=3.1, strict_ccm=False):
   '''CCM modified by David Nataf, private communication.'''
   x = 10000./ wave                ; #Convert to inverse microns 
   a = 0.0*x
   b = 0.0*x
  
   # **************************************************************
   good = num.greater(x, 0.3) * num.less(x, 1.1)   # Infrared
   if len(num.nonzero(good)) > 0: 
      alpha = 1.61 - 0.67*(R_V-3.1)
      a = num.where(good, 0.574 * num.power(x, alpha), a)
      b = num.where(good, -0.527 * num.power(x, alpha), b)
  
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
   if num.sometrue(good):
      y = 1.0*x
      good1 = num.greater(y, 5.9)
      F_a = y*0.0    ; F_b = y*0.0
      if num.sometrue(good1):
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
   if num.sometrue(good):
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


def unred_fm(wave, flux, ebv, R_V=3.1, z=0, avglmc=False, lmc2=False):
   '''Deredden by the Fitzpatrick (1999) Law'''
   R_lambda = fm(wave, R_V, avglmc=avglmc, lmc2=lmc2)
   A_lambda = ebv*R_lambda
   unred_flux = flux * num.power(10.0, 0.4*A_lambda)
   return unred_flux

def unred_fm07(wave, flux, ebv, R_V=3.1, z=0):
   '''Deredden by the Fitzpatrick & Massa (2007) Law'''
   R_lambda = fm07(wave, R_V)
   A_lambda = ebv*R_lambda
   unred_flux = flux * num.power(10.0, 0.4*A_lambda)
   return unred_flux


def unred(wave, flux, ebv, R_V = 3.1, z=0, redlaw='ccm', strict_ccm=0):
   '''
   de-redden (or redden, if you set ebv < 0) a spectrum using color
   excess E(B-V) of [ebv] and [R_V].  Optionally, you can redshift the
   spectrum by [z]. You can choose between the Cardelli, Clayton and Mathis
   (redlaw='ccm') or the Fitzpatrick and Massa (1999) law (redlaw='fm'). If
   using ccm, you can either use the default CCM modified by O'Donnel (1994)
   or use the script CCM by specifying strict_ccm=True.
   '''

   if z > 0:
      wave = wave*(1+z)
  
   if redlaw == 'ccm':
      a,b = ccm(wave, strict_ccm)
      # note CCM gives A_lambda/A_V = a + b/Rv
      # therefore, A_lambda = E(B-V)*R_V*(a + b/Rv) = ebv*(a*R_V + b)
      A_lambda = ebv * (a*R_V + b)
      unred_flux = flux * num.power(10.0, 0.4*A_lambda)
   elif redlaw == 'nataf':
      a,b = nataf(wave, strict_ccm)
      # note CCM gives A_lambda/A_V = a + b/Rv
      # therefore, A_lambda = E(B-V)*R_V*(a + b/Rv) = ebv*(a*R_V + b)
      A_lambda = ebv * (a*R_V + b)
      unred_flux = flux * num.power(10.0, 0.4*A_lambda)
   elif redlaw == 'fm' or redlaw == 'f99':
      R_lambda = fm(wave, R_V, avglmc=False, lmc2=False)
      A_lambda = ebv*R_lambda
      unred_flux = flux * num.power(10.0, 0.4*A_lambda)
      # This is a little trick to return an a and b like CCM
      a = R_lambda/R_V
      b = R_lambda*0
   elif redlaw == 'fm07':
      R_lambda = fm07(wave, R_V)
      A_lambda = ebv*R_lambda
      unred_flux = flux * num.power(10.0, 0.4*A_lambda)
      # This is a little trick to return an a and b like CCM
      a = R_lambda/R_V
      b = R_lambda*0
   else:
      raise ValueError, "Unkwown reddening law %s" % redlaw
   return(unred_flux, a , b)


def R_z(wave, z, R_V = 3.1, strict_ccm=0):
   '''Compute host reddening law for effective observed wavelength wave, at redshift
   z, assuming R_V.'''
   wave = wave/(1+z)
   a,b = ccm(wave, strict_ccm)
   R = R_V*a + b
   return(R, a , b)
