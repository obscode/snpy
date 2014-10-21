'''A module to compute the reddening law based on specified filters.

The basic problem is that we observe in a particular filter and are
correcting by the color in two other filters.  So, we want to make
the following correction:

   m1_true = m1_obs - R*E(m2 - m3)

   where E(m2 - m3) = m2_obs - m3_obs - (m2_true - m3_true)

   is the color excess through filters 2 and 3.  We want to be able
   to compute R and also take this value of R and convert back to 
   the CCM R_V parameter.

'''

from snpy.filters import fset
from snpy import kcorr
import deredden
import numpy as np
import scipy
import os
import pickle

base = os.path.dirname(__file__)
if os.path.isfile(os.path.join(base,'Ia_R_splines.pickle')):
   f = open(os.path.join(base,'Ia_R_splines.pickle'))
   Bspls = pickle.load(f)
   f.close()

def R_lambda(f, Rv, EBV, redlaw='ccm'):
   '''A fast implementation based on bivariate splines.'''
   scalar = (len(np.shape(Rv)) == 0 and len(np.shape(EBV)) == 0)
   Rv = np.atleast_1d(Rv)
   EBV = np.atleast_1d(EBV)
   #res = np.diag(Bspls[redlaw][f](Rv,EBV))
   res = Bspls[redlaw][f].ev(Rv,EBV)
   if scalar:
      return res[0]
   else:
      return res
   
def Rv_to_R(f1, f2, f3, Rv, EBV=0.01, day=0, redlaw='ccm',
      strict_ccm=0, version='H3'):
   '''Convert from R_V and optionally EBV to an observed R through
   filter f1, corrected by f2-f3 color.  You can choose which day
   the SN SED should be (default day 0, max).  You can also specify
   which reddening law to use:  redlaw='ccm' for Cardelli et al., or
   redlaw='fm' for Fitzpatric and Malla (1999), If redlaw='ccm', you 
   can specify whether we should use the strict CCM relation (default no).
   Finally, you an choose which SED sequence to use:  'H', 'H3', or '91bg'
   (default H3).'''

   # get the SNIa SED
   wave,flux = kcorr.get_SED(day, version=version)
   # Redden according to CCM (plus improvements if strict_ccm=0)
   rflux,a,b = deredden.unred(wave, flux, -EBV, Rv, redlaw=redlaw,
         strict_ccm=strict_ccm)
   
   # compute synthetic magnitudes of original and reddened filter
   m1 = fset[f1].synth_mag(wave, flux)
   m1_red = fset[f1].synth_mag(wave, rflux)
   m2 = fset[f2].synth_mag(wave, flux)
   m2_red = fset[f2].synth_mag(wave, rflux)
   m3 = fset[f3].synth_mag(wave, flux)
   m3_red = fset[f3].synth_mag(wave, rflux)

   # Compute extinctions
   A1 = m1 - m1_red;   A2 = m2 - m2_red;  A3 = m3 - m3_red

   # Return reddening law R
   return(A1/(A2 - A3))


def R_to_Rv(f1, f2, f3, R, EBV=0.01, day=0, strict_ccm=0, version='H3'):
    '''Same as Rv_toR, but work in reverse to find Rv, given R'''

    # A temporary function for which we want to find the zero
    f = lambda x:  Rv_to_R(f1, f2, f3, x, EBV, day, strict_ccm, version)-R

    # Use Brent's method to find the zero of the function:  Rv
    result = scipy.optimize.brentq(f, 0, 10.0)
    return(result)

