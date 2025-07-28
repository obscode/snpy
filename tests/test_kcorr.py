#!/usr/bin/env python

# This script creates a fake SN spectrum, puts it at several redshifts
# and compares with kcorr

import pytest
from snpy import kcorr,getSED,fset,mangle_spectrum
from numpy import *

def test_kcorr():
   filters = ['u','B','V','r','i','Y','J','H']

   wave0,flux0 = getSED(0, version="H3")
   # redden the flux a bit
   rflux = kcorr.redden(wave0,flux0, 0.1, 0.0, 0.0)
   # generate artificial photometry at redshift 0.1:
   mag0s = []
   magzs = []
   for f in filters:
      mag0s.append(fset[f].synth_mag(wave0, rflux, z=0.))
      magzs.append(fset[f].synth_mag(wave0, rflux, z=0.1)+2.5*log10(1.1))
   
   mag0s = array([mag0s])
   magzs = array([magzs])
   masks = greater(magzs, 0)
   k0s = array([kcorr.kcorr([0], f, f, 0.1)[0][0] for f in filters])
   ks,mask,Rts,m_opts = kcorr.kcorr_mangle([0], filters, magzs, masks, filters, 
         0.1, version='H3', full_output=True, gradient=True, normfilter='H')
   # difference in kcorrs with no mangling:
   delta1 = absolute((magzs[0,:]-mag0s)[0] - k0s)
   # difference in kcorrs with mangling:
   delta2 = absolute((magzs[0,:]-mag0s)[0] - ks[0])

   # The mangling should be better and precise
   assert all(delta2 < delta1) and all(delta2 < 0.01)
   
def test_SEDs():
   wave1,flux1 = kcorr.get_SED(0, version='H3')
   wave2,flux2 = kcorr.get_SED(0, version='H3+L')

   delta1 = fset['B'].synth_mag(wave1, flux1) - fset['B'].synth_mag(wave2,flux2)
   delta2 = fset['H'].synth_mag(wave1, flux1) - fset['H'].synth_mag(wave2,flux2)

   assert (delta1 < 0.001) and (delta2 - 0.785952 < 0.001)


