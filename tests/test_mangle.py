#!/usr/bin/env python

# This script creates a fake SN spectrum, warps it in several ways and tests
# the methods used to warp it back.

import pytest
from snpy import kcorr,getSED,fset,mangle_spectrum
from numpy import *

def test_mangle_mpfit_nogradient():
   # Test mangling with the mpfit non-linear least-squares fitter
   # using no gradient at endpoints
   filters = ['u','B','V','r','i','Y','J','H']

   wave0,flux0 = getSED(0, version="H3")
   # redden the flux a bit
   rflux = kcorr.redden(wave0,flux0, 0.1, 0.0, 0.0)
   mags = array([fset[f].synth_mag(wave0, rflux) for f in filters])
   gids = greater(wave0, fset['u'].wave[0])*less(wave0, fset['H'].wave[-1])

   for method in mangle_spectrum.mangle_functions:
      mflux,awaves,pars = mangle_spectrum.mangle_spectrum2(wave0,flux0, filters,
            mags, method=method, lstsq=False)
      p_diff = sum((mflux[0][gids]-rflux[gids])/rflux[gids])/sum(gids)*100
      assert(p_diff < 1.)

def test_mangle_mpfit_gradient():
   # Test mangling with the mpfit non-linear least-squares fitter
   # using gradient at endpoints
   filters = ['u','B','V','r','i','Y','J','H']

   wave0,flux0 = getSED(0, version="H3")
   # redden the flux a bit
   rflux = kcorr.redden(wave0,flux0, 0.1, 0.0, 0.0)
   mags = array([fset[f].synth_mag(wave0, rflux) for f in filters])
   gids = greater(wave0, fset['u'].wave[0])*less(wave0, fset['H'].wave[-1])

   for method in mangle_spectrum.mangle_functions:
      mflux,awaves,pars = mangle_spectrum.mangle_spectrum2(wave0,flux0, filters,
            mags, method=method, gradient=True, lstsq=False)
      p_diff = sum((mflux[0][gids]-rflux[gids])/rflux[gids])/sum(gids)*100
      assert(p_diff < 0.1)

def test_mangle_lstsq():
   filters = ['u','B','V','r','i','Y','J','H']

   wave0,flux0 = getSED(0, version="H3")
   # redden the flux a bit
   rflux = kcorr.redden(wave0,flux0, 0.1, 0.0, 0.0)
   mags = array([fset[f].synth_mag(wave0, rflux) for f in filters])
   gids = greater(wave0, fset['u'].wave[0])*less(wave0, fset['H'].wave[-1])

   for method in mangle_spectrum.mangle_functions:
      mflux,awaves,pars = mangle_spectrum.mangle_spectrum2(wave0,flux0, filters,
            mags, method=method, gradient=True)
      p_diff = sum((mflux[0][gids]-rflux[gids])/rflux[gids])/sum(gids)*100
      assert(p_diff < 0.1)

