#!/usr/bin/env python

# This script creates a fake SN spectrum, puts it at several redshifts
# and compares with kcorr

import matplotlib
matplotlib.use('Agg')
from snpy import kcorr,getSED,fset,mangle_spectrum
from numpy import *
from matplotlib import pyplot as plt

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
      0.1, version='H3', full_output=True, normfilter='H')
print "Expected kcorrs:\n", (magzs[0,:]-mag0s)[0]
print "Computed kcorrs with no mangling\n", k0s 
print "Computed kcorrs with tspline mangling:\n", ks[0]
ks,mask,Rts,m_opts2 = kcorr.kcorr_mangle([0], filters, magzs, masks, filters, 
      0.1, version='H3', method='spline', full_output=True, normfilter='H')
print "Computed kcorrs with spline mangling:\n", ks[0]

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(wave0, flux0)
ax.plot(wave0, rflux/flux0, '-', label="true mangling")
mflux = mangle_spectrum.apply_mangle(wave0,flux0, **m_opts[0])[0]
# Normalize, since only colors count in mangling
norm = fset['H'].response(wave0,rflux)/fset['H'].response(wave0,mflux)
ax.plot(wave0, mflux*norm/flux0, '-', label="tspline mangling")
mflux = mangle_spectrum.apply_mangle(wave0,flux0, **m_opts2[0])[0]
norm = fset['H'].response(wave0,rflux)/fset['H'].response(wave0,mflux)
ax.plot(wave0, mflux*norm/flux0, '-', label="spline mangling")
ax.legend(loc='upper left', fontsize=12)
for filt in filters:
   ax.plot(fset[filt].wave,fset[filt].resp,'-', color='blue', alpha=0.6)
   ax.fill_between(fset[filt].wave, fset[filt].resp, color='blue', alpha=0.4)
plt.savefig('kcorr.pdf')
