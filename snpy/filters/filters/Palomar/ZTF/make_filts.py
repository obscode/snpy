#!/usr/bin/env python

from scipy.interpolate import splrep,splev
from scipy.signal import convolve, boxcar
from numpy import loadtxt, savetxt, array

# original filter functions from the SVO Filter Profile Service V2
wg,rg = loadtxt('ztf_g_band.dat', unpack=True)
wr,rr = loadtxt('ztf_g_band.dat', unpack=True)
wi,ri = loadtxt('ztf_g_band.dat', unpack=True)

# Atmosphere extinction
wa,atm = loadtxt('atm_ext_ctio.dat', unpack=True)
wl,lines = loadtxt('atm_ext_lines.dat', unpack=True)

# First, down-sample the ztf filters with boxcar average
window = boxcar(10)
window = window/sum(window)   # normalize

wg = wg[::10]
rg = convolve(rg, window, mode='same')[::10]
wr = wr[::10]
rr = convolve(rr, window, mode='same')[::10]
wi = wi[::10]
ri = convolve(ri, window, mode='same')[::10]

# We need to interpolate over atmospheric functions
tck_atm = splrep(wa, atm, k=3)   # Smooth function
tck_lines = splrep(wa, atm, k=1, s=0)  # less smooth

# g-band only needs atm, no lines
rg = rg * splev(wg, tck_atm)
savetxt('ztf_g_atm.dat', array([wg,rg]).T, fmt='%.8f')

# r and i bands need the lines
rr = rr * splev(wr, tck_atm) * splev(wr, tck_lines)
savetxt('ztf_r_atm.dat', array([wr,rr]).T, fmt='%.8f')
ri = ri * splev(wi, tck_atm) * splev(wi, tck_lines)
savetxt('ztf_i_atm.dat', array([wi,ri]).T, fmt='%.8f')
