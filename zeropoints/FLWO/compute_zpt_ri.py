#!/usr/bin/env python

from snpy import *
from numpy import *
import sys,os,string
from snpy.filters import standards
from snpy.filters import standard_mags
import scipy

# Filters to process
fs = ['rk','ik']

# dictionary of natural magnitudes
sm = {}

sloan_mags = standard_mags['Smith']
BD17 = standards['Smith']['bd17B']
# Convert Smith et al. magnitudes to CfA 3/4 Natural magnitudes. First, we
# need natural and standard V mags of bd17
V_bd17 = 9.43   # From Fukugita et al. (1996)
v_bd17 = fset['Vk'].synth_mag(BD17)
# Color term (period 1) from Hicken et al. (2011)
r_nat = v_bd17 - 1.0684*(V_bd17 - sloan_mags['bd17']['r_40'])
i_nat = v_bd17 - 1.0239*(V_bd17 - sloan_mags['bd17']['i_40'])

print "Assuming natural magnitdes for BD+17:"
print "   V_nat = ", v_bd17
print "   r_nat = ", r_nat
print "   i_nat = ", i_nat

print "Zero-points:"
print "zp(r_nat) = ", fset['rk'].compute_zpt(BD17, r_nat)
print "zp(i_nat) = ", fset['ik'].compute_zpt(BD17, i_nat)
