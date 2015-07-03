#!/usr/bin/env python

# test the fitting methods.
from snpy import get_sn
from snpy.utils import fit1dcurve
print type(str)

funcs = fit1dcurve.functions
print "%-11s %-11s %-11s" % ('Bmax', 'Tmax','dm15')

for func in funcs:
   s = get_sn('SN2006ax_kcorr.snpy')
   s.B.template(method=func)
   print "%6.3f(%s) %6.3f(%s) %6.3f(%s)" % \
         (s.B.Mmax, str(s.B.e_Mmax).split('.')[1][0:3],
          s.B.Tmax, str(s.B.e_Tmax).split('.')[1][0:3],
          s.B.dm15, str(s.B.e_dm15).split('.')[1][0:3])
   
