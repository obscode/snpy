#!/usr/bin/env python

from snpy import get_sn
from snpy.utils import fit1dcurve

# Test the loading of the objects, interpolating the light-curves, and
# computing k-corrections.

s = get_sn('SN2006ax.txt')
for f in s.data:
   s.data[f].template()
s.kcorr()
s.save('SN2006ax_kcorr.snpy')

methods = fit1dcurve.functions.keys()

for method in methods:
   s.B.template(method=method)
   s.B.plot()

