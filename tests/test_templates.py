#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from snpy import get_sn
from snpy.utils import fit1dcurve
from snpy.version import __version__

# Test the loading of the objects, interpolating the light-curves, and
# computing k-corrections.

s = get_sn('SN2006ax.txt')
for f in s.data:
   s.data[f].template()
s.kcorr()
s.save('SN2006ax_kcorr.snpy')
s.save('old_saves/SN2006ax_kcorr-%s.snpy' % __version__)

methods = fit1dcurve.functions.keys()

for method in methods:
   s.B.template(method=method)
   s.B.plot()

