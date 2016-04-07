#!/usr/bin/env python
import matplotlib
matplotlib.use('Tkagg')
from matplotlib.pyplot import *
from snpy import get_sn
from snpy.utils import fit1dcurve
from snpy.version import __version__

# Test computing of bolometric light-curve

s = get_sn('SN2006ax_color_model.snpy')
for f in s.data:
   s.data[f].template()

x2,y2,f2,lims2 = s.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='direct', interpolate='spline', SED='H3', verbose=True)
plot(x2,y2/1e45, 's', label='Direct: Hsiao')

x1,y1,f1,lims1 = s.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='SED', interpolate='spline', refband='B', lam1=lims2[0][0],
      lam2=lims2[0][1], SED='H3', verbose=True)
plot(x1,y1/1e45, 'o', label='SED: Hsiao Tempalte')

x3,y3,f3,lims3 = s.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='direct', interpolate='spline', SED=None, verbose=True, 
      outfile='test_bolo.dat')
plot(x2,y2/1e45, 'd', label='Direct: Vega')

legend()
xlabel('Days after B-maximum')
ylabel('Bolometric Luminosity ($10^{45} erg\cdot s^{-1}$)')

savefig('bolometric.eps')

