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
plot(x2,y2/1e43, 's', label='Direct: Hsiao')

x4,y4,f4,lims4 = s.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='direct', interpolate='model', extrapolate=True, SED='H3', 
      verbose=True)
plot(x4,y4/1e43, 's', label='Direct: Hsiao+extrap')

# Try to use consistent integration limits when comparing
l0 = min([lim[0] for lim in lims4])
l1 = max([lim[1] for lim in lims4])

x1,y1,f1,lims1 = s.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='SED', interpolate='spline', refband='B', lam1=l0,
      lam2=l1, SED='H3', verbose=True)
plot(x1,y1/1e43, 'o', label='SED: Hsiao Tempalte')

x3,y3,f3,lims3 = s.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='direct', interpolate='spline', SED=None, verbose=True, 
      outfile='test_bolo.dat')
plot(x3,y3/1e43, 'd', label='Direct: Vega')

legend()
xlabel('Days after B-maximum')
ylabel('Bolometric Luminosity ($10^{45} erg\cdot s^{-1}$)')

savefig('bolometric.eps')

