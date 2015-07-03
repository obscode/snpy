#!/usr/bin/env python

from fit_lightcurves import *
from pygplot import *
import numpy as num

f = open('stretch.dat','w')

for name in sys.argv[1:]:
   data = columns(name)
   data2 = columns(name.replace("template","sigma"))

   bands = ['B','V','R','I']
   s = super('test')
   for i in range(4):
      s.data[bands[i]] = {}
      s.data[bands[i]]['MJD'] = data[0]
      s.data[bands[i]]['mag'] = data[i+1]
      s.data[bands[i]]['e_mag'] = data2[i+1]

   s.mask['B'] = num.greater(s.data['B']['MJD'], 0)
   s.z = 0
   s.get_restbands()
   s.s = 1
   gids = num
   s.fit(['B'], dm15=1.1, EBVhost=0, Tmax=0)
   
   print >>f, name, s.s

f.close()
