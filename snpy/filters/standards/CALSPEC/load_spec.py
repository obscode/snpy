#!/usr/bin/env python

import pyfits
from glob import glob
from numpy import array,savetxt

fout = open('standards.dat', 'w')
files = glob('*.fits')
for fil in files:
   fts = pyfits.open(fil)
   target = fts[0].header['TARGETID']
   if 'SOURCE' in fts[0].header:
      source = fts[0].header['SOURCE']
   else:
      source = ''
   sid = fil.split('.')[0]
   fout.write("%-25s %-30s None CALSPEC %s %s\n" % (sid, fil, target,source))
   savetxt(fil.replace('.fits','.dat'), 
         array([fts[1].data['WAVELENGTH'], fts[1].data['FLUX']]).T)
   fts.close()
fout.close()   
