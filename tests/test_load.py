#!/usr/bin/env python
import glob
from snpy import get_sn

# Test to see if save files from previous versions of SNooPy can load.
for file in glob.glob('old_saves/*_kcorr_*.snpy'):
   try:
      s = get_sn(file)
   except:
      print "Failed to load old SNooPy file", file
