#!/usr/bin/env python
'''A module to make WEB querries to the NED Schlegel dust map calulator
and extract E(B-V)'''

import urllib
import re

debug = 0

BASE_URL = "http://nedwww.ipac.caltech.edu/cgi-bin/nph-calc?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&lon=%sd&lat=%sd&pa=0.0&out_csys=Equatorial&out_equinox=J2000.0"
pat = re.compile('E\(B-V\)\s*=\s*([0-9\.]+)\s*mag')


def get_dust_RADEC(ra, dec):
   '''Query NED with given ra and dec to get the E(B-V).  To remain compatible
   with dust_getval module, return a list (dust_getval results an array).'''
   if debug:
      print "get_dust_RADEC:  Querying URL:  ",BASE_URL % ("%.5f"%(ra),"%.5f"%(dec))
   u = urllib.urlopen(BASE_URL % ("%.5f"%(ra),"%.5f"%(dec)))
   if not u:
      print "Failed to connect to NED.  E(B-V) query failed"
      return([None],[1])
   lines = u.readlines()
   if debug:
      print "get_dust_RADEC:  Got back data:", lines
   u.close()
   result = None
   for line in lines:
      if pat.findall(line):
         result = float(pat.findall(line)[0])
   if result is None:
      print "E(B-V) query to NED failed.  Check RA/DEC and make" +\
            " sure you have internet access."
   return [result],[1]
