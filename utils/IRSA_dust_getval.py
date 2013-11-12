#!/usr/bin/env python
'''A module to make WEB querries to the IRSA Schlegel dust map calulator
and extract E(B-V)'''

import urllib
import re
from xml.dom.minidom import parse

debug = 0

BASE_URL = "http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=%.5f+%.5f"

def get_dust_RADEC(ra, dec, calibration="SF11"):
   '''Query IRSA with given ra and dec to get the E(B-V).  To remain compatible
   with dust_getval module, return a list (dust_getval results an array).
   You can specify calibration of "SF11" or "SFD98"'''
   if debug:
      print "get_dust_RADEC:  Querying URL:  ",BASE_URL % (ra,dec)
   try:
      u = urllib.urlopen(BASE_URL % (ra,dec))
   except:
      print "Failed to connect to IRSA.  E(B-V) query failed"
      return([None],[1])
   if not u:
      print "Failed to connect to IRSA.  E(B-V) query failed"
      return([None],[1])
   dom = parse(u)
   u.close()
   result = None
   tag = {'SF11':'meanValueSandF',
          'SFD98':'meanValueSFD'}
   try:
      EBVstr = dom.getElementsByTagName(tag[calibration])[0].childNodes[0].data
      result = float(EBVstr.strip().split()[0])
   except:
      print "E(B-V) query to IRSA failed.  Check RA/DEC and make" +\
            " sure you have internet access."
      return ([None],[1])
   return [result],[1]

def get_dust_sigma_RADEC(ra, dec, calibration="SF11"):
   '''Query IRSA with given ra and dec to get the error in E(B-V). 
   You can specify calibration of "SF11" or "SFD98"'''
   if debug:
      print "get_dust_RADEC:  Querying URL:  ",BASE_URL % (ra,dec)
   try:
      u = urllib.urlopen(BASE_URL % (ra,dec))
   except:
      print "Failed to connect to IRSA.  E(B-V) query failed"
      return([None],[1])
   if not u:
      print "Failed to connect to IRSA.  E(B-V) query failed"
      return([None],[1])
   dom = parse(u)
   u.close()
   result = None
   tag = {'SF11':'stdSandF',
          'SFD98':'stdSFD'}
   try:
      EBVstr = dom.getElementsByTagName(tag[calibration])[0].childNodes[0].data
      result = float(EBVstr.strip().split()[0])
   except:
      print "E(B-V) query to IRSA failed.  Check RA/DEC and make" +\
            " sure you have internet access."
      return None
   return result
