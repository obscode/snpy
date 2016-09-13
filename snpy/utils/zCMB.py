'''
Module for computing the CMB redshift from the heliocentric redshift and
coordinates. For now, uses NED calculatro, but maybe later will do all the 
math, so we don't need the interweb'''

import urllib, re
url_temp = "http://ned.ipac.caltech.edu/cgi-bin/velc?lon=%.6fd&in_csys=Equatorial&lat=%.6fd&in_equinox=J2000.0&vel=%.3f&vfrom=Heliocentric&vto=3K&alon=&a_csys=Equatorial&alat=&a_equinox=J2000.0&avel=0.0"
pat = re.compile(r'([0-9]+(\.[0-9]*)?)')

def z_cmb(z_hel, ra, dec):
   v_hel = z_hel*3e5
   u = urllib.urlopen(url_temp % (ra, dec, v_hel))
   lines = u.readlines()
   vcmb = None
   for line in lines:
      if line.find('Output:') >= 0:
         res = pat.search(line)
         if res is not None:
            vcmb = float(res.groups()[0])
   if vcmb is None:
      raise ValueError, "Could not parse output from NED"
   return vcmb/3e5
