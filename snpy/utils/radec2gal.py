'''
Module for computing the galactic coordinates from stored RA/DEC
coordinates. If we have astropy,use that, otherwise try NED calculator'''

import urllib, re
try:
   from astropy import coordinates,units
except:
   coordinates,units = None,None

url_temp = "http://ned.ipac.caltech.edu/cgi-bin/calc?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&lon=%.5fd&lat=%.5fd&pa=0.0&out_csys=Galactic&out_equinox=J2000.0"
pat = re.compile(r'([0-9]+(\.[0-9]*)?)')

def radec2gal(ra, dec):
   if coordinates is not None:
      c = coordinates.SkyCoord(ra=ra*units.degree, dec=dec*units.degree,
            frame='icrs')
      l = c.galactic.l.deg
      b = c.galactic.b.deg
      return l,b
   u = urllib.urlopen(url_temp % (ra, dec))
   lines = u.readlines()
   indx = None
   for i,line in enumerate(lines):
      if line.find('Output:') >= 0:
         indx = i+1
         break
   if index is None:
      raise ValueError, "Could not parse output from NED"
   vals = pat.findall(lines[indx])
   return float(vals[0]),float(vals[1])
