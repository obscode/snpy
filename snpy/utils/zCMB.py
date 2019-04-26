'''
Module for computing the CMB redshift from the heliocentric redshift and
coordinates. For now, uses NED calculatro, but maybe later will do all the 
math, so we don't need the interweb.

Yes, let's do the math. It's not hard:  using astropy.coordinates, you can
get the angular separation easily. We'll leave the NED version here for
posterity, though.'''

import six
if six.PY3:
   import urllib.request as urllib
else:
   import urllib
import re
from astropy.coordinates import SkyCoord
from math import sin,cos,pi,sqrt
url_temp = "http://ned.ipac.caltech.edu/cgi-bin/velc?lon=%.6fd&in_csys=Equatorial&lat=%.6fd&in_equinox=J2000.0&vel=%.3f&vfrom=Heliocentric&vto=3K&alon=&a_csys=Equatorial&alat=&a_equinox=J2000.0&avel=0.0"
pat = re.compile(r'([0-9]+(\.[0-9]*)?)')

# The following from Fixsen et al. (1996)
CMB_l = 264.14
CMB_b = 48.26
CMB_V = 371.0

def z_cmb(z_hel, ra, dec):
   '''heliocentric to CMB redshift conversion.
      z_hel:  heliocentric redshift of object
      ra,dec: coorinates in degrees.
   '''
   v_hel = z_hel*3e5
   cmb_dir = SkyCoord(l=CMB_l, b=CMB_b, frame='galactic', unit='deg')
   obj_dir = SkyCoord(ra, dec, frame='icrs', unit='deg')
   ang_sep = cmb_dir.separation(obj_dir).value*pi/180
   v_corr = CMB_V*cos(ang_sep)
   vcmb = v_hel + v_corr
   return vcmb/3e5

def v_LG(z_hel, ra, dec):
   '''heliocentric to LG velocity conversion.
      z_hel:  heliocentric redshift of object
      ra,dec: coorinates in degrees.
   '''
   v_hel = z_hel*3e5
   obj_dir = SkyCoord(ra, dec, frame='icrs', unit='deg')
   l = obj_dir.galactic.l.to('rad').value
   b = obj_dir.galactic.b.to('rad').value
   V_LG = v_hel - 79*cos(l)*cos(b) + 296*sin(l)*cos(b)-36*sin(b)
   return V_LG

def v_infall(Vo, ra, dec, Vfid, Va, ra_a, dec_a, gamma=2):
   '''Compute an infall velocity due to an attractor.  See Mould+2000
      ra,dec:  coordinates of object in deg
      Vo:      LG-frame velocity of object
      Vfid:    amplitude of attractor at LG
      Va:      distance to attractor in km/s
      la,lb:   galactic coordinates of attractor
      gamma:   slope of attractor's density profile
   '''
   obj_dir = SkyCoord(ra, dec, frame='icrs', unit='deg')
   att_dir = SkyCoord(ra_a, dec_a, frame='icrs', unit='deg')
   theta = att_dir.separation(obj_dir).to('rad').value
   # projected distance
   roa = sqrt(Vo**2 + Va**2 - 2*Va*Vo*cos(theta))
   Vinf = Vfid*cos(theta) + Vfid*((Vo - Va*cos(theta))/roa)* pow(roa/Va,1-gamma)
   return Vinf

def z_cmb_NED(z_hel, ra, dec):
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
      raise ValueError("Could not parse output from NED")
   return vcmb/3e5
