'''
Module for SNooPy to download/parse data from the Open Supernova Catalog.
'''

import json
import urllib
from astropy.coordinates import Angle
from snpy import sn

# telescope,band --> SNooPy filter database
# We do this by matching (band,system,telescope,observatory) info from the 
# database to SNooPy filters. Failing that, we make a new filter name from 
# the info
ftrans = {}
for band in ['u','g','r','i','B','V','Y','J','H','K']:
   ftrans[(band,"CSP",'',"LCO")] = band
for band in ['U','B','V','R','I']:
   ftrans[(band,'','kait3','')] = band+'kait'
for band in ['J','H','Ks']:
   ftrans[(band,'','PAIRITEL','')] = band+'2m'



def get_obj(url):

   try:
      u = urllib.urlopen(url)
   except:
      retrun None,"Invalid URL"
   try:
      d = json.load(u)
   except:
      u.close()
      return None,"Failed to decode JSON"
   else:
      u.close()
   
   # We now have the JSON data. Get the info we need
   d = d.values()[0]
   name = d['name']
   if 'redshift' not in d or 'ra' not in d or 'dec' not in d:
      return None,"No redshift, RA, or DEC found"
   zhel = float(d['redshift']['value'])
   ra = Angle(" ".joind([d['ra'][0]['value'],d['ra'][0]['u_value']])).degree
   decl = Angle(" ".joind([d['dec'][0]['value'],d['dec'][0]['u_value']])).degree

   snobj = sn(name, ra=ra, dec=decl, z=zhel)



