'''
Module for SNooPy to download/parse data from the Open Supernova Catalog.
'''

import json
import urllib
from astropy.coordinates import Angle
from snpy import sn,lc
from numpy import array

# telescope,band --> SNooPy filter database
# We do this by matching (band,system,telescope,observatory) info from the 
# database to SNooPy filters. 
ftrans = {}
ftrans_standard = {}
for band in ['u','g','r','i','B','V','Y','J','H','K']:
   ftrans[(band,"CSP",'',"LCO")] = band
for band in ['U','B','V','R','I']:
   ftrans[(band,'','kait3','')] = band+'kait'
for band in ['J','H','Ks']:
   ftrans[(band,'','PAIRITEL','')] = band+'2m'

# These are for data in (what I'm assuming) would be standard filters.
# We will issue a warning, though.
for band in ['U','B','V','R','I']:
   ftrans_standard[(band,'','','')] = band+"s"
for band in ['u','g','r','i','z']:
   ftrans_standard[(band,'','','')] = band+"_40"



def get_obj(url):

   try:
      u = urllib.urlopen(url)
   except:
      return None,"Invalid URL"
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
   zhel = float(d['redshift'][0]['value'])
   ra = Angle(" ".join([d['ra'][0]['value'],d['ra'][0]['u_value']])).degree
   decl = Angle(" ".join([d['dec'][0]['value'],d['dec'][0]['u_value']])).degree

   snobj = sn(name, ra=ra, dec=decl, z=zhel)

   # Next, the photometry.
   MJD = {}
   mags = {}
   emags = {}
   known_unknowns = []
   unknown_unknowns = []
   for p in d['photometry']:
      t = (p.get('band',''),p.get('system',''),p.get('telescope',''),
           p.get('observatory',''))
      if t in ftrans:
         b = ftrans[t]
      elif t in ftrans_standard:
         b = ftrans_standard[t]
         if t not in known_unknowns:
            known_unknowns.append(t)
            print "Warning:  no telescope/system info, assuming standard ",b
      else:
         # No idea
         unknown_unknowns.append(t)
         continue
      if b not in MJD:  
         MJD[b] = []
         mags[b] = []
         emags[b] = []
      MJD[b].append(float(p['time']))
      mags[b].append(float(p['magnitude']))
      emags[b].append(float(p['e_magnitude']))
   for b in MJD:
      snobj.data[b] = lc(snobj, b, array(MJD[b]), array(mags[b]), 
            array(emags[b]))
      snobj.data[b].time_sort()
   snobj.get_restbands()

   if len(unknown_unknowns) > 0:
      unknown_unknowns = list(set(unknown_unknowns))
      print "Warning:  the following photometry was not recognized by SNooPy"
      print "and was not imported:"
      for item in unknown_unknowns:
         print item
   return(snobj,'Success')
