'''
Module for SNooPy to download/parse data from the Open Supernova Catalog.
'''

import json
import urllib
from astropy.coordinates import Angle
from snpy import sn,lc
from numpy import array,log10

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
for band in ['B','V','R','I']:
   ftrans[(band,'','kait4', '')] = band+'kait'
for band in ['U','V','B']:
   ftrans[(band, 'Vega','Swift','')] = band+"_UVOT"
for band in ['UVW1','UVW2','UVM2']:
   ftrans[(band, 'Vega','Swift','')] = band
for band in ['g','r','i','z']:
   ftrans[(band, '', 'PS1','')] = "ps1_"+band

# These are for data in (what I'm assuming) would be standard filters.
# We will issue a warning, though.
for band in ['U','B','V','R','I']:
   ftrans_standard[(band,'','','')] = band+"s"
for band in ['u','g','r','i','z']:
   ftrans_standard[(band,'','','')] = band+"_40"
for band in ["u'","g'","r'","i'","z'"]:
   ftrans_standard[(band,'','','')] = band[0]+"_40"


warning_message = {
      'upperlims':'Warning: Data with upper-limits not imported',
      }

def get_obj(url, full_data=False):

   try:
      u = urllib.urlopen(url)
   except:
      if full_data:
         return None, "Invalid URL", None
      return None,"Invalid URL"
   try:
      d = json.load(u)
   except:
      u.close()
      if full_data:
         return None,"Failed to decode JSON",None
      return None,"Failed to decode JSON"
   else:
      u.close()
   
   # We now have the JSON data. Get the info we need
   d = d.values()[0]
   name = d['name']
   if 'redshift' not in d or 'ra' not in d or 'dec' not in d:
      if full_data:
         return None,"No redshift, RA, or DEC found",d
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
   warnings = []
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
      if 'time' in p and 'magnitude' in p and 'e_magnitude' in p:
         MJD[b].append(float(p['time']))
         mags[b].append(float(p['magnitude']))
         emags[b].append(float(p['e_magnitude']))
      elif 'time' in p and 'countrate' in p and 'e_countrate' in p \
            and 'zeropoint' in p:
         if float(p['countrate']) < 0: continue
         MJD[b].append(float(p['time']))
         mags[b].append(-2.5*log10(float(p['countrate'])) + \
               float(p['zeropoint']))
         emags[b].append(1.087*float(p['e_countrate'])/float(p['countrate']))
      else:
         if 'upperlims' not in warnings:
            warnings.append('upperlims')
   for b in MJD:
      if len(MJD[b]) > 0:
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
   if warnings:
      for warning in warnings:
         print warning_message[warning]
   if full_data:
      return(snobj, 'Success', d)
   return(snobj,'Success')
