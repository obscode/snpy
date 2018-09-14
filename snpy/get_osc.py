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
standard_warnings = {}
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
   standard_warnings[band] = "Johnson/Kron/Cousins "
for band in ['u','g','r','i','z']:
   ftrans_standard[(band,'','','')] = band+"_40"
   standard_warnings[band] = "Sloan (APO) "
for band in ["u'","g'","r'","i'","z'"]:
   ftrans_standard[(band,'','','')] = band[0]+"_40"
   standard_warnings[band] = "Sloan (USNO-40) "
for band in ["J","H","Ks"]:
   ftrans_standard[(band[0],'','','')] = band+"2m"
   standard_warnings[band[0]] = "2MASS "


warning_message = {
      'upperlims_noerr':'Warning: Data lacking errorbars or with upper-limits not imported',
      'upperlims':'Warning: Data with upper-limits not imported',
      }

def get_obj(url, full_data=False, allow_no_errors=False, missing_error=0.01):
   '''Attempt to build a SNooPy object from a Open Supernova Catalog server
   URL.'''

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
   # All primary sources
   all_sources_dict = [item for item in d['sources'] \
                           if not item.get('secondary',False)]
   all_sources = {}
   for source in all_sources_dict:
      all_sources[source['alias']] = (source.get('bibcode',''),
                                       source.get('reference',''))

   # Next, the photometry.
   used_sources = []
   MJD = {}
   mags = {}
   emags = {}
   sids = {}
   known_unknowns = []
   unknown_unknowns = []
   warnings = []
   for p in d['photometry']:
      t = (p.get('band',''),p.get('system',''),p.get('telescope',''),
           p.get('observatory',''))

      # Deal with source of photometry
      ss = p.get('source').split(',')
      this_source = None
      for s in ss:
         if s in all_sources:
            this_source = all_sources[s]
            break
      if this_source is None:
         print "Warning:  no primary source, skipping"
         continue

      if t in ftrans:
         b = ftrans[t]
      elif t in ftrans_standard:
         b = ftrans_standard[t]
         if t not in known_unknowns:
            known_unknowns.append(t)
            print "Warning:  no telescope/system info, assuming ", \
                  standard_warnings[b[0]], b[0]
      elif (t[0],"","","") in ftrans_standard:
         b = ftrans_standard[(t[0],"","","")]
         if t not in known_unknowns:
            known_unknowns.append(t)
            print "Warning: telescope/system defined by %s/%s/%s not "\
                  "recognized, assuming %s %s" %\
                  (t[1],t[2],t[3],standard_warnings[t[0]],t[0])
      else:
         # No idea
         if t not in unknown_unknowns:
            unknown_unknowns.append(t)
            print "Warning: telescope/system defined by %s/%s/%s not "\
                  "recognized and can't figure out the filter %s" % \
                  (t[1],t[2],t[3],t[0])
         unknown_unknowns.append(t)
         continue
      if b not in MJD:  
         MJD[b] = []
         mags[b] = []
         emags[b] = []
         sids[b] = []
      if 'time' in p and 'magnitude' in p:
         if not allow_no_errors and 'e_magnitude' not in p:
            if 'upperlims' not in warnings: warnings.append('upperlims')
            continue
         MJD[b].append(float(p['time']))
         mags[b].append(float(p['magnitude']))
         emags[b].append(float(p.get('e_magnitude',missing_error)))
      elif 'time' in p and 'countrate' in p and 'zeropoint' in p:
         if not allow_no_errors and 'e_countrate' not in p:
            if 'upperlims' not in warnings: warnings.append('upperlims')
            continue
         if float(p['countrate']) < 0: continue
         MJD[b].append(float(p['time']))
         mags[b].append(-2.5*log10(float(p['countrate'])) + \
               float(p['zeropoint']))
         ec = p.get('e_countrate',None)
         if ec is not None:
            emags[b].append(1.087*float(p['e_countrate'])/float(p['countrate']))
         else:
            emags[b].append(missing_error)
      else:
         if 'upperlims_noerr' not in warnings:
            warnings.append('upperlims_noerr')
         continue
      if this_source not in used_sources:
         used_sources.append(this_source)
      # At this point we're actually using the photometry, so find source
      sid = used_sources.index(this_source)
      sids[b].append(sid)
   for b in MJD:
      if len(MJD[b]) > 0:
         snobj.data[b] = lc(snobj, b, array(MJD[b]), array(mags[b]), 
               array(emags[b]), sids=array(sids[b], dtype=int))
         snobj.data[b].time_sort()
   snobj.sources = used_sources
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
