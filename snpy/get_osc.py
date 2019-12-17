'''
Module for SNooPy to download/parse data from the Open Supernova Catalog.
'''

from __future__ import print_function
import six
import json
if six.PY3:
   import urllib.request as urllib
else:
   import urllib
from astropy.coordinates import Angle
from snpy import sn,lc,fset
from numpy import array,log10
import astropy.units as u
from snpy.filters import spectrum
from snpy.specobj import timespec

def CfAbands(filt, MJD):
   if MJD < 51913.0:
      return filt[0]+'s'  # standard photometry
   elif 51913.0 < MJD < 55058:
      if filt[0] == 'U': return 'U4sh'
      if filt[0] == 'I': return 'I4sh'
      if filt[0] == 'R': return 'R4sh'
      return filt[0]+'k1' # natural photometry CfA3 + CfA4 period 1
   else:
      if filt[0] == 'U': return 'U4sh'
      if filt[0] == 'I': return 'I4sh'
      if filt[0] == 'R': return 'R4sh'

      return filt[0]+'k2' # natural photometry CfA4 period 2

# Some well-known publications and their mappings:
pubs = {
   '1999AJ....117..707R':  # Riess et al. (1999) Standard Photometry
      CfAbands,
   '2006AJ....131..527J':  # Jha et al. (2006) Standard Photometry
      CfAbands,
   '2009ApJ...700..331H':  # Hicken et al. (2009) CfA3 Natural Photometry
      CfAbands,
   '2012ApJS..200...12H':  # Hicken et al. (2012) CfA4 Natural Photometry
      CfAbands
      }

# telescope,band --> SNooPy filter database
# We do this by matching (band,system,telescope,observatory) info from the 
# database to SNooPy filters. 
ftrans = {}
ftrans_standard = {}
standard_warnings = {}
for band in ['u','g','r','i','B','V','Y','J','H','K']:
   ftrans[(band,"CSP",'',"LCO")] = band
for band in ['U','B','V','R','I']:
   ftrans[(band,'','kait2','')] = band+'kait'
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

# Our own photometric systems:
def CSP_systems(filt, MJD):
   '''Given a filter name and MJD date, output the correct telescope and
   system information.'''
   if filt == "V":
      if MJD < 53748.0:
         return (dict(telescope='Swope',instrument='Site2',band='V-3014',
            zeropoint="{:.4f}".format(fset['V0'].zp)))
      elif MJD < 53759.0:
         return (dict(telescope='Swope',instrument='Site2',band='V-3009',
            zeropoint="{:.4f}".format(fset['V1'].zp)))
      elif MJD < 56566.0:
         return (dict(telescope='Swope',instrument='Site2',band='V-9844',
            zeropoint="{:.4f}".format(fset['V'].zp)))
      else:
         return (dict(telescope='Swope',instrument='e2v',band='V-9844',
            zeropoint="{:.4f}".format(fset['V2'].zp)))
   if filt == "Jrc2":
         return (dict(telescope='Swope',instrument='RetroCam',band='J',
            zeropoint="{:.4f}".format(fset[filt].zp)))
   if filt in ['u','g','r','i','B']:
      if MJD < 56566.0:
         return (dict(telescope='Swope',instrument='Site2',band=filt,
            zeropoint="{:.4f}".format(fset[filt].zp)))
      else:
         return (dict(telescope='Swope',instrument='e2v',band=filt,
            zeropoint="{:.4f}".format(fset[filt+'2'].zp)))
   if filt in ['Y','J','H']:
      if MJD < 55743.0:
         return (dict(telescope='Swope',instrument='RetroCam',band=filt,
            zeropoint="{:.4f}".format(fset[filt].zp)))
      else:
         return (dict(telescope='DuPont',instrument='RetroCam',band=filt,
            zeropoint="{:.4f}".format(fset[filt+'d'].zp)))
   return({})

MJD_offsets = {
      'MJD':0,
      'JD':-2400000.5
      }


warning_message = {
      'upperlims_noerr':'Warning: Data lacking errorbars or with upper-limits not imported',
      'upperlims':'Warning: Data with upper-limits not imported',
      }

OSC_template = '''https://sne.space/astrocats/astrocats/supernovae/output/json/{}.json'''

def get_obj(url, full_data=True, allow_no_errors=False, missing_error=0.01):
   '''Attempt to build a SNooPy object from a Open Supernova Catalog server
   URL.'''

   if url.find('osc:') == 0:
      # Try to construct a url based only on a name.
      url = OSC_template.format(url.split(':')[1])

   try:
      uf = urllib.urlopen(url)
   except:
      return None,"Invalid URL"
   try:
      d = json.load(uf)
   except:
      uf.close()
      if full_data:
         return None,"Failed to decode JSON",None
      return None,"Failed to decode JSON"
   else:
      uf.close()
   
   # We now have the JSON data. Get the info we need
   d = list(d.values())[0]
   name = d['name']
   if 'redshift' not in d or 'ra' not in d or 'dec' not in d:
      return None,"No redshift, RA, or DEC found"
   zhel = float(d['redshift'][0]['value'])
   ra = Angle(" ".join([d['ra'][0]['value'],d['ra'][0]['u_value']])).degree
   decl = Angle(" ".join([d['dec'][0]['value'],d['dec'][0]['u_value']])).degree

   snobj = sn(name, ra=ra, dec=decl, z=zhel)
   # All primary sources
   all_sources_dict = [item for item in d['sources'] \
                           if not item.get('secondary',False)]
   all_sources_dict2 = [item for item in d['sources'] \
                           if item.get('secondary',False)]
   all_sources = {}
   for source in all_sources_dict:
      all_sources[source['alias']] = (source.get('bibcode',''),
                                       source.get('reference',''))
   all_sources2 = {}
   for source in all_sources_dict2:
      all_sources2[source['alias']] = (source.get('bibcode',''),
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
   photometry = d.get('photometry', [])
   for p in photometry:
      if p.get('upperlimit',False):
         continue
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
         for s in ss:
            if s in all_sources2:
               this_source = all_sources2[s]

      if this_source is None:
         print("Warning:  no primary source, skipping")
         continue

      bibcode = this_source[0]
      if bibcode in pubs:
         b = pubs[bibcode](t[0],float(p['time']))
      elif t in ftrans:
         b = ftrans[t]
      elif t in ftrans_standard:
         b = ftrans_standard[t]
         if t not in known_unknowns:
            known_unknowns.append(t)
            print("Warning:  no telescope/system info, assuming ", \
                  standard_warnings[b[0]], b[0])
      elif (t[0],"","","") in ftrans_standard:
         b = ftrans_standard[(t[0],"","","")]
         if t not in known_unknowns:
            known_unknowns.append(t)
            print("Warning: telescope/system defined by %s/%s/%s not "\
                  "recognized, assuming %s %s" %\
                  (t[1],t[2],t[3],standard_warnings[t[0]],t[0]))
      else:
         # No idea
         if t not in unknown_unknowns:
            unknown_unknowns.append(t)
            print("Warning: telescope/system defined by %s/%s/%s not "\
                  "recognized and can't figure out the filter %s" % \
                  (t[1],t[2],t[3],t[0]))
         unknown_unknowns.append(t)
         continue
      if b not in MJD:  
         MJD[b] = []
         mags[b] = []
         emags[b] = []
         sids[b] = []
      if 'time' in p and 'magnitude' in p:
         if not allow_no_errors and 'e_magnitude' not in p and\
               'e_lower_magnitude' not in p and 'e_upper_magnitude' not in p:
            if 'upperlims' not in warnings: warnings.append('upperlims')
            continue
         MJD[b].append(float(p['time']))
         mags[b].append(float(p['magnitude']))
         if 'e_magnitude' in p:
            emags[b].append(float(p['e_magnitude']))
         elif 'e_lower_magnitude' in p and 'e_upper_magnitude' in p:
            emags[b].append((float(p['e_lower_magnitude']) +\
                  float(p['e_upper_magnitude']))/2)
         else:
            emags[b].append(missing_error)
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
      print("Warning:  the following photometry was not recognized by SNooPy")
      print("and was not imported:")
      for item in unknown_unknowns:
         print(item)
   if warnings:
      for warning in warnings:
         print(warning_message[warning])

   # lastly, the spectroscopy
   if d.get('spectra',None) is not None:
      spectra = []
      dates = []
      sids = []
      for s in d['spectra']:
         wu = s.get('u_wavelengths', 'Agnstrom')
         fu = s.get('u_fluxes', 'Uncalibrated')
         
         try:
            wu = u.Unit(wu)
         except ValueError:
            print("Warning:  unrecognized unit for wavelength: {}".format(wu))
            print("  assuming Angstroms")
            wu = u.Angstrom
 
         if fu == 'Uncalibrated':
            fluxed = False
            fu = u.dimensionless_unscaled
         else:
            try:
               fu = u.Unit(fu)
               fluxed = True
            except ValueError:
               print("Warning:  unrecognized unit for flux: {}".format(fu))
               fluxed = False
               fu = u.dimensionless_unscaled
 
         tu = s.get('u_time', 'MJD')
         t = float(s['time'])
         if tu not in MJD_offsets:
            print("Warning:  unrecognized time unit: {}".format(tu))
            if len(s['time'].split('.')[0]) == 7 and s['time'][0] == '2':
               print("   assuming JD")
               t = t - 2400000.5
            elif len(s['time'].split('.')[0]) == 5 and s['time'][0] == '5':
               print("   assuming MJD")
            else:
               print("   skipping this spectrum.")
               continue
 
         w = array([float(item[0]) for item in s['data']])*wu
         f = array([float(item[1]) for item in s['data']])*fu
         dr = s.get('deredshifted', False)
         if dr:
            w = w*(1+zhel)
 
         # At this point, we should be able to convert to the units we want
         w = w.to('Angstrom').value
         if fluxed: f = f.to('erg / (s cm2 Angstrom)')
         f = f.value
 
         # source reference
         srcs = s.get('source','').split(',')
         this_source = None
         for src in srcs:
            if src in all_sources:
               this_source = all_sources[src]
               break
         if this_source is None:
            print("Warning: spectrum has no source")
 
         if this_source not in used_sources:
            used_sources.append(this_source)
         # At this point we're actually using the spectroscopy, so find source
         sid = used_sources.index(this_source)
         sids.append(sid)
         
         spectra.append(spectrum(wave=w, flux=f, fluxed=fluxed, 
            name="Spectrum MJD={:.1f}".format(t)))
         dates.append(t)
      snobj.sdata = timespec(snobj, dates, spectra)
      snobj.sdata.sids = sids

   if full_data:
      # make a dictionary of the remaining OSC meta data and make it a member 
      # variable
      snobj.osc_meta = {}
      for key in d.keys():
         if key not in ['name','redshift','ra','dec','sources','photometry',
               'spectra']:
            snobj.osc_meta[key] = d[key]

   return(snobj,'Success')

def to_osc(s, ref=None, bibcode=None, source=1):
   '''Given a supernova object, s, output to JSON format suitable for upload to
   the OSC.'''
   
   data = {s.name:{"name":s.name}}

   if ref or bibcode:
      sources = [dict(bibcode=bibcode, name=ref, alias=str(source))]
      data['sources'] = sources
   phot = []
   for filt in s.data:
      for i in range(len(s.data[filt].MJD)):
         datum = dict(survey='CSP', observatory='LCO')
         datum.update(CSP_systems(filt, s.data[filt].MJD[i]))
         datum['time'] = "{:.3f}".format(s.data[filt].MJD[i])
         datum['u_time'] = "MJD"
         datum['magnitude'] = "{:.3f}".format(s.data[filt].mag[i])
         flux,eflux = s.data[filt].flux[i],s.data[filt].e_flux[i]
         datum['flux'] = "{:.5f}".format(flux)
         datum['u_flux'] = "s^-1 cm^-2"
         datum['e_flux'] = "{:.5f}".format(eflux)
         datum['e_upper_magnitude'] = "{:.3f}".format(
               -2.5*log10((flux-eflux)/flux))
         datum['e_lower_magnitude'] = "{:.3f}".format(
               -2.5*log10(flux/(flux+eflux)))
         datum['source'] = "{}".format(source)
         phot.append(datum)
   data['photometry'] = phot

   return json.dumps(data, indent=4)

