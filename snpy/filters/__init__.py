'''This module provices two classes:  spectrum and filter (filter is a
sub-class of spectrum.  These classes are designed to make certain tasks more
convenient, especially for filters.  Aside from being containers of the
spectrum/filter data, they provide the following funcionality:

   spectrum:
      - read a wavelength-flux two-column file.
      - make safe copies of the data
      - several member variables for bookkeeping:
         o max/min wavelenghts
         o name, coment
   filter:
      - everthing in spectrum, plus:
      - computation of effective wavelength for a given spectrum
      - compute zero-point of a filter based on reference spectrum
        and supplied magnitude of the reference
      - compute the convolution of the filter with a spectrum.
      - compute a synthetic magnitude based on a supplied filter.

This module also supplies a dictionary called filters with some common filters
and zero-points built in.'''

from __future__ import print_function
import os,sys,types
import numpy as num
import scipy.integrate
import scipy.interpolate
from glob import glob
import string
from snpy.utils.deredden import unred
#from astropy import units as u
from astropy.constants import c,h

interp_method = 'spline'
integ_method = 'simpsons'
subsample = 1

base = os.path.abspath(globals()['__file__'])
base = os.path.dirname(base)
filter_base = os.path.join(base,'filters')
stand_base = os.path.join(base,'standards')

#h = 6.626068e-27  # erg s
#c = 2.997925e18   # Angstrom/s
ch = c * h
ch = ch.to('erg Angstrom').value

class spectrum:
   '''This class defines a spectrum.  It contains the response as Numeric 
   arrays.  It has the following member data:
      name:      string describing the filter (eg. 'B')
      file:      where to find the filter response
      wave:      Numeric array of wavelengths (Angstroms)
      resp:      Numeric array of response
      min:       minimum wavelength defined
      max:       maximum wavelength defined
      comment:   They're useful, you know!
    
    There are also some useful functions:
      read():                   Read in the data and compute member data
   '''
   def __init__(self, name=None, filename=None, comment=None, wave=None,
         flux=None, fluxed=True, load=1):
      '''Creates a spectrum instance.  Required parameters:  name and file.  
      Can also specify the zero point (instead of using the comptute_zpt() 
      function do do it).'''
      self.name = name
      self.file = filename     # location of the filter response
      self.wave_data = wave    # wavelength of response
      self.resp_data = flux    # response
      self.comment = comment   # any words?
      self.fluxed = fluxed     # indicates the spectrum is in physical units
      if filename is not None and load==1:  self.read()

   def __str__(self):
      return "%s:  %s" % (self.name, self.comment)
   def __repr__(self):
      return "%s:  %s" % (self.name, self.comment)

   def read(self):
      '''Reads in the response for file and updates several member functions.'''
      if self.file is not None:
         f = open(self.file)
         lines = f.readlines()
         self.wave_data = num.array([float(line.split()[0]) \
               for line in lines if line[0] != "#"])
         self.resp_data = num.array([float(line.split()[1]) \
               for line in lines if line[0] != "#"])
         f.close()

   def copy(self):
      return(spectrum(self.name, self.file))

   def waverange(self):
      if self.wave is not None:
         return(num.minimum.reduce(self.wave), num.maximum.reduce(self.wave))
      else:
         return(None)

   def __getattr__(self, name):
      if name in ['wave','resp','flux'] and self.wave_data is None:
         self.read()
      if name == 'wave':
         return self.wave_data
      elif name == 'resp':
         return self.resp_data
      elif name == 'flux':
         return self.resp_data
      elif name == "wavemax":
         return self.waverange()[1]
      elif name == "wavemin":
         return self.waverange()[0]
      elif name == "ave_wave":
         if self.wave is not None and self.resp is not None:
            return num.sum(self.wave*self.resp)/num.sum(self.resp)
         else:
            return None
      else:
         raise AttributeError("Error:  attribute %s not defined" % (name))


class filter(spectrum):
   '''This class defines a filter.  It contains the response as Numeric arrays.  It has
   the following member data:
      name:      string describing the filter (eg. 'B')
      file:      where to find the filter response
      zp:        the photometric zero point (in vega only for now)
      wave:      Numeric array of wavelengths (Angstroms)
      resp:      Numeric array of response
      ave_wave:  The effective wavelength for a flat spectrum
      min:       minimum wavelength defined
      max:       maximum wavelength defined
      comment:   They're useful, you know!
    
    There are also some useful functions:
      read():                   Read in the data and compute member data
      compute_zpt(spec, mag):   Compute zero point based on spectrum instance spec
                                and standard magnitude mag.
      response(wave, flux):     Compute the convolution with spectrum given by
                                wave,flux
      eff_wave(wave, flux):     Compute the effective wavelength given a spectrum
      synth_mag(wave, flux):    Compute a synthetic magnitude, given a spectrum.
      R(Rv, z=0, strict_ccm=0): Compute the filters' Reddening coefficient for 
                                assumed value of Rv and redshift z.'''

   def __init__(self, name, file=None, zp=None, comment=None):
      '''Creates a filter instance.  Required parameters:  name and file.  Can 
      also specify the zero point (instead of using the comptute_zpt() 
      function do do it).'''
      spectrum.__init__(self, name, file)
      self.zp = zp
      self.comment = comment
      self.read()
      self.tck = None     # Used for interpolating the filter response
      self.mint = None    #    "

   def read(self):
      '''Reads in the response for file and updates several member functions.'''
      spectrum.read(self)

   def compute_zpt(self, spectrum, mag, zeropad=0):
      '''Compute the photometric zero point.  If spectrum is a list of spectra, then
      a zero point is computed for each and returned as a Numeric array (which you
      can then average, median, whatever.'''
      # get the response if needed:
      if self.wave is None: self.read()

      if type(spectrum) is not list:
         spectrum = [spectrum]
         only1 = 1
      else:
         only1 = 0

      zpts = []
      for spec in spectrum:
         # Check to see if 
         if spec.wave is None:  spec.read()
         if not spec.fluxed:
            raise ValueError("spectrum must be in erg/s/cm^2/Angstrom")
   
         # Compute the integral spec1*spec2*(lambda/ch):
         result = self.response(spec, zeropad=zeropad)

         # Now use the spectrum's magnitude to compute zero point:
         zpt = 2.5*num.log10(result) + mag
         zpts.append(zpt)
   
      if only1:
         return(zpts[0])
      else:
         return(zpts)

   def eval(self, wave):
      '''evaluate the filter on the sequence of wavelengths.'''
      if interp_method == "spline":
         if self.tck is None:
            self.tck = scipy.interpolate.splrep(self.wave, 
                  self.resp, k=1, s=0)
         fresp_int = scipy.interpolate.splev(wave, self.tck)
      else:
         if self.mint is None:
            self.mint = scipy.interpolate.interp1d(self.wave, 
                  self.resp, kind=interp_method)
         fresp_int = self.mint(trim_wave)

      fresp_int = num.where(num.less(wave,self.wave.min()), num.nan, fresp_int)
      fresp_int = num.where(num.greater(wave,self.wave.max()), num.nan, 
                  fresp_int)
      return(fresp_int)

   def response(self, specwave, flux=None, z=0, zeropad=0, photons=1):
      '''Get the response of this filter over the specified spectrum.  This
      spectrum can be defined as a spectrum instance, in which case you simply
      need to specify [specwave]  Or, you can specify a wavelength
      and flux vector, in which case, you need to specify both [specwave] (which
      is now taken to be the wavelength vector) and the flux as [flux].  If z is
      supplied, first redshift the spectrum by this amount.  If zeropad is true,
      then the spectrum is assumed to be zero where the filter extends beyond
      its definition, otherwise -1 is returned if the filter extends beyond the
      spectrum's definition.  If photons=1, the integrand is multiplied by the
      wavelength vector and divided by c*h, i.e., the photon flux is
      computed..'''

      # Handle the intput parameters
      if flux is None:
         # We must have a spectrum object:
         if not isinstance(specwave, spectrum):
            raise TypeError("If specifying just specwave, it must be a spectrum object")
         # if this object is not fluxed, then return -1
         if not getattr(specwave, 'fluxed', True):
            return -1.0
         wave = specwave.wave
         spec = specwave.flux
      else:
         if type(specwave) is not num.ndarray or type(flux) is not num.ndarray:
            raise TypeError("If specifying both specwave and flux, they must be arrays")
         if len(num.shape(specwave)) != 1 or len(num.shape(flux)) != 1:
            raise TypeError("specwave and flux must be 1D arrays")
         wave = specwave
         spec = flux

      if z > 0:
         swave = wave*(1.+z)
      elif z < 0:
         swave = wave/(1.+z)
      else:
         swave = wave
      if (self.wavemin < swave[0] or self.wavemax > swave[-1]) and not zeropad:
            return(-1.0)

      # Now figure out the limits of the integration:
      x_min = num.minimum.reduce(self.wave)
      x_max = num.maximum.reduce(self.wave)
      try:
         i_min = num.nonzero(num.greater(swave - x_min, 0))[0][0]
      except:
         i_min = 0
      try:
         i_max = num.nonzero(num.greater(swave - x_max, 0))[0][0]
      except:
         i_max = len(swave)-1
   
      if i_min >= 5:
         i_min -= 5
      else:
         i_min = 0
      if i_max <= len(swave)-6:
         i_max += 5
      else:
         i_max = len(swave) - 1

      trim_spec = spec[i_min:i_max+1:subsample]
      trim_wave = swave[i_min:i_max+1:subsample]
      # Now, we need to resample the response wavelengths to the spectrum:
      if interp_method == "spline":
         if self.tck is None:
            self.tck = scipy.interpolate.splrep(self.wave, 
                  self.resp, k=1, s=0)
         fresp_int = scipy.interpolate.splev(trim_wave, self.tck)
      else:
         if self.mint is None:
            self.mint = scipy.interpolate.interp1d(self.wave, self.resp, 
                  kind=interp_method)
         fresp_int = self.mint(trim_wave)
      # Zero out any places beyond the definition of the filter:
      fresp_int = num.where(num.less(trim_wave, x_min), 0, fresp_int)
      fresp_int = num.where(num.greater(trim_wave, x_max), 0, fresp_int)

      integrand = fresp_int*trim_spec
      if photons:
         integrand = integrand*trim_wave/ch

      if integ_method=='simpsons':
         result = scipy.integrate.simps(integrand, x=trim_wave)
      elif integ_method=='trapz':
         result = scipy.integrate.trapezoid(integrand, x=trim_wave)
      else:
         result = (trim_wave[-1] - trim_wave[0])/(len(trim_wave)-1)*\
               sum(integrand)

      return(result)

   def ABoff(self):
      '''Compute the AB offset for this filter. Due to the way SNooPy stores
      the zero-points, this only depends on filter function shape.'''
      return 65.4469-48.6-self.zp + \
            2.5*num.log10(scipy.integrate.trapezoid(self.flux/self.wave,self.wave))


   def synth_mag(self, specwave, flux=None, z=0, zeropad=0):
      '''Compute the synthetic magnitude based on the input spectrum defined by
      (specwave) or (specwave,flux).  If z is supplied, first redshift the 
      input spectrum by this amount.'''

      # First check to make sure the spectrum is fluxed
      if isinstance(specwave, spectrum):
         if not specwave.fluxed: return(num.nan)
      res = self.response(specwave,flux=flux,z=z,zeropad=zeropad)
      if res <= 0:
         return(num.nan)
      else:
         return(-2.5*num.log10(res) + self.zp)

   def synth_abmag(self, specwave, flux=None, z=0, zeropad=0):
      '''Compute the synthetic AB magnitude of the input spectrum defined by
      (specwave) or (specwave, flux).  If z is supplied, first blueshift
      the filter by this amount (ie, you are observing a redshifed spectrum).'''

      # First check to make sure the spectrum is fluxed
      if isinstance(specwave, spectrum):
         if not specwave.fluxed: return(num.nan)

      numer = self.response(specwave, flux=flux, z=z, zeropad=zeropad,
            photons=1)*ch
      # numer is in erg*Angstrom/s/cm^2
      if numer <= 0:
         return(num.nan)

      if not isinstance(specwave, spectrum):
         # 3631 Jy*c --> erg*Angstrom/s/cm^2
         denom = self.response(specwave, 3631*1.e-23*c/specwave, photons=0)
      else:
         wave = specwave.wave
         denom = self.response(wave, 3631*1.e-23*c/wave, photons=0)

      result = -2.5*num.log10(numer/denom)# - 48.6
      return(result)

   def mag2flux(self, mag, specwave=None, flux=None, z=0):
      '''Convert a magnitude in this filter to the flux in erg/s/cm^2.'''
      if len(num.shape(mag)) == 0:
         scalar = True
         mag = num.array([mag])
      else:
         scalar = False
         mag = num.asarray(mag)
      if specwave is None:
         wave,flux = standards['Vega']['VegaB'].wave,\
                     standards['Vega']['VegaB'].resp
      elif isinstance(specwave, spectrum):
         if not specwave.fluxed:
            raise ValueError("spectrum must be  in erg/s/cm^2/Angstrom")
         wave,flux = specwave.wave,specwave.flux
      else:
         if flux is None:
            raise TypeError("specwave must either be a spectrum instance or "\
                  "an array of wavelengths and flux must be specified")
         wave,flux = specwave,flux
      flam = num.power(10, -0.4*(mag-self.zp))   # in photons/s/cm^2
      flam = flam/self.response(wave, flux, z=z)  # now in ergs/s/cm^2
      
      # now weighted average over the filter
      flam = flam * self.response(wave,flux,z=z,photons=0)
      flam = flam / self.response(wave, wave*0+1, z=z, photons=0)

      if scalar:
         return flam[0]
      return flam

   def eff_wave(self, specwave, flux=None, z=0, zeropad=0):
      '''Compute the effective wavelength for this filter, given the
      spectrum defined by (specwave) or (specwave, flux).  If z is 
      non-zero, first redfhift the spectrum by this amount.'''
      if not isinstance(specwave, spectrum):
         s_wave = specwave
         s_flux = flux
      else:
         s_wave = specwave.wave
         s_flux = specwave.flux

      numer = self.response(s_wave, flux=s_flux*s_wave, z=z, zeropad=zeropad,
            photons=0)
      denom = self.response(s_wave, flux=s_flux, z=z, zeropad=zeropad,
            photons=0)
      if numer <= 0 or denom <=0:
         return(num.nan)
      return(numer/denom)


   def copy(self):
      '''Return a copy of this instance.'''
      return(filter(self.name, self.file, self.zp, self.comment))

   def R(self, Rv=3.1, wave=None, flux=None, z=0.0, EBV=0.001, redlaw='ccm',
         strict_ccm=False):
      '''For a given reddening law Rv (default 3.1), find the ratio of total
      to selective absorption for this filter:  R = A/E(B-V).  You can 
      specify a specific spectrum by supplying a wave and flux and redshift
      (default is defined by filters.reference_wave and 
      filters.refernce_flux at z=0).  You can also specify E(B-V) (EBV) which can
      change the value of R if the spectrum is significantly non-stellar. You
      can specify redlaw='fm' if you prefer a Fitzpatric (1999) reddening 
      law.'''
      global standards

      if wave is None:
         wave,flux = standards['Vega']['VegaB'].wave,\
                     standards['Vega']['VegaB'].resp
      flux0 = self.response(wave, flux, z, photons=1)
      if flux0 <= 0:
         return(num.nan)
      redf = unred(wave, flux, -EBV, Rv, z, redlaw=redlaw, 
            strict_ccm=strict_ccm)[0]
      fluxr = self.response(wave, redf, z, photons=1)
      return(-2.5*num.log10(fluxr/flux0)/EBV)



class system:
   '''An object that contains a photometric system of standards.'''

   def __init__(self, name):
      self.name = name
      self.SEDs = {}

   def add_SED(self, SED):
      if not isinstance(SED, spectrum):
         raise ValueError("SED must be a spectrum instance")
      self.SEDs[SED.name] = SED

   def list_SEDs(self):
      for SED in self.SEDs:
         print("\t"+self.SEDs[SED].name)

   def keys(self):
      return list(self.SEDs.keys())

   def values(self):
      return list(self.SEDs.values())

   def __contains__(self, item):
      return self.SEDs.__contains__(item)

   def __iter__(self):
      return self.SEDs.__iter__()

   def __getattr__(self, attr):
      if attr in self.__dict__['SEDs']:
         return self.__dict__['SEDs'][attr]
      else:
         raise AttributeError

   def __getitem__(self, key):
      if key in self.SEDs:
         return self.SEDs[key]
      else:
         raise AttributeError

   def __str__(self):
      ret = "system %s with standards:  " % (self.name)
      for key in self.SEDs:  ret += "%s, " % (key)
      return ret

   def __repr__(self):
      return self.__str__()

class standard_set:
   '''An object that will contain all the standard SEDs.  The
   standard set contains a dictionary of system objects.  Each
   system object contains a dictionary of spectrum objects.
   So standards.Vega.Bohlin04  would refer to the Bohlin & Gllliand
   2004 SED.  You can also refer to the spectrum with a unique ID 
   as, e.g., standards['VegaB'].'''

   def __init__(self):
      self.systems = {}
      self.spectra = {}

   def add_system(self, name):
      self.systems[name] = system(name)

   def list_systems(self):
      for syst in self.systems:
         print(syst+": "+self.systems[syst].name)

   def list_SEDs(self):
      for syst in self.systems:
         print(self.systems[syst].name)
         self.systems[syst].list_SEDs()

   def cache_spectra(self):
      for syst in list(self.systems.values()):
         for sed in list(syst.SEDs.values()):
            if sed.name in self.spectra:
               print("Warning!  Encountered multiple filter IDs for %s" %\
                        (sed))
            self.spectra[sed.name] = sed

   def __getattr__(self, attr):
      if attr in self.__dict__['systems']:
         return self.__dict__['systems'][attr]
      elif attr in self.__dict__['spectra']:
         return self.__dict__['spectra'][attr]
      else:
         raise AttributeError

   def __getitem__(self, key):
      if key in self.spectra:
         return self.spectra[key]
      elif key in self.systems:
         return self.systems[key]
      else:
         raise KeyError("spectrum ID %s not found" % (key))

   def __setitem__(self, key, value):
      self.spectra[key] = value

   def __contains__(self, key):
      if key in self.spectra:
         return True
      else:
         return False

class filter_set:
   '''An object that will contain all the filter instances.  The
   filter set contains a dictionary of observatory objects.  Each
   observatory object contains a dictionary of telescope objects.
   And each telescope object contains a dictionary of filter
   objects.  So filters.LCO.Swope.B  would refer to the B filter
   on the Swope telescope at the LCO observatory.  You can also
   refer to a filter with a unique ID as, e.g., filters['Bswo'].'''

   def __init__(self):
      self.observatories = {}
      self.filters = {}

   def add_observatory(self, name):
      self.observatories[name] = observatory(name)

   def cache_filters(self):
      for obs in list(self.observatories.values()):
         for tel in list(obs.telescopes.values()):
            for filt in list(tel.filters.values()):
               if filt.name in self.filters:
                  print("Warning!  Encountered multiple filter IDs for %s" %\
                        (filt))
               self.filters[filt.name] = filt

      
   def list_observatories(self):
      for obs in self.observatories:
         print(self.observatories[obs].name)

   def list_telescopes(self):
      for obs in self.observatories:
         print(self.observatories[obs].name)
         self.observatories[obs].list_telescopes()

   def list_filters(self):
      for obs in self.observatories:
         print(self.observatories[obs].name)
         self.observatories[obs].list_filters()

   def __getattr__(self, attr):
      if attr in self.__dict__['observatories']:
         return self.__dict__['observatories'][attr]
      elif attr in self.__dict__['filters']:
         return self.__dict__['filters'][attr]
      else:
         raise AttributeError

   def __getitem__(self, key):
      if key in self.filters:
         return self.filters[key]
      else:
         raise KeyError("filter ID %s not found" % (key))

   def __setitem__(self, key, value):
      self.filters[key] = value

   def __contains__(self, key):
      if key in self.filters:
         return True
      else:
         return False

class observatory:
   '''An object that contains telescope objects.  We could also
   add other info like Lat, long, altitude, etc...'''

   def __init__(self, name):
      self.name = name
      self.telescopes = {}

   def add_telescope(self, name):
      self.telescopes[name] = telescope(name)

   def list_telescopes(self):
      for tel in self.telescopes:
         print("\t"+self.telescopes[tel].name)

   def list_filters(self):
      for tel in self.telescopes:
         print("\t"+self.telescopes[tel].name)
         self.telescopes[tel].list_filters()

   def __getattr__(self, attr):
      if attr in self.__dict__['telescopes']:
         return self.__dict__['telescopes'][attr]
      else:
         raise AttributeError

   def __str__(self):
      ret = "observatory %s with telescopes:  " % (self.name)
      for key in self.telescopes:  ret += "%s, " % (key)
      return ret

   def __repr__(self):
      return self.__str__()

class telescope:
   '''An object that contains filter objects.  It is a child to
   the observatory class, which is a child to the filter_set.'''
   
   def __init__(self, name):
      self.name = name
      self.filters = {}

   def add_filter(self, filter_object):
      if not isinstance(filter_object, filter):
         raise TypeError("Error: filter_object must be a filter type")
      self.filters[filter_object.name] = filter_object

   def list_filters(self):
      for f in self.filters:
         print("\t\t'%s':  %s" % (f, self.filters[f].comment))

   def __getattr__(self, attr):
      if attr in self.__dict__['filters']:
         return self.__dict__['filters'][attr]
      else:
         raise AttributeError

   def __str__(self):
      ret = "telescope %s with filters: " % (self.name)
      for k in list(self.filters.keys()):
         ret += "%s, " % (k)
      return ret

   def __repr__(self):
      return self.__str__()

# Now load in the standard spectra and filter set.
standards = standard_set()
standard_mags = {}
dirs = glob(stand_base+'/*')
dirs = [dir for dir in dirs if os.path.isdir(dir) and \
      os.path.isfile(os.path.join(dir,'standards.dat'))]

for dir in dirs:
   sname = os.path.basename(dir)
   standards.add_system(sname)
   standard_mags[sname] = {}
   f = open(os.path.join(stand_base, sname, 'standards.dat'))
   lines = f.readlines()
   for line in lines:
      if line[0] == "#":  continue
      l = line.split()
      standards[sname].add_SED(spectrum(l[0], 
         os.path.join(stand_base,sname,l[1]), " ".join(l[3:]), load=0))
      standard_mags[sname][l[0]] = {}
      if os.path.isfile(os.path.join(stand_base,sname,l[2])):
         f2 = open(os.path.join(stand_base,sname,l[2]))
         #lines2 = f2.readlines()
         #lines2 = list(map(string.split, lines2))
         lines2 = [line.split() for line in f2.readlines()]
         for i in range(len(lines2)):
            if lines2[i][0] == "#":  continue
            standard_mags[sname][l[0]][lines2[i][0]] = float(lines2[i][1])
         f2.close()
   f.close()
standards.cache_spectra()

vegaK = standards.Vega.VegaK
vegaH = standards.Vega.VegaH01
vegaH85 = standards.Vega.VegaH85
vegaB = standards.Vega.VegaB
vega = vegaB
bd17 = standards.Smith.bd17


fset = filter_set()
obsdirs = glob(os.path.join(filter_base,'*'))
obsdirs = [obs for obs in obsdirs if os.path.isdir(obs)]
for obs in obsdirs:
   obs_name = os.path.basename(obs)
   fset.add_observatory(os.path.basename(obs_name))
   teldirs = glob(os.path.join(filter_base,obs,'*'))
   teldirs = [dir for dir in teldirs if os.path.isdir(dir) and \
         os.path.isfile(os.path.join(dir,'filters.dat'))]
   for dir in teldirs:
      tel_name = os.path.basename(dir)
      fset.observatories[obs_name].add_telescope(tel_name)
      f = open(os.path.join(dir, 'filters.dat'))
      lines = f.readlines()
      for line in lines:
         l = line.split()
         if l[2].find('=') >= 0:
            # We have a std=mag format
            std,mag = [item.strip() for item in l[2].split('=')]
            #std,mag = list(map(string.strip, l[2].split('=')))
            if std  in standards:
               try:
                  m = float(mag)
               except:
                  raise ValueError("Could not convert standard magnitude for filter %s" %\
                        l[0])
               newf = filter(l[0], os.path.join(dir,l[1]), 0.0, 
                  " ".join(l[3:]))
               newf.zp = newf.compute_zpt(standards[std], m)
               fset.observatories[obs_name].telescopes[tel_name].add_filter(newf)
            else:
               raise ValueError("Could not find standard %s for filter %s" % (std,l[0]))


         elif l[2] == 'AB':
            # We have an AB system, so in principle there is no standard. The
            # zero-point is derived from the filter function alone. See
            # documentation.
            newf = filter(l[0], os.path.join(dir,l[1]), 0.0, 
               " ".join(l[3:]))
            newf.zp = 16.84692 + 2.5*num.log10(
                  scipy.integrate.trapezoid(newf.resp/newf.wave, x=newf.wave))
            fset.observatories[obs_name].telescopes[tel_name].add_filter(newf)
         else:
            fset.observatories[obs_name].telescopes[tel_name].add_filter(
                         filter(l[0], os.path.join(dir,l[1]),
                             float(l[2]), " ".join(l[3:])))
      f.close()
fset.cache_filters()
