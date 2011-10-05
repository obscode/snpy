# A python object that holds spectra for destiny observations

from numpy import *
from pygplot import *
from snpy import fset
from snpy.filters import ch
import scipy


class timespec:

   def __init__(self, parent, MJD, waves, fluxes, e_fluxes):
      self.parent = parent
      self.MJD = asarray(MJD)
      self.waves = waves
      self.fluxes = spectra
      self.e_fluxes = errors

   def sort_time(self):
      '''If the time sequence of the spectra are not in ascending order, shuffle
         them accordingly.'''
      sids = argsort(self.MJD)
      self.MJD = take(self.MJD, sids)
      self.spectra = [self.spectra[i] for i in sids]
      self.waves = [self.waves[i] for i in sids]
      self.errors =  [self.errors[i] for i in sids]

   def dayindex(self, day):
      '''Returns the index of the data corresponding to day.'''
      id = argmin(absolute(self.MJD - day))
      return(id)
   
   def plot(self, day, p=None, **popts):
      if 'xlabel' not in popts: popts['xlabel'] = 'Wavelength (Angstroms)'
      if 'ylabel' not in popts: popts['ylabel'] = 'Flux (erg/s/cm/A)'

      if p is None:
         p = Plot(**popts)
      id = self.dayindex(day)
      p.line(self.waves[id], self.spectra[id])
      p.line(self.waves[id], self.errors[id], color='red')
      p.plot()
      p.close()
      return(p)

   def synth_mag(self, filter, days, z=0):
      '''Produce a synthetic magnitude for filter on the days specified, with an
      optional redshift (ie, the filter is redshifted).'''
      if type(filter) is type(""):
         filter = fset[filter]
      elif type(filter) is not type(fset['Bs']):
         raise TypeError, "Error:  filter must be a string or filter object"
      mags = []
      for day in days:
         id = self.dayindex(day)
         mags.append(filter.synth_mag(self.waves[id], self.spectra[id],z))
      mags = array(mags)
      return(array(mags))

   def synth_emag(self, filter, days, z=0):
      '''Produce the error in a synthetic magnitude based on errors in the spectra.
      This assumes errors in the spectra are not correlated... probably not right.'''
      if type(filter) is type(""):
         filter = fset[filter]
      elif type(filter) is not type(fset['Bs']):
         raise TypeError, "Error:  filter must be a string or filter object"
      fwave = filter.wave/(1.0 + z)
      fresp = filter.resp

      emags = []
      for day in days:
         id = self.dayindex(day)
         # re-sample the filter on the spectrum resolution
         tck = scipy.interpolate.splrep(fwave,fresp,0*fwave+1, k=1, s=0)
         fresp2 = scipy.interpolate.splev(self.waves[id], tck)
         delta_wave = average(self.waves[id][1:] - self.waves[id][:-1])
         F = filter.response(self.waves[id], self.spectra[id], z)
         varF = filter.response(self.waves[id], 
               power(self.errors[id],2)*fresp2*self.waves[id]/ch*delta_wave, 
               photons=1)
         em = 1.086*sqrt(varF)/F
         emags.append(em)
      return(array(emags))

   def synth_lc(self, filter, z=0):
      '''Produce a synthetic light-curve for the filter using all the days in the
      spectrum object.  Returns (days,mags)'''
      mags = self.synth_mag(filter, self.days, z)
      emags = self.synth_emag(filter, self.days, z)
      return(self.days, mags, emags)

