'''Along with photometric data, let's create a spectral sequence that can
be included in a SNooPy object.'''

from numpy import *
from snpy import fset
from snpy.filters import ch
from snpy.filters import spectrum
from snpy.filters import filt
from matplotlib import pyplot as plt
import scipy

def smoothspec(spec, boxsize):
   '''Do a boxcar average of a spectrum.'''
   # copy mirror-images on the ends to avoid transient signals
   s = r_[spec[boxsize-1:0:-1],spec,spec[-2:-boxsize-1:-1]]
   # window function is a box
   w = ones(boxsize, 'd')
   return convolve(w/w.sum(), s, mode='valid')

class timespec:
   '''This class defines a time series of spectra. It includes the time
   of each spectrum as MJD, and a list of spectrum objects (from the filters
   subpackage). We should be able to do things like plot spectra, fit lines,
   and do synthetic light-curves. It is meant to be contained by the
   :class:`snpy.sn` class.
   
   Args:
      parent (snpy.sn instance), the parent SNooPy sn object.
      MJD (float array): date of observation. While we don't enforce MJD
                         as the epoch zero-point, it should be consistent with
                         the photometry (snpy.lc).
      spectra (list of snpy.filters.spectum objects): the spectra.
      '''

   def __init__(self, parent, MJD, spectra):
      
      self.parent = parent
      self.MJD = asarray(MJD)
      if not isinstance(spectra, list):
         raise ValueError("spectra must be a list of snpy.filters.spectrum")
      for spec in spectra:
         if not isinstance(spec, spectrum):
            raise ValueError("spectra must be a list of snpy.filters.spectrum")
      self.spectra = spectra

   def __getitem__(self, i):
      return self.spectra.__getitem__(i)

   def sort_time(self):
      '''If the time sequence of the spectra are not in ascending order, shuffle
         them accordingly.'''
      sids = argsort(self.MJD)
      self.spectra = [self.spectra[i] for i in sids]

   def dayindex(self, day):
      '''Returns the index of the data corresponding to day.'''
      id = argmin(absolute(self.MJD - day))
      return(id)
   
   def plot(self, day=None, idx=None, ax=None, **popts):
      xlabel = popts.get('xlabel', "Wavelength (Angstroms)")
      ylabel = popts.get('ylabel', "Flux (erg/s/cm/A)")

      if ax is None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      if idx is None:
         if day is None:
            idx = 0
         else:
            idx = self.dayindex(day)
      ax.plot(self.spectra[idx].wave, self.spectra[idx].flux, '-')
      ax.set_xlabel(xlabel)
      ax.set_ylabel(ylabel)
      return(ax)

   def synth_lc(self, band, z=0, zeropad=False):
      '''Produce a synthetic light-curve for filter ``band`` with an
      optional redshift (ie, the spectrum is redshifted before computing).
      If a spectrum doesn't cover the filter definition, it will not be
      included in the light-curve.
      
      Args:
         band (str or filters.filt instance): The filter to compute.
         z (float):  Redshift the spectrum by this much before computing
                     light-curve
         zeropad (bool): Assume the spectrum has zero flux outside the
                     filter definitions (force an answer)
                     
      Returns:
         (MJD, mags):  MJD(float array): epochs of the light-curve
                       mags(floag array): synthetic magnitudes
      '''
      if not isinstance(band, filt) and not isinstance(band, basestring):
         raise ValueError, "band must be string or filters.filt instance"
      if isinstance(band, basestring):
         band = fset[band]
      mags = array([band.synth_mag(spec,z=z,zeropad=zeropad) \
            for spec in self.spectra])
      gids = ~isnan(mags)
      gids = gids * array([isinstance(spec.flux, u.quantity.Quantity) 
         for spec in self.spectra])
      return (self.MJD[gids], mags[gids])

   def synth_emag(self, filt, days, z=0):
      '''Produce the error in a synthetic magnitude based on errors in the 
      spectra.  This assumes errors in the spectra are not correlated... 
      probably not right.'''
      if type(filt) is type(""):
         filt = fset[filt]
      elif type(filt) is not type(fset['Bs']):
         raise TypeError("Error:  filter must be a string or filter object")
      fwave = filt.wave/(1.0 + z)
      fresp = filt.resp

      emags = []
      for day in days:
         id = self.dayindex(day)
         # re-sample the filter on the spectrum resolution
         tck = scipy.interpolate.splrep(fwave,fresp,0*fwave+1, k=1, s=0)
         fresp2 = scipy.interpolate.splev(self.waves[id], tck)
         delta_wave = average(self.waves[id][1:] - self.waves[id][:-1])
         F = filt.response(self.waves[id], self.spectra[id], z)
         varF = filt.response(self.waves[id], 
               power(self.errors[id],2)*fresp2*self.waves[id]/ch*delta_wave, 
               photons=1)
         em = 1.086*sqrt(varF)/F
         emags.append(em)
      return(array(emags))


