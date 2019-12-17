'''Along with photometric data, let's create a spectral sequence that can
be included in a SNooPy object.'''

from numpy import *
from snpy import fset
from snpy.filters import ch
from snpy.filters import spectrum
from snpy.filters import filter
from matplotlib import pyplot as plt
from matplotlib.widgets import Button,CheckButtons
import scipy
import six

cyc = plt.rcParams['axes.prop_cycle']
colors = cyc.by_key()['color']

elements = {
      'H':(656.45377, 4861.3615, 4340.462, 4101.74),
      'NaI':(5889.950,5895.924,8194.824),
      'CaII':(3933,3968,8498,8542,8662),
      'SiII':(4130, 5972, 6355)}
elnames = elements.keys()

def smoothspec(spec, boxsize):
   '''Do a boxcar average of a spectrum.'''
   # copy mirror-images on the ends to avoid transient signals
   s = r_[spec[boxsize-1:0:-1],spec,spec[-2:-boxsize-1:-1]]
   # window function is a box
   w = ones(boxsize, 'd')
   return convolve(w/w.sum(), s, mode='valid')[boxsize/2:-boxsize/2+1]

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

   def __getstate__(self):
      odict = self.__dict__.copy()
      return odict

   def __setstate__(self, state):
      self.__dict__.update(state)


   def sort_time(self):
      '''If the time sequence of the spectra are not in ascending order, shuffle
         them accordingly.'''
      sids = argsort(self.MJD)
      self.MJD = self.MJD[sids]
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
         band (str or filters.filter instance): The filter to compute.
         z (float):  Redshift the spectrum by this much before computing
                     light-curve
         zeropad (bool): Assume the spectrum has zero flux outside the
                     filter definitions (force an answer)
                     
      Returns:
         (MJD, mags):  MJD(float array): epochs of the light-curve
                       mags(floag array): synthetic magnitudes
      '''
      if not isinstance(band, filter) and not \
            isinstance(band, six.string_types):
         raise ValueError("band must be string or filters.filter instance")
      if isinstance(band, six.string_types):
         band = fset[band]
      mags = array([band.synth_mag(spec,z=z,zeropad=zeropad) \
            for spec in self.spectra])
      gids = ~isnan(mags)
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

class InteractiveSpec:

   def __init__(self, series):
      self.series = series
      self.curID = 0
      self.fig = None
      self.ax = None

      self._bindings = {}
      self._binding_help = {}

      self._box = 1

      # Lists of artists that need to be cleared each cycle
      self._pline = []
      self._pfill = []
      self._l = []

      self.initplot()

   def initplot(self):
      self.fig,self.ax = plt.subplots(figsize=(10,5))
      self.ax.set_xlabel('Wavenlength (Angstroms)')
      self.ax.set_ylabel('Flux (erg/s/cm^2/Angstrom)')
      plt.tight_layout()
      plt.subplots_adjust(right=0.7)

      labs = ['|<','<','>','>|']
      binds = [self._first, self._previous, self._next, self._last]
      x0,x1,y0,y1 = 0.75,0.95,0.9,0.95
      dx = (x1-x0)/4
      self._buttons = []
      for i in range(4):
         butax = plt.axes([x0 + dx*i, y0, dx, y1-y0])
         self._buttons.append(Button(butax, labs[i]))
         self._buttons[-1].on_clicked(binds[i])

      #for el in elements:
      #   for line in elements[el]:

      cax = plt.axes([x0, 0.5, x1-x0, 0.3])
      arg = tuple([False for el in elnames])
      self._welements = CheckButtons(cax, elnames, arg)
      for i in range(len(self._welements.labels)):
         self._welements.labels[i].set_color(colors[i])

      self._welements.on_clicked(self._elements)

      self.plot(0)

      self.set_bindings()

      # The element lines
      self._ellines = {}
      self._ellines_z = {}
      z = self.series.parent.z
      for i,el in enumerate(elnames):
         self._ellines[el] = []
         self._ellines_z[el] = []
         for li in elements[el]:
            self._ellines[el].append(
               self.ax.axvline(li, color=colors[i], linestyle='--', 
                  visible=False, alpha=0.5))
            self._ellines_z[el].append(
               self.ax.axvline(li*(1+z), color=colors[i], visible=False,
                  alpha=0.5))

   def plot(self, ID, clear=True):
      '''Plot the spectrum number ID.'''
      if ID < 0 or ID > len(self.series.spectra)-1:
         # Out of range
         return

      wave,flux = self._getcurspec()
      MJD = self.series.MJD[ID]
      label = '{:.1f}'.format(MJD)
      if self.series.parent.Tmax > 0:
         label += " ({:.1f})".format(MJD-self.series.parent.Tmax)
      if clear:
         ll = self._l + self._pline + self._pfill
         for l in ll:
            l.remove()
         self.ax.set_prop_cycle(None)   # reset colors/styles, etc
         self._l = self.ax.plot(wave, flux, '-', label=label)
      else:
         self._l += self.ax.plot(wave, flux, '-', label=label)

      x0,x1 = wave.min()*0.95, wave.max()*1.05
      y0,y1 = flux.min()*0.95, flux.max()*1.05
      self.ax.set_xlim(x0,x1)
      self.ax.set_ylim(y0,y1)
      self.ax.legend(loc='upper right')

   def _first(self, event):
      self.curID = 0
      self.plot(0)

   def _last(self, event):
      self.curID = len(self.series.spectra)-1
      self.plot(self.curID)

   def _next(self, event):
      if self.curID == len(self.series.spectra)-1:
         return
      self.curID += 1
      self.plot(self.curID)

   def _previous(self, event):
      if self.curID == 0:
         return
      self.curID -= 1
      self.plot(self.curID)

   def set_bindings(self):
      # get rid of key bindings
      ids = list(self.fig.canvas.callbacks.callbacks['key_press_event'].keys())
      for id in ids:
         self.fig.canvas.mpl_disconnect(id)
      self.bind_key('q', lambda event: plt.close(self.fig), 'quit viewer')
      self.bind_key('?', self.help, 'Get help on key bindings')
      self.bind_key('p', self._pequiv, 'Measure pseudo-equivalent width')
      self.bind_key('s', self._smooth, 'Decrease smoothing')
      self.bind_key('S', self._smooth, 'Increase smoothing')

      self._key = self.fig.canvas.mpl_connect('key_press_event', 
            self._key_press)

   def _key_press(self, event):
      # portal callback that hands off to registered events
      if event.key in self._bindings:
         self._bindings[event.key](event)

   def help(self, event):
      '''produce a help message with registered bindings.'''
      keys = self._bindings.keys()
      keys.sort()
      print("Available key commands:")
      for key in keys:
         print("{:5s}: {}".format(key, self._binding_help[key]))

   def bind_key(self, key, callback, helpmsg):
      '''Bind the function [callback] when [key] is pressed. [helpmsg]
      is the help message to pring when self.help() is called.'''
      if not isinstance(key, six.string_types) or len(key) != 1:
         raise TypeError("key must be a single character")
      if not callable(callback):
         raise TypeError("callback must be callable")
      if not isinstance(helpmsg, six.string_types):
         raise ValueError("callback must be a string")
      self._binding_help[key] = helpmsg
      self._bindings[key] = callback

   def _getcurspec(self):
      '''Get the current spectrum, smoothed if necessary.'''
      spec = self.series.spectra[self.curID]
      x,y = spec.wave, spec.flux
      if self._box > 1:
         y = smoothspec(y, self._box)
      return x,y

   def _pequiv(self, event):
      # the data we're working with
      spec = self.series.spectra[self.curID]
      w,f = event.xdata, event.ydata

      wave,flux = self._getcurspec()
      dists = power(wave-w,2)+power(flux-f,2)
      id = argmin(dists)
      if getattr(self, '_pstate', None) is None:
         # First point on the plot
         # Find the closest point
         self._pstate = (wave[id], flux[id])
         print("Choose second point")
      else:
         # Second point on the plot
         x0,y0 = self._pstate
         x1,y1 = wave[id],flux[id]
         if x0 > x1:
            x0,x1 = x1,x0
            y0,y1 = y1,y0
         
         self._pline += self.ax.plot([x0, x1],[y0,y1], '-')
         gids = greater(wave,x0)*less(wave,x1)
         yy = lambda x: y0 + (y1-y0)/(x1-x0)*(x-x0)
         self._pfill.append(self.ax.fill_between(wave[gids], flux[gids],
               yy(wave[gids]), facecolor=self._pline[-1].get_color(), 
               edgecolor='none', alpha=0.5))
         pcont = scipy.integrate.trapz(yy(wave[gids]), x=wave[gids])
         integ = scipy.integrate.trapz(flux[gids], x=wave[gids])
         line = pcont - integ
         imin = argmin(flux[gids])
         xmin = spec.wave[gids][imin]

         pequiv = line/yy(xmin)
         print("Pseudo-continuum:  ", pcont)
         print("Line area:", line)
         print("Pseudo-equivalent width:", pequiv)
         self._pstate = None
   
   def _smooth(self, event):
      if event.key == 's':
         # decrease smoothing
         if self._box == 1:
            return
         self._box -= 2
      else:
         self._box += 2
      wave,flux = self._getcurspec()
      self._l[-1].set_ydata(flux)

   def _elements(self, label):
      if label not in self._ellines:
         return
      for line in self._ellines[label]:
         line.set_visible(not line.get_visible())
      for line in self._ellines_z[label]:
         line.set_visible(not line.get_visible())
      plt.draw()
