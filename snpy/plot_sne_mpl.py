'''plot_sne_pg.py

   A collection of plotting routines that SNPY uses to display data.  This
   version uses the matplotlib library.
'''
from numpy import *
import matplotlib
#matplotlib.use('TkAgg')
from snpy import myplotlib
from matplotlib import pyplot,rcParams
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator,AutoLocator
from snpy import kcorr
import types, string
from snpy.filters import fset
from snpy import mangle_spectrum
#from snpy.utils import InteractiveFit

#rcParams['font.family'] = 'serif'
rcParams['font.size'] = 12
dpi = rcParams['figure.dpi']

# Set these so we're never more than this size in pixels.
max_width = 1024
max_height = 768

class color_dict(dict):
   colors = {'B':'blue', 'V':'orange','R':'red','I':'#fa8072','G':'green',
         'U':'purple', 'Y':'#4682b4', 'J':'#008b8b','H':'orange','K':'#8b4513'}
   def __contains__(self, key):
      if key[0].upper() in self.colors:
         return 1
      else:
         return dict.__contains__(self, key)

   def __getitem__(self, key):
      try:
         value = dict.__getitem__(self, key)
         return value
      except KeyError:
         k = key[0].upper()
         if k in self.colors:
            return self.colors[k]
         else:
            return 'k'

class symbol_dict(dict):
   symbols = {'B':'o', 'V':'s','R':'^','I':'D','G':'v','U':'<','Y':'>', 
         'J':'p', 'H':'h', 'K':'*'}
   def __contains__(self, key):
      if key[0].upper() in self.symbols:
         return 1
      else:
         return dict.__contains__(self, key)
   def __getitem__(self, key):
      try:
         value = dict.__getitem__(self, key)
         return value
      except KeyError:
         k = key[0].upper()
         if k in self.symbols:
            return self.symbols[k]
         else:
            return 'o'

# Some default colors and labels to use in plotting.
#default_colors = color_dict()
#default_symbols = symbol_dict()

class click_line:
   '''an interactive chord.  Given an axis instance and a start
   point (x0, y0), draw a dynamic chord that follows the mouse
   until the close() function is called (which returns the 
   coordinates of the final rectangle).'''

   def __init__(self, ax, x0, y0):
      self.ax = ax
      self.x0 = x0
      self.y0 = y0
      self.line = Line2D([x0],[y0], linestyle='-', color='red', marker='x')
      ax.add_artist(self.line)

   def connect(self):
      self.cidmotion = self.line.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

   def on_motion(self, event):
      if event.inaxes != self.line.axes:  return

      self.line.set_data([self.x0, event.xdata], [self.y0, event.ydata])
      self.ax.figure.canvas.draw()

   def close(self):
      self.line.figure.canvas.mpl_disconnect(self.cidmotion)
      res = self.line.get_data()
      self.line.remove()
      self.line.axes.figure.canvas.draw()
      return(res[0][0], res[0][1], res[1][0], res[1][1])



class click_window:
   '''An interactive window.  Given an axis instance and a start point
   (x0,y0), draw a dynamic rectangle that follows the mouse until
   the close() function is called (which returns the coordinates of
   the final rectangle.  Useful or selecting out square regions.'''

   def __init__(self, ax, x0, y0):
      self.ax = ax
      self.x0 = x0
      self.y0 = y0
      self.rect = Rectangle((x0,y0), width=0, height=0, alpha=0.1)
      ax.add_artist(self.rect)

   def connect(self):
      self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

   def on_motion(self, event):
      # Have we left the axes?
      if event.inaxes != self.rect.axes:  return

      self.rect.set_width(event.xdata - self.x0)
      self.rect.set_height(event.ydata - self.y0)
      self.ax.figure.canvas.draw()

   def close(self):
      self.rect.figure.canvas.mpl_disconnect(self.cidmotion)
      extent = self.rect.get_bbox().get_points()
      self.rect.remove()
      self.ax.figure.canvas.draw()
      return(list(ravel(extent)))

class click_xrange:
   '''An interactive xrange selector.  Given an axis and a starting
   x0 location, draw a full-height rectange that follows the mouise.
   Similar to click_window, but more appropriate for selecting out
   an x-range.'''

   def __init__(self, ax, x0):
      self.ax = ax
      self.x0 = x0
      y0,y1 = ax.get_ybound()
      self.rect = Rectangle((x0,y0), width=0, height=(y1-y0), alpha=0.1)
      ax.add_artist(self.rect)

   def connect(self):
      self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

   def on_motion(self, event):
      # Have we left the axes?
      if event.inaxes != self.rect.axes:  return

      self.rect.set_width(event.xdata - self.x0)
      self.ax.figure.canvas.draw()

   def close(self):
      self.rect.figure.canvas.mpl_disconnect(self.cidmotion)
      self.rect.remove()
      self.ax.figure.canvas.draw()
      return(self.x0, self.rect.get_x()+self.rect.get_width())

class click_yrange:
   '''An interactive yrange selector.  Given an axis and a starting
   y0 location, draw a full-width rectange that follows the mouise.
   Similar to click_window, but more appropriate for selecting out
   a y-range.'''

   def __init__(self, ax, y0):
      self.ax = ax
      self.y0 = y0
      x0,x1 = ax.get_xbound()
      self.rect = Rectangle((x0,y0), width=(x1-x0), height=0, alpha=0.1)
      ax.add_artist(self.rect)

   def connect(self):
      self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

   def on_motion(self, event):
      # Have we left the axes?
      if event.inaxes != self.rect.axes:  return

      self.rect.set_height(event.ydata - self.y0)
      self.ax.figure.canvas.draw()

   def close(self):
      self.rect.figure.canvas.mpl_disconnect(self.cidmotion)
      self.rect.remove()
      self.ax.figure.canvas.draw()
      return(self.y0, self.rect.get_y()+self.rect.get_height())


class ButtonClick:
   '''Given a figure instance and an optional dictionary of key-bindings,
   oversee the key bindings in a figure.  There are a set of default
   bindings:  x:  select xrange,  y:  select yrange,  b:  select box
   range, 'm':  mask/unmask data.  These can be overridden using the 
   bindings keyword.  Simply set it to a dictionary with the button as 
   the key and the call-back function as the value.'''

   def __init__(self, figure, bindings={}):
      self.figure = figure
      for ax in self.figure.get_axes():
         ax.xlim_0 = ax.get_xbound()
         ax.ylim_0 = ax.get_ybound()
      self.id = None
      self.id2 = None
      self.pending = None     # Set to axis which has a pending event
      self.pending_key = None
      self.bindings = {}
      self.mouse_bindings = {}
      for key in bindings:
         if key[0:5] == 'mouse':
            self.mouse_bindings[key] = bindings[key]
         else:
            self.bindings[key] = bindings[key]

   def __del__(self):
      self.disconnect()

   def connect(self):
      self.id = self.figure.canvas.mpl_connect('key_press_event', self.keypress)
      if len(list(self.mouse_bindings.keys())) > 0:
         self.id2 = self.figure.canvas.mpl_connect('button_press_event', self.buttonpress)
      return self.id

   def disconnect(self):
      if self.id is not None:
         self.figure.canvas.mpl_disconnect(self.id)
      if self.id2 is not None:
         self.figure.canvas.mpl_disconnect(self.id2)

   def buttonpress(self, event):
      if event.inaxes is None:  return
      if self.pending_key is not None:  return
      if 'mouse'+str(event.button) in self.mouse_bindings:
         return self.mouse_bindings['mouse'+str(event.button)](*(event,))
      #elif event.button == 1:
      #   print "%f %f" % (event.xdata, event.ydata)


   def keypress(self, event):
      #print "key press event:", event.inaxes,event.key
      if event.inaxes is None:  return
      if event.key in self.bindings:
         return self.bindings[event.key](*(event,))

      if self.pending_key is not None:
         if self.pending_key != event.key:  return

      if event.key == 'x':
         # change the xrange
         if self.pending is None:
            # Make a new xrange selector
            self.pending = click_xrange(event.inaxes, event.xdata)
            self.pending.connect()
         else:
            if event.inaxes != self.pending.ax:  return
            x0,x1 = self.pending.close()
            if absolute(x1-x0) == 0:
               self.pending.ax.set_xbound(self.pending.ax.xlim_0)
            else:
               self.pending.ax.set_xbound((x0,x1))
            self.pending = None

      if event.key == 'y':
         # change the xrange
         if self.pending is None:
            # Make a new xrange selector
            self.pending = click_yrange(event.inaxes, event.ydata)
            self.pending.connect()
         else:
            if event.inaxes != self.pending.ax:  return
            y0,y1 = self.pending.close()
            if absolute(y1-y0) == 0:
               self.pending.ax.set_ybound(self.pending.ax.ylim_0)
            else:
               self.pending.ax.set_ybound((y0,y1))
            self.pending = None

      if event.key == 'b':
         # change the xrange
         if self.pending is None:
            # Make a new xrange selector
            self.pending = click_window(event.inaxes, event.xdata, event.ydata)
            self.pending.connect()
         else:
            if event.inaxes != self.pending.ax:  return
            x0,y0,x1,y1 = self.pending.close()
            if absolute(x1-x0) == 0 and absolute(y1-y0) == 0:
               self.pending.ax.set_xbound(self.pending.ax.xlim_0)
               self.pending.ax.set_ybound(self.pending.ax.ylim_0)
            else:
               self.pending.ax.set_xbound((x0,x1))
               self.pending.ax.set_ybound((y0,y1))
            self.pending = None
      if event.key == 'm':
         lc = getattr(event.inaxes, 'lc_inst', None)
         if lc is None:  return
         # Look for the first Line2D that has picker not None
         l = [line for line in event.inaxes.lines if \
               line.get_picker() is not None]
         if not l:  return
         l = l[0]
         xs,ys = l.get_data()
         id = argmin(power(event.xdata-xs,2) + power(event.ydata-ys,2))
         lc.mask[id] = ~lc.mask[id]
         replot_lc(lc)
         return

      if event.key == 'q':
         self.disconnect()
         return

      self.figure.canvas.draw()

def change_SED(event):
   fig = event.inaxes.figure
   ax = event.inaxes
   if event.key == 'd':
      if fig.current_day < 70:  
         fig.current_day += 1
         ax.set_title('Filter responses + %s SED (day %d)' % \
               (fig.version,fig.current_day))
   elif event.key == 'D':
      if fig.current_day > -19:  
         fig.current_day -= 1
         ax.set_title('Filter responses + %s SED (day %d)' % \
               (fig.version,fig.current_day))
   wave,flux = kcorr.get_SED(fig.current_day, version='H3')
   fig.sed_line.set_data(wave*(1+fig.z), flux/flux.max())
   if fig.sed_fill is not None:  
      fig.sed_fill.remove()
      fig.sed_fill = fig.axes[0].fill_between(wave*(1+fig.z), flux/flux.max(),
         facecolor='black', alpha=0.1)
   fig.canvas.draw()



def plot_filters(self, bands=None, day=0, fill=0, outfile=None):
   '''Plot the filt responses, both the observed ones and the restbands
   defined in the sn instance in the rest frame of the SN.  If bands is specified,
   only plot those.  Specify which SED to plot using day.'''
   p = pyplot.figure(112)
   p.clear()
   ax = p.add_subplot(111, 
         title='Filter responses + %s SED (day %d)' % (self.k_version,day),
         xlabel='Wavelength ($\AA$)', ylabel='Relative Flux',
         autoscale_on=False, ylim=(0,1.1))

   p.current_day = day
   p.version = self.k_version
   p.z = self.z

   # Next, for each filt, plot the response and rest band closest to it:
   if bands is None:  bands = list(self.data.keys())
   if type(bands) is bytes:
      bands = [bands]

   minw = 1e12
   maxw = -1
   for band in bands:
      f1 = fset[band]
      f2 = fset[self.restbands[band]]
      minw = min(minw, f1.wave[0], f2.wave[0])
      maxw = max(maxw, f1.wave[-1], f2.wave[-1])
      maxresp = f1.resp.max()
      max_wave = f1.wave[argmax(f1.resp)]
      ax.plot(f1.wave, f1.resp/maxresp, '-', color='blue')
      if fill:
         ax.fill_between(f1.wave, f1.resp/maxresp, facecolor='blue', alpha=0.1)
      ax.annotate(band, (max_wave,1.05), color='blue', 
            horizontalalignment='center', verticalalignment='center')
      maxresp = maximum.reduce(f2.resp)
      max_wave = f2.wave[argmax(f2.resp)]*(1+self.z)
      ax.plot(f2.wave*(1+self.z), f2.resp/maxresp, '-', color='red')
      if fill:
         ax.fill_between(f2.wave*(1+self.z), f2.resp/maxresp, facecolor='red', alpha=0.1)
      ax.annotate(self.restbands[band], (max_wave, 1.05), color='red',
            horizontalalignment='center', verticalalignment='center')

   # first, get the SED:
   sed_wav,sed_flux = kcorr.get_SED(day, self.k_version)
   # normalize:
   sed_flux = sed_flux/sed_flux.max()
   p.sed_line = ax.plot(sed_wav*(1+self.z),sed_flux, scalex=False, color='black')[0]
   if fill:
      p.sed_fill = ax.fill_between(sed_wav*(1+self.z),sed_flux,
            facecolor='black', alpha=0.1)
   else:
      p.sed_fill = None

   ax.set_xbound((minw, maxw))

   pyplot.tight_layout()
   pyplot.draw()
   p.cb = ButtonClick(p, bindings={'d':change_SED, 'D':change_SED})
   p.cb.connect()
   if outfile is not None:
      p.savefig(outfile)
   return p


def plot_SN_panel(obj, ax, filt, delt, symbol, color, Toff, **kwargs):
   '''plot a single SN photometry panel'''
   label = filt
   ax.mylabels = []
   single = kwargs.get('single', False)
   flux = kwargs.get('flux', False)
   msize = kwargs.get('msize', 6)
   linewidth=kwargs.get('linewidth', 1)
   SNR = kwargs.get('SNR_flag', None)

   #if not single:
   #   ax.mylabels.append(ax.text(0.9, 0.9, label, transform=ax.transAxes, 
   #      horizontalalignment='right', fontsize=kwargs.get('fsize', 12), 
   #      verticalalignment='top'))
   if kwargs.get('mask', False):
      x = obj.data[filt].MJD[obj.data[filt].mask]
      if not flux:
         y = obj.data[filt].mag[obj.data[filt].mask] + delt
         ey = obj.data[filt].e_mag[obj.data[filt].mask]
      else:
         y = obj.data[filt].flux[obj.data[filt].mask]*power(10, -0.4*delt)
         ey = obj.data[filt].e_flux[obj.data[filt].mask]
   else:
      x = obj.data[filt].MJD
      if not flux:
         y = obj.data[filt].mag + delt
         ey = obj.data[filt].e_mag
      else:
         y = obj.data[filt].flux*power(10, -0.4*(delt))
         ey = obj.data[filt].e_flux


   label = filt
   if single:
      label = label + '+'+'%.1f' % delt
      
   ax.errorbar(x-Toff, y, yerr=ey, barsabove=True, capsize=0,
         elinewidth=1, fmt=symbol, ms=msize, 
         mfc=color, mec=color, ecolor=color, label=label, linestyle='None')
           #ecolor='black')
   if kwargs.get('label_bad', False):
      gids = equal(obj.data[filt].mask, 0)
      if any(gids):
         x = obj.data[filt].MJD[gids] - Toff
         if not flux:
            y = obj.data[filt].mag[gids] + delt
         else:
            y = obj.data[filt].flux[gids]*power(10, -0.4*(delt))
         ax.plot(x, y, marker='x', mec='red', ms=12, mew=1, linestyle='')
   if not single and SNR is not None:
      cs = ['orange','red']
      for i in range(2): 
         gids = less(obj.data[filt].SNR, SNR[i])
         if any(gids):
            x = obj.data[filt].MJD[gids] - Toff
            if not flux:
               y = obj.data[filt].mag[gids] + delt
            else:
               y = obj.data[filt].flux[gids]*power(10, -0.4*(delt))
            ax.plot(x, y, symbol, mfc=cs[i], ms=msize)


   # Now check to see if there is a model to plot:
   if kwargs.get('plotmodel', True) and obj.model.Tmax is not None and \
         filt in obj.model._fbands:
      t = arange(-10,70,1.0) + obj.Tmax
      mag,err,gids = obj.model(filt, t)
      if any(gids):
         if not flux:
            y = mag + delt
         else:
            zp = fset[filt].zp
            y = power(10, -0.4*(mag - zp + delt))
         ax.plot(compress(gids,t-Toff), compress(gids,y), 
               color='k', linewidth=linewidth)
         l = ax.plot(compress(gids,t-Toff), compress(gids,y+err), 
               '--',color='k', linewidth=linewidth, scaley=False)
         #l[0].autoscale=False
         l = ax.plot(compress(gids,t-Toff), compress(gids,y-err), 
               '--',color='k', linewidth=linewidth, scaley=False)
         #l[0].autoscale=False
   elif obj.data[filt].interp is not None:
      d = obj.data[filt].interp.domain()
      t = arange(d[0], d[1]+1.0, 1.0)
      mag,gids = obj.data[filt].eval(t, t_tol=-1)
      if not flux:
         y = mag + delt
      else:
         zp = fset[filt].zp
         y = power(10, -0.4*(mag - zp + delt))
      ax.plot(t[gids] - Toff, y[gids], color='k',
            linewidth=linewidth)

   #if single and kwargs.get('legend', True):
   ax.legend(loc='upper right', numpoints=1,ncol=2, prop={'size':'small'})


def plot_sn(self, **kwargs):
   '''Plot out the supernova data in a nice format.  There are several 
   options:
      - cols,rows:  Number of columns and rows for panel plot (single=False)
      - xrange,yrange:  specify the ranges to plot as lists [xmin,xmax], 
        [ymin,ymax]
      - title:  optional title
      - interactive:  allows for an interactive plot.
      - single:  plot out as a single (rather than panelled) plot?
      - offset:  automatically offset between the lightcurves (for single plots)
      - fsize:  override the font size used to plot the graphs
      - linewidth:  override the line width
      - symbols:  dictionary of symbols, indexed by band name.
      - colors:  dictionary of colors to use, indexed by band name.
      - relative:  plot only relative magnitudes (normalized to zero)?
      - legend:  do we plot the legend?
      - mask:  Omit plotting masked out data?
      - label_bad:  label the masked data with red x's?
      - Nxticks:  maximum number of x-axis tick marks.
      - Nyticks:  maximum number of y-axis tick marks.
      - JDoffset: If true, compute a JD offset and put it in the x-axis label
                  (useful if x-labels are crowded)
      - SNR_flag: If a tuple, SNR levels to flag in the plot
      - filtsep:  If not none, filters with wavelengths that differ by less
                  than this value are grouped into a single panel (for
                  multi-panel plots).
      - combfilt: If true, combine filters that have the same prefix letter
                  (e.g., g, g2, g_p will all be plotted in the same panel)
      - overplot: specify another SN object to plot along with this one.
                  restbands will be used to match filters if there isn't
                  a one-to-one correspondence.
   '''

   xrange = kwargs.get('xrange', self.xrange)
   yrange = kwargs.get('yrange',self.yrange)
   symbols = kwargs.get('symbols',symbol_dict())
   colors = kwargs.get('colors', color_dict())
   linewidth = kwargs.get('linewidth', 1)
   single = kwargs.get('single', False)
   offset = kwargs.get('offset', False)
   flux = kwargs.get('flux', False)
   oobj = kwargs.pop('overplot', None)
   Nxticks = kwargs.pop('Nxticks', None)
   Nyticks = kwargs.pop('Nxticks', None)
   epoch = kwargs.get('epoch', False)
   JDoff = kwargs.get('JDoffset', False)
   cols = kwargs.get('cols', None)
   rows = kwargs.get('rows', None)
   prunex = kwargs.get('prunex', None)
   pruney = kwargs.get('pruney', None)
   filtsep = kwargs.get('filtsep', None)
   combfilt = kwargs.get('combfilt', False)

   if single: 
      grp = False  # Doesn't make sense
   elif filtsep is not None or combfilt:
      grp = True
   else:
      grp = False

   # See  what filters we're going to use:
   if self.filter_order is not None:
      bands = self.filter_order
   else:
      bands = list(self.data.keys())
      eff_wavs = []
      for filt in bands:
         eff_wavs.append(fset[filt].ave_wave)
      eff_wavs = asarray(eff_wavs)
      ids = argsort(eff_wavs)
      self.filter_order = [bands[i] for i in ids]
      bands = self.filter_order

   # Now let's make a sensible time offset so that we can get rid of
   # overlapping xlabels
   Toff = self.data[bands[0]].MJD.min()
   if epoch:
      if self.Tmax is None:  
         Tmax = 0
      else:
         Tmax = self.Tmax
      Toff = Tmax
   elif JDoff:
      if JDoff is True or JDoff == 'auto':
         Toff = int(min(concatenate([l.MJD for l in list(self.data.values())]))/10.0)
         Toff = Toff*10
      else:
         try:
            Toff = float(JDoff)
         except:
            raise ValueError("JDOffset must be Boolean, 'auto', or number")

   # Figure out logic of where stuff gets plotted. Some filters may be
   # combined because they have the same rest-band (so in some sense are
   # "close").
   if grp:   # group by filter kind
      rbs = []    # list of "rest-bands"
      bidx = {}   # index of panel for each filter
      if combfilt:
         for band in bands:
            if band[0] not in rbs:
               rbs.append(band[0])
            bidx[band] = rbs.index(band[0])
      elif filtsep is not None:
         waves = array([fset[band].ave_wave for band in bands])
         sids = argsort(waves)
         filts = [bands[idx] for idx in sids]
         seps = waves[sids][1:] - waves[sids][:-1]
         rbs.append(filts[0])
         idx = 0
         bidx[filts[0]] = 0
         for i in range(1,len(filts)):
            if seps[i-1] > filtsep:
               rbs.append(filts[i])
               idx += 1
            bidx[filts[i]] = idx

      else:
         for band in bands:
            if band not in self.restbands:
               rbs.append(band)
               bidx[band] = rbs.index(band)
            else:
               if self.restbands[band] not in rbs:
                  rbs.append(self.restbands[band])
               bidx[band] = rbs.index(self.restbands[band])
      n_plots = len(rbs)
   else:
      bidx = {}
      for i in range(len(bands)):
         bidx[bands[i]] = i
      n_plots = len(bands)

   for b in bands:
      if not single and not grp:
         colors[b] = 'blue'
         symbols[b] = 'o'
      elif not single and grp:
         colors[b] = None
         symbols[b] = 'o'
      else:
         if b not in colors:  colors[b] = 'black'
         if b not in symbols:  symbols[b] = 'o'
   if kwargs.get('relative', False):
      if self.data[bands[0]].model is not None:
         rel_off = min(self.data[bands[0]].model)
      else:
         rel_off = min(self.data[bands[0]].mag)
   else:
      rel_off = 0

   # Should we plot the y-axis upside-down?
   if not flux:
      flip = 1
      ylabel = 'magnitude'
   else:
      flip = 0
      ylabel = 'relative flux'
   if not single:
      if cols is None and rows is None:
         cols = int(round(sqrt(n_plots)))
         rows = (n_plots // cols)
      elif cols is None:
         cols = n_plots // rows
      else:
         rows = n_plots // cols
      if n_plots % cols:  rows += 1
      figwidth = min(8*dpi,max_width)*1.0/dpi
      figheight = min(8.*rows//cols*dpi, max_height)*1.0/dpi
      p = myplotlib.PanelPlot(cols, rows, num=110, figsize=(figwidth,figheight))
   else:
      p = myplotlib.SimplePlot(num=110, figsize=(8,8))

   title = kwargs.get('title', self.name)
   p.title(title)
   if kwargs.get('epoch', False):
      p.xlabel('Days after B maximum')
   else:
      p.xlabel('JD - %d' % Toff)
   p.ylabel(ylabel)
   for ax in p.axes:
      if flip:  
         if not ax.yaxis_inverted():  ax.invert_yaxis()
      if xrange is not None:
         ax.set_autoscalex_on(False)
         ax.set_xlim((xrange[0]-Toff,xrange[1]-Toff))
      if yrange is not None:
         #print yrange
         ax.set_autoscaley_on(False)
         ax.set_ylim(yrange)
      if Nxticks is not None:
         ax.xaxis.set_major_locator(MaxNLocator(Nxticks, prune=prunex))
      else:
         ax.xaxis.set_major_locator(AutoLocator())

      if Nyticks is not None:
         ax.yaxis.set_major_locator(MaxNLocator(Nyticks, prune=pruney))
      else:
         ax.xaxis.set_major_locator(AutoLocator())

   offsets = self.lc_offsets()

   for i,filt in enumerate(bands):
      # Add extra space, if needed
      #delt = max(i*dm, i*dm + maxes[0] - maxes[i])
      #delt = round(delt, 1)
      delt = rel_off
      if offset and single:
         delt = delt + offsets[i]
      if not single:
         ax = p.axes[bidx[filt]]
         # make a reference to the lightcruve we are plotting.
         ax.lc = self.data[filt]
      else:
         ax = p.axes[0]
      plot_SN_panel(self, ax, filt, delt, symbols[filt], colors[filt], Toff,
            **kwargs)

   if oobj is not None:
      if single:
         raise RuntimeError("Ploting multple objects in a single panel is not supported")
      # PLotting a second object on top of this one.
      # See  what filters we're going to use:
      if oobj.filter_order is not None:
         obands = oobj.filter_order
      else:
         obands = list(oobj.data.keys())
         eff_wavs = []
         for filt in obands:
            eff_wavs.append(fset[filt].ave_wave)
         eff_wavs = asarray(eff_wavs)
         ids = argsort(eff_wavs)
         oobj.filter_order = [obands[i] for i in ids]
         obands = oobj.filter_order
      if epoch and getattr(oobj,'Tmax', None) is not None:
         oToff = oobj.Tmax
      else:
         oToff = Toff
      for filt in obands:
         # Figure out the pairing
         if filt in bands:
            # one-to-one
            #i = bands.index(filt)
            i = bidx[filt]
         elif filt in list(self.restbands.dict.values()):
            for band in self.restbands:
               if self.restbands[band] == filt and band in bands:
                  i = bidx[filt]
                  break
         else:
            i = -1

         if i < 0:  continue
         plot_SN_panel(oobj, p.axes[i], bands[i], rel_off, 's', 'r', oToff,
               **kwargs)
   
   #p.draw()
   p.set_limits(dox=(xrange is None), doy=(yrange is None), all_equal=1)
   p.draw()
   outfile = kwargs.get('outfile', None)
   if outfile is not None:
      p.fig.savefig(outfile)
   return(p)

def plot_lira(t, t2, t_maxes, BV, eBV, BV2, tmin, tmax, c):
   p = pyplot.figure(113)
   p.clear()
   ax = p.add_subplot(111)
   ax.set_xlabel('Epoch (Vmax)')
   ax.set_ylabel('B-V')

   ax.errorbar(t-t_maxes[0], BV, yerr=eBV, fmt='o', capsize=0, color='black')
   ax.plot(t2, BV2, 'o', mfc='red', label='B-V (deredshifted)')
   ax.plot([tmin, tmax], [0.732-0.0095*(tmin-55), 0.732-0.0095*(tmax-55)], '-',
         color='blue', label='Lira Law')
   ax.plot([tmin, tmax], [c[0]+c[1]*(tmin-55), c[0]+c[1]*(tmax-55)], '-', color='red',
         label='Fit')
   ax.legend(prop={'size':12})
   p.canvas.draw()
   return p

def plot_sBV(self, t, BV, eBV, tBVmax, etBVmax, ts, ys, eys, Tmax):
   p = pyplot.figure(115)
   p.clear()
   ax = p.add_subplot(111)
   ax.set_xlabel('t - t(Bmax) (days)')
   ax.set_ylabel('B-V')

   ax.errorbar(t-Tmax, BV, yerr=eBV, fmt='o', capsize=0, color='black')
   ax.plot(ts-Tmax, ys, '-', color='blue')
   ax.fill_between(ts-Tmax, ys-eys, ys+eys, facecolor='blue', alpha=0.2, 
         edgecolor='none')
   ax.axvline(tBVmax-Tmax, color='red')
   ax.axvspan(tBVmax-Tmax-etBVmax, tBVmax-Tmax+etBVmax, color='red', alpha=0.5)
   p.tight_layout()
   p.canvas.draw()
   return p


def plot_color(self, f1, f2, epoch=True, deredden=True, interp=False, 
      dokcorr=False, outfile=None, clear=True):
   '''Plot the color ([f1]-[f2]) evolution curve for the SN.  If  [epoch]
   is True and Bmax is defined, plot relative to T(Bmax).  If [deredden]
   is True, remove MW reddening.  Specify [outfile] to save to file.'''
   p = pyplot.figure(114)
   if clear: p.clear()
   ax = p.add_subplot(111)
   ax.set_xlabel('JD - JD(Bmax)')
   ax.set_ylabel('%s-%s' % (f1,f2))

   MJD,BV,eBV,flag = self.get_color(f1, f2, interp=interp, use_model=0, 
         dokcorr=dokcorr)
   if epoch:
      if self.Tmax:
         t0 = self.Tmax
      elif 'B' in self.data and getattr(self.data['B'],'Tmax',None) is not None:
         t0 = self.B.Tmax
      else:
         t0 = 0
   else:
      t0 = 0

   if deredden:
      Ia_w,Ia_f = kcorr.get_SED(0, 'H3')
      R1 = fset[f1].R(wave=Ia_w, flux=Ia_f, Rv=3.1)
      R2 = fset[f2].R(wave=Ia_w, flux=Ia_f, Rv=3.1)
      BV = BV - (R1-R2)*self.EBVgal
      eBV = sqrt(power(eBV,2) + (R1-R2)**2*0.1**2*self.EBVgal**2)

   bids = equal(flag, 0)
   rids = equal(flag, 1)
   ax.errorbar(MJD[bids]-t0, BV[bids], yerr=eBV[bids], fmt='o', capsize=0, 
         color='black', label='obs')
   if any(rids):
      ax.errorbar(MJD[rids]-t0, BV[rids], yerr=eBV[rids], fmt='o', capsize=0, 
         color='red', mfc='red', label='inter')
   ax.legend(prop={'size':12})
   pyplot.tight_layout()
   pyplot.draw()
   #p.canvas.draw()
   if outfile is not None:
      p.savefig(outfile)
   return p

def plot_lc(self, epoch=1, flux=0, symbol=4, outfile=None, use_model=True):
   # clear out any previous bindings...  gotta be a better way to do this...
   if pyplot.fignum_exists(111):
      for lc in list(self.parent.data.values()):
         mp = getattr(lc, 'mp', None)
         if mp is not None:
            if getattr(mp, 'bc',None) is not None:
               mp.bc.disconnect()
   if flux: 
      flipaxis = 0
   else:
      flipaxis = 1
   if self.interp is not None or self.band in self.parent.model._fbands:
      self.mp = myplotlib.PanelPlot(1,2,pheights=[0.2,0.8], num=111)
      self.mp.axes[0].set_ylabel('residuals')
      if flux:
         self.mp.axes[1].set_ylabel('flux (photons/s/cm$^2$)')
      else:
         self.mp.axes[1].set_ylabel('mag')
      p = self.mp.axes[1]
      p2 = self.mp.axes[0]
      # some references so that we can get at the data from within callbacks
      p.epoch = epoch
      p.flux = flux
      p.lc_inst = self
      p2.lc_inst = self
      if flipaxis:
         p.invert_yaxis()
         p2.invert_yaxis()
   else:
      self.mp = myplotlib.PanelPlot(1,1, num=111)
      if flux:
         self.mp.axes[0].set_ylabel('flux (photons/s/cm$^2$)')
      else:
         self.mp.axes[0].set_ylabel('mag')
      p = self.mp.axes[0]
      p.lc_inst = self
      p.epoch = epoch
      p.flux = flux
      if flipaxis:  p.invert_yaxis()
   if self.parent.Tmax is None:
      Tmax = 0
   else:
      Tmax = self.parent.Tmax

   self.mp.xlabel('Epoch (days)')
   self.mp.title('%s %s lightcurve\n' % (self.parent.name, self.band))
      
   if flux:
      y = self.flux
      ey = self.e_flux
   else:
      y = self.mag
      ey = self.e_mag
   # Plot the actual data (in the top panel)
   p.errorbar(self.MJD - epoch*Tmax, y, yerr=ey, barsabove=True, 
         capsize=0, elinewidth=1, fmt='o', mfc='blue', linestyle='None',
           ecolor='black', picker=True)
   if not all(self.mask):
      # Plot any masked out data s red filled symbols
      p._x, = p.plot(self.MJD[~self.mask] - epoch*Tmax, y[~self.mask], 'o', 
            color='red')
   else:
      p._x = None

   # Order of preference:  plot the model, else plot the spline, else
   #  plot nothing. If both and use_model is True, use the model
   if self.band in self.parent.model._fbands and use_model:
      t = arange(-10,70,1.0) + Tmax
      m,em,mask = self.parent.model(self.band, t)
      m_m,m_em,m_mask = self.parent.model(self.band, self.MJD)
      x = self.MJD
      if flux and self.parent.model.model_in_mags:
         p._model = p.plot(t[mask]-epoch*Tmax, power(10, -0.4*(m[mask]-self.filter.zp)), 
               '-', color='black')[0]
         y = self.flux - power(10,-0.4*(m_m - self.filter.zp))
         dy = self.e_flux
      elif not flux and not self.parent.model.model_in_mags:
         p._model = p.plot(t[mask]-epoch*Tmax, -2.5*log10(m[mask]) + \
               self.filter.zp, '-')[0]
         y = self.mag + 2.5*log10(m_m) + self.filter.zp
         dy = self.e_mag
      elif flux:
         p._model = p.plot(t[mask]-epoch*Tmax,m[mask], '-',
               color='black')[0]
         y = self.flux - m_m
         dy = self.e_flux
      else:
         p._model = p.plot(t[mask]-epoch*Tmax,m[mask],'-', color='black')[0]
         y = self.mag - m_m
         dy = self.e_mag

      # The residuals plot in the lower panel
      p2._errb = p2.errorbar(x - epoch*Tmax, y, yerr=dy, fmt='o', linestyle='None',
            mfc='blue', capsize=0, elinewidth=1, barsabove=True, ecolor='black',
            picker=True)
      if not all(self.mask):
         p2._x, = p2.plot(x[~self.mask] - epoch*Tmax, y[~self.mask], 'o', 
               color='red')
      else:
         p2._x = None
      p2.axhline(y=0)

   elif self.interp is not None:
      if self.tmin is not None and self.tmax is not None:
         t = arange(self.tmin, self.tmax+1, 1.0)
      else:
         t = arange(self.MJD[0], self.MJD[-1] + 1, 1.0)
      self._t = t
      self._epoch = epoch
      m,mask = self.eval(t, t_tol=-1)
      m_model,m_mask = self.eval(self.MJD, t_tol=-1)
      if flux:
         p._model = p.plot(compress(mask,t - epoch*Tmax), 
                compress(mask, power(10, -0.4*(m - self.filter.zp))),'-',
                color='black')[0]
         x = self.MJD
         y = self.flux - power(10, -0.4*(m_model-self.filter.zp))
         dy = self.e_flux
      else:
         p._model = p.plot(compress(mask,t - epoch*Tmax), compress(mask,m),'-',
               color='black')[0]
         x = self.MJD
         y = self.mag - m_model
         dy = self.e_mag
      # The residuals
      p2._errb = p2.errorbar(x - epoch*Tmax,y, yerr=dy, barsabove=True,
         capsize=0, elinewidth=1, fmt='o', mfc='blue', linestyle='None',
         ecolor='black', picker=True)
      if not all(self.mask):
         p2._x, = p2.plot(x[~self.mask] - epoch*Tmax, y[~self.mask], 'o',
               color='red')
      else:
         p2._x = None

   self.mp.set_limits()
   if len(self.mp.axes) > 1:
      self.mp.axes[0].axhline(0, color='black')
      if flux:
         self.mp.axes[0].set_ylim((y[self.mask*m_mask].min(),
            y[self.mask*m_mask].max()))
      else:
         self.mp.axes[0].set_ylim((y[self.mask*m_mask].max(),
            y[self.mask*m_mask].min()))
   self.mp.draw()
   self.mp.bc = ButtonClick(self.mp.fig)
   self.mp.bc.connect()
   if outfile is not None:
      self.mp.fig.savefig(outfile)
   return(self.mp)

def replot_lc(self):
   '''After mucking around, we want to update the fit, if needed.'''
   # First, transfer the current mask to the interp's mask
   if not pyplot.fignum_exists(111):  return
   fig = pyplot.figure(111)
   if fig is not self.mp.fig:
      return

   if self.interp is None or getattr(self, '_t', None) is None:
      for p in self.mp.axes:
         # Only need to deal with possible mask
         if p._x:  p._x.remove()
         if not all(self.mask):
            xx,yy = p.lines[0].get_data()
            p._x, = p.plot(xx[~self.mask], yy[~self.mask], 'o', color='red')
            self.mp.fig.canvas.draw()
      return

   self.interp.mask = self.mask
   self.interp.setup = False
   m,mask = self.eval(self._t, t_tol=-1)
   m_model,m_mask = self.eval(self.MJD, t_tol=-1)
   p = self.mp.axes[1]
   p2 = self.mp.axes[0]
   t = self._t
   epoch = self._epoch

   #clear out the previous stuff
   p._model.remove()
   p2._errb[0].remove()
   p2._errb[2][0].remove()
   if p._x:  p._x.remove()
   if p2._x:  p2._x.remove()

   if self.parent.Tmax is not None:
      Tmax = self.parent.Tmax
   else:
      Tmax = 0
   if p.flux:
      p._model = p.plot(compress(mask,t - epoch*Tmax), 
             compress(mask, power(10, -0.4*(m - self.filter.zp))),'-',
             color='black')[0]
      x = self.MJD
      y = self.flux - power(10, -0.4*(m_model-self.filter.zp))
      yy = self.flux
      dy = self.e_flux
   else:
      p._model = p.plot(compress(mask,t - epoch*Tmax), compress(mask,m),'-',
            color='black')[0]
      x = self.MJD
      y = self.mag - m_model
      yy = self.mag
      dy = self.e_mag
   p2._errb = p2.errorbar(x - epoch*Tmax,y, yerr=dy, barsabove=True,
      capsize=0, elinewidth=1, fmt='o', mfc='blue', linestyle='None',
      ecolor='black')
   if not all(self.mask):
      p._x, = p.plot(self.MJD[~self.mask] - epoch*Tmax, yy[~self.mask], 'o', 
            color='red')
      p2._x, = p2.plot(x[~self.mask] - epoch*Tmax, y[~self.mask], 'o',
            color='red')
   else:
      p._x = None
      p2._x = None
   if p.flux:
      p2.set_ylim(y[self.mask*m_mask].min(), 
            y[self.mask*m_mask].max())
   else:
      p2.set_ylim(y[self.mask*m_mask].max(), 
            y[self.mask*m_mask].min())
   self.mp.fig.canvas.draw()


def plot_kcorrs(self, colors=None, symbols=None, outfile=None):
   '''Plot the k-corrections, both mangled and un-mangled.'''
   # See  what filters we're going to use:
   bands = list(self.ks.keys())
   if len(bands) == 0:
      return

   if colors is None:  colors = color_dict()
   if symbols is None:  symbols = symbol_dict()

   eff_wavs = []
   for filt in bands:
      eff_wavs.append(fset[filt].ave_wave)
   eff_wavs = asarray(eff_wavs)
   ids = argsort(eff_wavs)
   bands = [bands[i] for i in ids]
      
   n_plots = len(bands)
   p = myplotlib.PanelPlot(1, n_plots, num=112, figsize=(6,n_plots))
        
   p.title("Use 'm' to plot mangled SED for any point")
   p.xlabel('Epoch (days)')
   p.ylabel('K-corrections')

   for i in range(len(bands)):
      #p.axes[i].yaxis.set_major_locator(MaxNLocator(5))
      b = bands[i]
      x = self.data[b].MJD - self.Tmax
      days = (arange(int(x[0]), int(x[-1]), 1)) #/(1+self.z)/self.ks_s
      rest_days = days/(1+self.z)/self.ks_s
      k,k_m = list(map(array, kcorr.kcorr(rest_days.tolist(), self.restbands[b], 
         b, self.z, self.EBVgal, 0.0, version=self.k_version)))
      k_m = equal(k_m, 1)
      p.axes[i].filt = b
      p.axes[i].inst = self
      p.axes[i].plot(days[k_m],k[k_m], '-', color=colors[b])
      p.axes[i].plot(x[self.ks_mask[b]], self.ks[b][self.ks_mask[b]], 
            symbols[b], color=colors[b])
      p.axes[i].text(0.9, 0.9, b, transform=p.axes[i].transAxes,
            verticalalignment='top', horizontalalignment='right')

   p.set_limits(all_equal=0)
   p.draw()
   p.bc = ButtonClick(p.fig, bindings={'m':plot_mangled_SED})
   p.bc.connect()
   if outfile is not None:
      p.fig.savefig(outfile)
   return(p)

def plot_mangled_SED(event):

   f = pyplot.figure(num=113)
   f.clear()
   ax = f.add_subplot(111)
   ax.set_xlabel('Wavelength (Angstroms)')
   ax.set_ylabel('Flux')
   band = event.inaxes.filt
   self = event.inaxes.inst

   x = event.xdata;  y = event.ydata
   id = argmin(power(x - self.data[band].t,2) + power(y - self.ks[band],2))
   day = int(self.data[band].t[id])
   wave,flux = kcorr.get_SED(day, version='H3')
   ax.plot(wave,flux, label='Original SED', color='black')
   if band in self.ks_mopts:
      if 'state' in self.ks_mopts[band][id]:
         man_flux = mangle_spectrum.apply_mangle(wave,flux, 
               **self.ks_mopts[band][id])[0]
      else:
         args = {}
         args['state'] = dict(ave_waves=self.ks_mopts[band][id]['sw'])
         args['pars'] = self.ks_mopts[band][id]['sf']
         man_flux = mangle_spectrum.apply_mangle(wave,flux,init=False,**args)[0]
   else:
      man_flux = flux
   ax.plot(wave, man_flux, label='Mangled SED', color='darkgreen')
   ax.plot(wave, man_flux/flux*ax.get_ylim()[1], color='orange')

   f1 = fset[band]
   f2 = fset[self.restbands[band]]
   ax.plot(f1.wave, f1.resp/f1.resp.max()*ax.get_ylim()[1], color='red')
   ax.plot(f2.wave*(1+self.z), f2.resp/f2.resp.max()*ax.get_ylim()[1], 
         color='blue')
   #wmax = max(f1.wave.max(), f2.wave.max()*(1+self.z))
   #wmin = max(f1.wave.min(), f2.wave.min()*(1+self.z))
   #dw = wmax - wmin
   #wmax += dw*0.05
   #wmin -= dw*0.05

   ax.legend()
   #ax.set_xlim(wmin,wmax)

   pyplot.tight_layout()
   f.canvas.draw()

def bind_c(event):
   if event.key == 'c':
      # re-compute the lc params.
      if not event.inaxes:  return
      event.inaxes.lc.compute_lc_params()
      draw_lc_params(event.inaxes.lc, event.inaxes)
      event.inaxes.figure.canvas.draw()

def draw_lc_params(self, ax):
   # draw lines that indicate the light-curve parameters
   for t in getattr(self, '_lc_labs', []):
      try:
         t.remove()
      except ValueError:
         pass
   self._lc_labs = []
   if self.Mmax and self.Tmax:
      t = ax.axhline(self.Mmax, color='0.5');  self._lc_labs.append(t)
      t = ax.text(self.Tmax, self.mag.max(), "$T_{max} = %.1f$" % self.Tmax,
         rotation=90, va='bottom', ha='right', color='0.5')
      self._lc_labs.append(t)
      t = ax.axvline(self.Tmax, color='0.5');  self._lc_labs.append(t)
      t = ax.text(self.Tmax+15*(1+self.parent.z), self.Mmax + 0.1, "$M_{max} = %.2f$" % self.Mmax,
         va='top', ha='left', color='0.5')
      self._lc_labs.append(t)
   if self.Mmax and self.dm15:
      t = ax.axhline(self.Mmax + self.dm15, color='0.5');  self._lc_labs.append(t)
      self._lc_labs.append(t)
      t = ax.text(self.Tmax+15*(1+self.parent.z), self.Mmax+self.dm15+0.1,
         "$\Delta m_{15} = %.2f$" % self.dm15, va='top', ha='right', 
         color='0.5')
      self._lc_labs.append(t)

def launch_int_fit(self, fitflux=False):
   '''Launch an interacive interpolator.'''
   from snpy.utils import InteractiveFit

   if not fitflux:
      ylabel = 'mag'
   else:
      ylabel = 'flux'

   mfit = InteractiveFit.InteractiveFit(self.interp, fignum=111, draw=False,
         xlabel='Epoch', ylabel=ylabel, extra_bindings={'q':lambda event: self.compute_lc_params()})
   if not fitflux:
      mfit.mp.axes[0].invert_yaxis()
      mfit.mp.axes[1].invert_yaxis()
   mfit.mp.axes[0].lc = self
   mfit.mp.axes[1].lc = self
   mfit.mp.fig.canvas.mpl_connect('key_press_event', bind_c)
   mfit.bind_help.append(('c','(re)compute light-curve parameters'))
   mfit.draw()
   return mfit

