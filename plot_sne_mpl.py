'''plot_sne_pg.py

   A collection of plotting routines that SNPY uses to display data.  This
   version uses the matplotlib library.
'''
from numpy import *
import matplotlib
#matplotlib.use('TkAgg')
import myplotlib
from matplotlib import pyplot,rcParams
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from snpy import kcorr
import types, string
from snpy.filters import fset
from snpy import mangle_spectrum

#rcParams['font.family'] = 'serif'
rcParams['font.size'] = 12

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
            raise KeyError, k

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
            raise KeyError, k

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

#class MaskData:
#   '''Given a figure instance, oversee the key bindings in a figure that
#   deal with masking and unmasking data.  This is done by manipulating the
#   lightcurve (lc) instance's mask array directly.  The default bindings are
#   'x' to both mask and unmask data.  These can be overridden by specifying
#   the maskkey and unmaskkey attributes.'''
#
#   def __init__(self, figure, maskkey='x', unmaskkey='x'):
#      self.figure = figure
#      self.mk = maskkey
#      self.umk = umaskkey


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
      if len(self.mouse_bindings.keys()) > 0:
         self.id2 = self.figure.canvas.mpl_connect('button_press_event', self.buttonpress)
      return self.id

   def disconnect(self):
      if self.id is None:  return
      self.figure.canvas.mpl_disconnect(self.id)

   def buttonpress(self, event):
      if event.inaxes is None:  return
      if self.pending_key is not None:  return
      if 'mouse'+str(event.button) in self.mouse_bindings:
         return apply(self.mouse_bindings['mouse'+str(event.button)], (event,))
      #elif event.button == 1:
      #   print "%f %f" % (event.xdata, event.ydata)


   def keypress(self, event):
      #print "key press event:", event.inaxes,event.key
      if event.inaxes is None:  return
      if event.key in self.bindings:
         return apply(self.bindings[event.key], (event,))

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
         lc.mask[id] = -lc.mask[id]
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



def plot_filters(self, bands=None, day=0, fill=0):
   '''Plot the filter responses, both the observed ones and the restbands
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

   # Next, for each filter, plot the response and rest band closest to it:
   if bands is None:  bands = self.data.keys()
   if type(bands) is types.StringType:
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

   pyplot.draw()
   p.cb = ButtonClick(p, bindings={'d':change_SED, 'D':change_SED})
   p.cb.connect()
   return p


def plot_sn(self, xrange=None, yrange=None, device=None, 
      title=None, interactive=0, single=0, dm=1, fsize=12., linewidth=1,
      symbols=None, colors=None, relative=0, legend=1, mask=1, label_bad=0,
      flux=0, epoch=1, msize=6, **pargs):
   '''Plot out the supernova data in a nice format.  There are several 
   options:
      - xrange,yrange:  specify the ranges to plot as lists [xmin,xmax], 
        [ymin,ymax]
      - title:  optional title
      - device:  which device to use.  Useful for output to file
      - interactive:  allows for an interactive plot.
      - single:  plot out as a single (rather than panelled) plot?
      - dm:  offset in magnitudes between the lightcurves (for single plots)
      - fsize:  override the font size used to plot the graphs
      - linewidth:  override the line width
      - symbols:  dictionary of symbols, indexed by band name.
      - colors:  dictionary of colors to use, indexed by band name.
      - relative:  plot only relative magnitudes (normalized to zero)?
      - legend:  do we plot the legend?
      - mask:  Omit plotting masked out data?
      - label_bad:  label the masked data with red x's?
   '''

   if not xrange:  xrange=self.xrange
   if not yrange:  yrange=self.yrange
   if not symbols:  symbols = symbol_dict()
   if not colors:  colors = color_dict()
   if fsize is None:  fsize = 12
   if linewidth is None:  linewidth=1

   # See  what filters we're going to use:
   if self.filter_order is not None:
      bands = self.filter_order
   else:
      bands = self.data.keys()
      eff_wavs = []
      for filter in bands:
         eff_wavs.append(fset[filter].ave_wave)
      eff_wavs = asarray(eff_wavs)
      ids = argsort(eff_wavs)
      self.filter_order = [bands[i] for i in ids]
      bands = self.filter_order
      
   for b in bands:
      if not single:
         colors[b] = 'blue'
         symbols[b] = 'o'
      else:
         if b not in colors:  colors[b] = 'black'
         if b not in symbols:  symbols[b] = 'o'
   if relative:
      if self.data[bands[0]].model is not None:
         rel_off = min(self.data[bands[0]].model)
      else:
         rel_off = min(self.data[bands[0]].mag)
   else:
      rel_off = 0
   
   n_plots = len(bands)
   # SHould we plot the y-axis upside-down?
   if not flux:
      flip = 1
      ylabel = 'magnitude'
   else:
      flip = 0
      ylabel = 'relative flux'
   if not single:
      cols = int(round(sqrt(n_plots)))
      rows = (n_plots / cols)
      if n_plots % cols:  rows += 1
      p = myplotlib.PanelPlot(cols, rows, num=110, figsize=(8,8.*rows/cols))
   else:
      p = myplotlib.SimplePlot(num=110, figsize=(8,8))

   if title is not None:  
      p.title(title)
   else:
      p.title(self.name)
   if epoch:
      p.xlabel('Days after B maximum')
   else:
      p.xlabel('Date')
   p.ylabel(ylabel)
   for ax in p.axes:
      if flip:  
         if not ax.yaxis_inverted():  ax.invert_yaxis()
      if xrange is not None:
         ax.set_autoscalex_on(False)
         ax.set_xlim(xrange)
      if yrange is not None:
         #print yrange
         ax.set_autoscaley_on(False)
         ax.set_ylim(yrange)

   i = 0
   offsets = self.lc_offsets()
   for filter in bands:
      # Add extra space, if needed
      #delt = max(i*dm, i*dm + maxes[0] - maxes[i])
      #delt = round(delt, 1)
      delt = offsets[i]
      if not single:
         ax = p.axes[i]
         # make a reference to the lightcruve we are plotting.
         ax.lc = self.data[filter]
      else:
         ax = p.axes[0]
      label = filter
      ax.mylabels = []
      if not single:
         ax.mylabels.append(ax.text(0.9, 0.9, label, transform=ax.transAxes, 
            horizontalalignment='right', fontsize=fsize, 
            verticalalignment='top'))
      if mask:
         x = self.data[filter].MJD[self.data[filter].mask]
         if not flux:
            y = self.data[filter].mag[self.data[filter].mask] + \
                  single*delt - relative*rel_off
            ey = self.data[filter].e_mag[self.data[filter].mask]
         else:
            y = self.data[filter].flux[self.data[filter].mask]* \
                  power(10, -0.4*(single*delt-relative*rel_off))
            ey = self.data[filter].e_flux[self.data[filter].mask]
      else:
         x = self.data[filter].MJD
         if not flux:
            y = self.data[filter].mag + single*delt - relative*rel_off
            ey = self.data[filter].e_mag
         else:
            y = self.data[filter].flux* \
                  power(10, -0.4*(single*delt-relative*rel_off))
            ey = self.data[filter].e_flux

      if self.Tmax is None:  
         Tmax = 0
      else:
         Tmax = self.Tmax
      ax.errorbar(x-Tmax*epoch, y, yerr=ey, barsabove=True, capsize=0,
            elinewidth=1, fmt=symbols[filter], ms=msize, 
            mfc=colors[filter], label=filter+'+'+'%.1f' % delt, linestyle='None',
              ecolor='black')
      if label_bad:
         gids = equal(self.data[filter].mask, 0)
         if sometrue(gids):
            x = self.data[filter].MJD[gids] - Tmax*epoch
            if not flux:
               y = self.data[filter].mag[gids] + \
                     single*delt - relative*rel_off
            else:
               y = self.data[filter].flux[gids]* \
                     power(10, -0.4*(single*delt-relative*rel_off))
            ax.plot(x, y, marker='x', mec='red', ms=12, mew=1, linestyle='')

      # Now check to see if there is a model to plot:
      if self.model.Tmax is not None and filter in self.model._fbands:
         t = arange(-10,70,1.0) + self.Tmax
         mag,err,gids = self.model(filter, t)
         if not flux:
            y = mag + single*delt - relative*rel_off
         else:
            zp = fset[filter].zp
            y = power(10, -0.4*(mag - zp + single*delt - \
                                relative*rel_off))
         ax.plot(compress(gids,t-self.Tmax*epoch), compress(gids,y), 
               color='k', linewidth=linewidth)
         l = ax.plot(compress(gids,t-self.Tmax*epoch), compress(gids,y+err), 
               '--',color='k', linewidth=linewidth)
         l[0].autoscale=False
         l = ax.plot(compress(gids,t-self.Tmax*epoch), compress(gids,y-err), 
               '--',color='k', linewidth=linewidth)
         l[0].autoscale=False
      elif self.data[filter].interp is not None:
         d = self.data[filter].interp.domain()
         t = arange(d[0], d[1]+1.0, 1.0)
         mag,gids = self.data[filter].eval(t, t_tol=-1)
         if not flux:
            y = mag + single*delt - relative*rel_off
         else:
            zp = fset[filter].zp
            y = power(10, -0.4*(mag - zp + single*delt - \
                      relative*rel_off))
         ax.plot(t[gids] - self.Tmax*epoch, y[gids], color='k',
               linewidth=linewidth)
      i = i + 1

      if single and legend:
         ax.legend(loc='lower right', numpoints=1, ncol=1, prop={'size':'small'})
         ax.text(0.05,0.05,self.name,transform=ax.transAxes)
   
   #p.draw()
   p.set_limits(dox=(xrange is None), doy=(yrange is None), all_equal=1)
   p.draw()
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

def plot_lc(self, device='/XSERVE', epoch=1, flux=0, symbol=4):
   # clear out any previous bindings...  gotta be a better way to do this...
   if flux: 
      flipaxis = 0
   else:
      flipaxis = 1
   if self.interp is not None or self.band in self.parent.model._fbands:
      self.mp = myplotlib.PanelPlot(1,2,pheights=[0.2,0.8], num=111)
      self.mp.axes[0].set_ylabel('residuals')
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
   if not alltrue(self.mask):
      # Plot any masked out data s red filled symbols
      p._x, = p.plot(self.MJD[-self.mask] - epoch*Tmax, y[-self.mask], 'o', 
            color='red')
   else:
      p._x = None

   # Order of preference:  plot the model, else plot the spline, else
   #  plot nothing
   if self.band in self.parent.model._fbands:
      t = arange(-10,70,1.0) + Tmax
      m,em,mask = self.parent.model(self.band, t)
      m_m,m_em,m_mask = self.parent.model(self.band, self.MJD)
      x = self.MJD
      if flux and self.parent.model.model_in_mags:
         p._model = p.plot(t[mask]-epoch*Tmax, power(10, -0.4*(m-self.filter.zp)), '-',
               color='black')[0]
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
      if not alltrue(self.mask):
         p2._x, = p2.plot(x[-self.mask] - epoch*Tmax, y[-self.mask], 'o', 
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
      if not alltrue(self.mask):
         p2._x, = p2.plot(x[-self.mask] - epoch*Tmax, y[-self.mask], 'o',
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
   return(self.mp)

def replot_lc(self):
   '''After mucking around, we want to update the fit, if needed.'''
   # First, transfer the current mask to the interp's mask
   if not pyplot.fignum_exists(111):  return
   fig = pyplot.figure(111)
   if fig is not self.mp.fig:
      return

   if self.interp is None:
      p = self.mp.axes[0]
      # Only need to deal with possible mask
      if p._x:  p._x.remove()
      if not alltrue(self.mask):
         xx,yy = p.lines[0].get_data()
         p._x, = p.plot(xx[-self.mask], yy[-self.mask], 'o', color='red')
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
   if not alltrue(self.mask):
      p._x, = p.plot(self.MJD[-self.mask] - epoch*Tmax, yy[-self.mask], 'o', 
            color='red')
      p2._x, = p2.plot(x[-self.mask] - epoch*Tmax, y[-self.mask], 'o',
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


def plot_kcorrs(self, device='13/XW', colors=None, symbols=None):
   '''Plot the k-corrections, both mangled and un-mangled.'''
   # See  what filters we're going to use:
   bands = self.ks.keys()
   if len(bands) == 0:
      return

   if colors is None:  colors = color_dict()
   if symbols is None:  symbols = symbol_dict()

   eff_wavs = []
   for filter in bands:
      eff_wavs.append(fset[filter].ave_wave)
   eff_wavs = asarray(eff_wavs)
   ids = argsort(eff_wavs)
   bands = [bands[i] for i in ids]
      
   for b in bands:
      if b not in colors:  colors[b] = 1
      if b not in symbols:  symbols[b] = 4
   n_plots = len(bands)
   cols = int(round(sqrt(n_plots)))
   rows = (n_plots / cols)
   if n_plots % cols:  rows += 1
   p = myplotlib.PanelPlot(1, n_plots, num=112, figsize=(6,n_plots),
        nymax=5, prunex=None)
   p.title("Use 'm' to plot mangled SED for any point")
   p.xlabel('Epoch (days)')
   p.ylabel('K-corrections')

   for i in range(len(bands)):
      #p.axes[i].yaxis.set_major_locator(MaxNLocator(5))
      b = bands[i]
      x = self.data[b].MJD - self.Tmax
      days = (arange(int(x[0]), int(x[-1]), 1)) #/(1+self.z)/self.ks_s
      rest_days = days/(1+self.z)/self.ks_s
      k,k_m = map(array, kcorr.kcorr(rest_days.tolist(), self.restbands[b], 
         b, self.z, self.EBVgal, 0.0, version=self.k_version))
      k_m = equal(k_m, 1)
      p.axes[i].filter = b
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
   return(p)

def plot_mangled_SED(event):

   f = pyplot.figure(num=113)
   f.clear()
   ax = f.add_subplot(111)
   ax.set_xlabel('Wavelength (Angstroms)')
   ax.set_ylabel('Flux')
   band = event.inaxes.filter
   self = event.inaxes.inst

   x = event.xdata;  y = event.ydata
   id = argmin(power(x - self.data[band].t,2) + power(y - self.ks[band],2))
   day = int(self.data[band].t[id])
   wave,flux = kcorr.get_SED(day, version='H3')
   ax.plot(wave,flux, label='Original SED', color='black')
   if band in self.ks_mopts:
      man_flux = mangle_spectrum.apply_mangle(wave,flux, **self.ks_mopts[band][id])
   ax.plot(wave, man_flux, label='Mangled SED', color='darkgreen')

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

   f.canvas.draw()
