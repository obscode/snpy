'''plot_sne_pg.py

   A collection of plotting routines that SNPY uses to display data.  This
   version uses the matplotlib library.
'''
from numpy import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot,rcParams
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import myplotlib
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


class ButtonClick:
   '''Given a figure instance and an optional dictionary of key-bindings,
   oversee the key bindings in a figure.  There are a set of default
   bindings:  x:  select xrange,  y:  select yrange,  b:  select box
   range.  These can be overridden using the bindings keyword.  Simply
   set it to a dictionary with the button as the key and the call-back
   function as the value.'''

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
      if len(self.mouse_bindings.keys()) > 1:
         self.id2 = self.figure.canvas.mpl_connect('button_press_event', self.buttonpress)

   def disconnect(self):
      if self.id is None:  return
      self.figure.canvas.mpl_disconnect(self.id)

   def buttonpress(self, event):
      if event.inaxes is None:  return
      if self.pending_key is not None:  return
      if 'mouse'+str(event.button) in self.mouse_bindings:
         return apply(self.bindings['mouse'+str(event.button)], (event,))
      #elif event.button == 1:
      #   print "%f %f" % (event.xdata, event.ydata)


   def keypress(self, event):
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
      if event.key == 'q':
         self.disconnect()
         return

      self.figure.canvas.draw()

def change_SED(event):
   fig = event.inaxes.figure
   if event.key == 'd':
      if fig.current_day > -19:  fig.current_day += 1
   elif event.key == 'D':
      if fig.current_day < 70:  fig.current_day -= 1
   wave,flux = kcorr.get_SED(fig.current_day, version='H3')
   fig.sed_line.set_data(wave*(1+fig.z), flux/flux.max())
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
   cb = ButtonClick(p, bindings={'d':change_SED, 'D':change_SED})
   cb.connect()
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
         print yrange
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
         ax.plot(compress(gids,t-self.Tmax*epoch), compress(gids,y+err), 
               '--',color='k', linewidth=linewidth)
         ax.plot(compress(gids,t-self.Tmax*epoch), compress(gids,y-err), 
               '--',color='k', linewidth=linewidth)
      elif self.data[filter].tck is not None:
         tck = self.data[filter].tck
         t = arange(tck[0][0], tck[0][-1], 1.0)
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

def mask_data(self):
   p = self.plot(single=0, label_bad=1, epoch=1, mask=0)
   for b in self.data:
      # a place to store red X's
      self.data[b].pts = {}
   bc = ButtonClick(p, bindings={'mouse1':bind_mask_data, 'u':bind_unmask_data})
   bc.connect()

def bind_mask_data(event):
   '''Mask the data that has been selected.'''
   if event.inaxes is None:  return
   lc = event.inaxes.lc

   if not sometrue(lc.mask):
      return
   # find distances to valid (unmasked) data:
   dists = power(lc.mag[lc.mask] - event.ydata, 2) + \
           power(lc.t[lc.mask] - event.xdata, 2)
   i = argmin(dists)
   # need to find index in original (masked+unmasked) data
   id = indices(lc.mask.shape)[0][lc.mask][i]
   x = lc.t[id]
   y = lc.mag[id]
   lc.mask[id] = 0
            
   lc.pts[id] = event.inaxes.plot(x, y, marker='x', mec='red', ms=12, mew=1, 
         linestyle='None')[0]

def bind_unmask_data(event):
   '''Unmask the data that has been masked.'''
   if event.inaxes is None:  return
   lc = even.inaxes.lc

   if alltrue(lc.mask):
      return
   gids = logical_not(lc.mask)
   dists = power(lc.mag[gids] - event.ydata, 2) + \
           power(lc.t[gids] - event.xdata, 2)
   i = argmin(dists)
   id = indices(lc.mask.shape)[0][gids][i]
   lc.mask[id] = 1
   lc.pts[id].remove()
   del lc.pts[id]

def plot_lira(t, t2, t_maxes, BV, eBV, BV2, tmin, tmax, c):
   p = pyplot.figure(113)
   p.clear()
   ax = p.add_subplot(111)
   ax.xlabel('Epoch (Vmax)')
   ax.ylabel('B-V')

   ax.plot(t-t_maxes[0], BV, 'o', color='black')
   ax.errorbar(t-t_maxes[0], BV, yerr=eBV, capsize=0, color='black')
   ax.plot(t2, BV2, 'o', markercolor='red', label='B-V (deredshifted)')
   ax.plot([tmin, tmax], [0.732-0.0095*(tmin-55), 0.732-0.0095*(tmax-55)], '-',
         color='blue', label='Lira Law')
   ax.line([tmin, tmax], [c[0]+c[1]*(tmin-55), c[0]+c[1]*(tmax-55)], '-', color='red',
         label='Fit')
   ax.legend(loc='upper right')
   p.canvas.draw()
   return p

def lc_data_conv(p, x, y):
   '''Given a plot instance p and (x,y) coordinates, convert to the real
   data (which can be different because of plotting as magnitudes, flux,
   epoch or no, etc.'''
   # first, if plotting in epochs:
   realx = x + p.epoch*p.lc_inst.parent.Tmax
   if not p.flux:
      realy = power(10, -0.4*(y - p.lc_inst.filter.zp))
   return(realx,realy)


def add_knot(event):
   '''Add a knot point to the spline.  This is called through the interactive interface).'''
   if event.inaxes is None: return
   ax = event.inaxes
   self = ax.lc_inst
   if self.model_type is None or self.model_type != 'spline':
      print "Adding knots not supported for this model"
      return
   x,y = lc_data_conv(ax, event.xdata, event.ydata)
   if x <= self.tck[0][0] or x >= self.tck[0][-1]:
      print "Error, can't add knots outside the current range."
   old_knots = self.tck[0]*1.0
   knots = concatenate([[x],self.tck[0]])
   knots = sort(knots)
   self.tck = (knots, self.tck[1], self.tck[2])
   try:
      update_plot(self, task=-1)
   except:
      self.tck = (old_knots, self.tck[1], self.tck[2])
      update_plot(self, task=-1)

def move_knot(event):
   '''Move the point closest to x to a new location.'''
   if event.inaxes is None:  return
   ax = event.inaxes
   self = ax.lc_inst
   if self.model_type is None or self.model_type != 'spline':
      print "Moving knots not supported for this model"
      return
   xx,yy = ax._spl.get_data()
   if ax.pending_move is None:
      deltas = absolute(xx - event.xdata)
      id = argmin(deltas)
      ax.pending_move = click_line(ax, xx[id], yy[id])
      ax.pending_move.connect()
   else:
      x0,x1,y0,y1 = ax.pending_move.close()
      ax.pending_move = None
      deltas = absolute(xx - x0)
      id = argmin(deltas)
      k = self.tck[2]
      if id <= k or id >= len(deltas) - k -1 :
         print "You can only move interior points"
         return
      x,y = lc_data_conv(ax, x1, y1)
      old_knots = self.tck[0]*1.0
      self.tck[0][id] = x
      try:
         update_plot(self, task=-1)
      except:
         self.tck = (old_knots, self.tck[1], self.tck[2])
         update_plot(self, task=-1)

def delete_knot(event):
   if event.inaxes is None:  return
   ax = event.inaxes
   self = ax.lc_inst
   if self.model_type is None or self.model_type != 'spline':
      print "Deleting knots not supported for this model"
      return
   old_knots = self.tck[0]*1.0
   deltas = absolute(ax._spl.get_data()[0] - event.xdata)
   id = argmin(deltas)
   k = self.tck[2]
   print id
   if id <= k or id >= len(deltas) - k -1 :
      print "You can only delete interior points"
      return
   list = self.tck[0].tolist()
   del list[id]
   self.tck = (array(list), self.tck[1], self.tck[2])
   try:
      update_plot(self, task=-1)
   except:
      self.tck = (old_knots, self.tck[1], self.tck[2])
      update_plot(self, task=-1)

def change_s(event):
   if event.inaxes is None:  return
   ax = event.inaxes
   self = ax.lc_inst
   if self.model_type == 'spline2':
      print "Changing s not supported for this model"
      return
   if self.tck is not None:
      if event.key == '-':
         self.s -= sqrt(2*self.s)
         if self.s < 0:  self.s = 0
      elif event.key == '+':
         self.s += sqrt(2*self.s)
      else:
         self.s = len(self.MJD)
      update_plot(self, task=0)
   else:
      if self.line.N is not None:  self.line.N = None
      if self.line.sigma is None:  self.line.sigma=3.0
      if event.key == '+':
         self.line.sigma *= 1.1
      else:
         self.line.sigma /= 1.1
      update_plot(self, task=-2)

def change_N(event):
   if event.inaxes is None:  return
   ax = event.inaxes
   self = ax.lc_inst
   if self.line.sigma is not None:  self.line.sigma = None
   if self.line.N is None:  self.line.N = 5
   if event.key == 'n':
      if self.line.N > 1:  self.line.N -= 1
   else:
      self.line.N += 1
   update_plot(self, task=-2)

def update_plot(self, task):
   '''Update the plot with new spline info.'''
   inst = self.mp
   if len(inst.axes) != 2:  return

   k = self.tck[2]
   if task == 0 or task == 1:
      self.spline_fit(task=task, k=k, fitflux=self.model_flux,
            s=self.s, tmin=self.tmin, tmax=self.tmax, method='spline')
   elif task == -1:
      t = self.tck[0][k+1:-(k+1)]
      self.spline_fit(task=task, knots=t, k=k, fitflux=self.model_flux,
            s=self.s, tmin=self.tmin, tmax=self.tmax, method='spline')
   m,mask = self.eval(self._t, t_tol=None)
   m_model,m_mask = self.eval(self.MJD, t_tol=None)

   if task == -2:
      if self.line.N is not None:  N = self.line.N
      else: N = -1
      if self.line.sigma is not None:  sigma = self.line.sigma
      else: sigma = -1
      inst._title.set_text(string.split(inst._title.get_text(), '\n')[0]+\
            '\n' + 'N = %d   sigma = %.1f   r-chi-sq = %.2f' % \
            (N,sigma,self.line.chisq()))
   else:
      inst._title.set_text(string.split(inst._title.get_text(), '\n')[0]+\
            '\n' + 'Np = %d  Nk = %d  s = %.1f r-chi-sq = %.2f' % \
            (len(self.mag), len(self.tck[0]), self.s, self.rchisq))

   if not inst.axes[1].yaxis_inverted():
      y = compress(self.mask, self.flux - power(10, -0.4*(m_model-self.filter.zp)))
      inst.axes[1]._model.set_ydata(power(10, -0.4*(m - self.filter.zp)))
      set_errorbar_ydata(inst.axes[0]._errb, y)
      if task != -2:
         inst.axes[1]._spl.remove()
         inst.axes[1]._spl = inst.axes[1].plot(self.tck[0] - self._epoch*self.parent.Tmax, 
               power(10, -0.4*(self.eval(self.tck[0], t_tol=None)[0] - \
               self.filter.zp)), 'o', mfc='red')[0]
   else:
      y = compress(self.mask, self.mag - m_model)
      inst.axes[1]._model.set_ydata(m)
      set_errorbar_ydata(inst.axes[0]._errb, y)
      if task != -2:
         inst.axes[1]._spl.remove()
         inst.axes[1]._spl = inst.axes[1].plot(self.tck[0] - self._epoch*self.parent.Tmax, 
               self.eval(self.tck[0], t_tol=None)[0], 'o', mfc='red')[0]
   inst.draw()

def set_errorbar_ydata(err, ydata):
   '''Given the compound obect err representing an errorbar return value,
   set the ydata of all the components.'''
   line = err[0]
   old_ydata = line.get_ydata()
   line.set_ydata(ydata)
   dy = ydata - old_ydata
   for line in err[1]:
      line.set_ydata(line.get_ydata() + dy)
   for col in err[2]:
      paths = col.get_paths()
      for i in range(len(paths)):
         paths[i].vertices[:,1] += dy[i]

def plot_lc(self, device='/XSERVE', interactive=0, epoch=1, flux=0, gloes=0,
          symbol=4):
   # clear out any previous bindings...  gotta be a better way to do this...
   if flux: 
      flipaxis = 0
   else:
      flipaxis = 1
   if gloes or self.model_type is not None or self.band in self.parent.model._fbands:
      self.mp = myplotlib.PanelPlot(1,2,pheights=[0.2,0.8], num=111)
      self.mp.axes[0].set_ylabel('residuals')
      self.mp.axes[1].set_ylabel('mag')
      p = self.mp.axes[1]
      p.pending_move = None
      p2 = self.mp.axes[0]
      # some references so that we can get at the data from within callbacks
      p.lc_inst = self
      p2.lc_inst = self
      p.epoch = epoch
      p.flux = flux
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
      y = power(10, -0.4*(self.mag - self.filter.zp))
      ey = y*self.e_mag/1.0857
   else:
      y = self.mag
      ey = self.e_mag
   p.errorbar(self.MJD - epoch*Tmax, y, yerr=ey, barsabove=True, 
         capsize=0, elinewidth=1, fmt='o', mfc='blue', linestyle='None',
           ecolor='black')

   # Order of preference:  plot the model, else plot the spline, else
   #  plot nothing, unless we explicitly ask for GLoEss.
   if self.band in self.parent.model._fbands:
      t = arange(-10,70,1.0) + Tmax
      m,em,mask = self.parent.model(self.band, t)
      m_m,m_em,m_mask = self.parent.model(self.band, self.MJD)
      x = self.MJD[self.mask*m_mask]
      if flux and self.parent.model.model_in_mags:
         p._model = p.plot(t[mask]-epoch*Tmax, power(10, -0.4*(m-self.filter.zp)), '-')[0]
         y = compress(self.mask*m_mask, 
               self.flux - power(10,-0.4*(m_m - self.filter.zp)))
         dy = self.e_flux[self.mask*m_mask]
      elif not flux and not self.parent.model.model_in_mags:
         p._model = p.plot(t[mask]-epoch*Tmax, -2.5*log10(m[mask]) + self.filter.zp, '-')[0]
         y = compress(self.mask*m_mask, self.mag + 2.5*log10(m_m) + self.filter.zp)
         dy = self.e_mag[self.mask*m_mask]
      elif flux:
         p._model = p.plot(t[mask]-epoch*Tmax,m[mask], '-')[0]
         y = self.flux[self.mask*m_mask] - m_m[self.mask*m_mask]
         dy = self.e_flux[self.mask*m_mask]
      else:
         p._model = p.plot(t[mask]-epoch*Tmax,m[mask],'-')[0]
         y = self.mag[self.mask*m_mask] - m_m[self.mask*m_mask]
         dy = self.e_mag[self.mask*m_mask]

      p2._errb = p2.errorbar(x - epoch*Tmax, y, yerr=dy, fmt='o', linestyle='None',
            mfc='blue', capsize=0, elinewidth=1, barsabove=True, ecolor='black')
      p2.axhline(y=0)

   elif gloes or self.model_type is not None:
      if self.tmin is not None and self.tmax is not None:
         t = arange(self.tmin, self.tmax+1, 1.0)
      else:
         t = arange(self.MJD[0], self.MJD[-1] + 1, 1.0)
      self._t = t
      self._epoch = epoch
      m,mask = self.eval(t, t_tol=None)
      m_model,m_mask = self.eval(self.MJD, t_tol=None)
      if flux:
         p._model = p.plot(compress(mask,t - epoch*Tmax), 
                compress(mask, power(10, -0.4*(m - self.filter.zp))),'-')[0]
         if self.tck is not None:
            p._spl = p.plot(self.tck[0] - epoch*Tmax, power(10, -0.4*(self.eval(self.tck[0], t_tol=None)[0] \
               - self.filter.zp)), 'o', mfc='red')[0]
         x = compress(self.mask*m_mask, self.MJD)
         y = compress(self.mask*m_mask, self.flux - power(10, -0.4*(m_model-self.filter.zp)))
         dy = compress(self.mask*m_mask, self.e_flux)
         p2._errb = p2.errorbar(x - epoch*Tmax,y, yerr=dy, barsabove=True,
               capsize=0, elinewidth=1, fmt='o', mfc='blue', linestyle='None',
               ecolor='black')
      else:
         p._model = p.plot(compress(mask,t - epoch*Tmax), compress(mask,m),'-')[0]
         if self.tck is not None:
            p._spl = p.plot(self.tck[0] - epoch*Tmax, 
                  self.eval(self.tck[0], t_tol=None)[0], 'o', mfc='red')[0]
         x = compress(self.mask*m_mask, self.MJD)
         y = compress(self.mask*m_mask, self.mag - m_model)
         dy = compress(self.mask*m_mask, self.e_mag)
         p2._errb = p2.errorbar(x - epoch*Tmax,y, yerr=dy, barsabove=True,
               capsize=0, elinewidth=1, fmt='o', mfc='blue', linestyle='None',
               ecolor='black')
      # useful stats:
      if self.tck is not None:
         extra_title = 'Np = %d  Nk = %d  s = %.1f r-chi-sq = %.2f' % \
            (len(self.mag), len(self.tck[0]), self.s, self.rchisq) 
      else:
         if self.line.N is not None:  N = self.line.N
         else: N = -1
         if self.line.sigma is not None:  sigma = self.line.sigma
         else: sigma = -1
         extra_title = 'N = %d   sigma = %.1f   r-chi-sq = %.2f' % \
               (N,sigma,self.line.chisq())
      self.mp._title.set_text(string.split(self.mp._title.get_text(), '\n')[0]+\
                              '\n' + extra_title)
   self.mp.set_limits()
   if len(self.mp.axes) > 1:
      self.mp.axes[0].axhline(0)
   self.mp.draw()
   self.mp.bc = ButtonClick(self.mp.fig, bindings={'a':add_knot, 'd':delete_knot, 
      'm':move_knot, '-':change_s, '+':change_s, '=':change_s,
            'n':change_N, 'N':change_N})
   self.mp.bc.connect()
   return(self.mp)
   
   
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
   bc = ButtonClick(p.fig, bindings={'mouse1':plot_mangled_SED})
   bc.connect()
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
   ax.set_xlim(wmin,wmax)

   f.canvas.draw()
