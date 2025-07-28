'''plot_sne_pg.py

   A collection of plotting routines that SNPY uses to display data.  This
   version uses the PGPLOT library.
'''
from numpy import *
import pygplot
from snpy import kcorr
import types
from snpy.filters import fset
import mangle_spectrum

# Some default colors and labels to use in plotting.
default_colors = {'Bs':'blue', 'Vs':'green', 'Rs':'yellow', 'Is':'red',
                  'B':'blue', 'V':'green',
                  'oBs':'blue', 'oVs':'green',
                  'R4m':'yellow','I4m':'red',
                  'Yc':'orange', 'Jc':'purple', 'g_m':'blue', 'r_m':'green',
                  'Y':'cyan', 'J':'green', 'H':'red', 'K':'purple',
                  'i_m':'yellow', 'z_m':'red', 'u_s':'purple','g_s':'green', 'r_s':'yellow',
                  'z_s':'red',
                  'i':'yellow', 'u':'purple', 'g':'green', 'r':'red',
                  'oi':'yellow', 'ou':'purple', 'og':'green', 'or':'red',
                  'i_s':'red'}
default_symbols = {'Bs':4, 'Vs':6, 'Rs':7, 'Is':11, 'B':4, 'V':6,
                  'i':7, 'u':11, 'u_s':11, 'g':4, 'r':6,
                  'oi':7, 'ou':11, 'og':4, 'or':6, 'oBs':4, 'oVs':6,
                   'R4m':6, 'I4m':7,
                   'Yc':4, 'Jc':6, 
                   'Y':4, 'J':6, 'H':7, 'K':11, 
                   'g_m':4, 'r_m':6, 'i_m':7, 'z_m':11, 'g_s':4, 'r_s':6, 'i_s':7,'z_s':11}

def plot_filters(self, bands=None, day=0):
   '''Plot the filter responses, both the observed ones and the restbands
   defined in the sn instance in the rest frame of the SN.  If bands is specified,
   only plot those.  Specify which SED to plot using day.'''
   p = pygplot.Plot(xlabel='Wavelength (\A)', ylabel='relative flux', 
        title='Filter responses + %s SED (day %d)' % (self.k_version,day),
        device="11/XWIN")

   # first, get the SED:
   sed_wav,sed_flux = kcorr.get_SED(day, self.k_version)
   # normalize:
   sed_flux = sed_flux/maximum.reduce(sed_flux)
   p.sed_line = p.line(sed_wav*(1+self.z),sed_flux, autorange=0)
   p.current_day = day
   p.version = self.k_version
   p.z = self.z

   # Next, for each filter, plot the response and rest band closest to it:
   if bands is None:  bands = list(self.data.keys())
   if type(bands) is bytes:
      bands = [bands]

   for band in bands:
      f1 = fset[band]
      f2 = fset[self.restbands[band]]
      maxresp = maximum.reduce(f1.resp)
      max_wave = f1.wave[argmax(f1.resp)]
      p.line(f1.wave, f1.resp/maxresp, color='blue')
      p.label(max_wave,1.05, band, color='blue')
      maxresp = maximum.reduce(f2.resp)
      max_wave = f2.wave[argmax(f2.resp)]*(1+self.z)
      p.line(f2.wave*(1+self.z), f2.resp/maxresp, color='red')
      p.label(max_wave, 1.05, self.restbands[band], color='red')

   p.plot()
   p.interact(bindings={'D':advance_day, 'd':advance_day})
   pygplot.ppgplot.pgask(0)
   p.close()

def advance_day(x, y, key, p):
   '''A call-back function used by plot_filters() to chnage the eqpoch of the
   SNIa SED plotted in the background.'''

   if key == 'd':
      if p.current_day < 70:
         p.current_day += 1
         p.title = 'Filter responses + %s SED (day %d)' % (p.version,
               p.current_day)
         p.sed_line.x,p.sed_line.y = kcorr.get_SED(p.current_day, p.version)
         p.sed_line.y = p.sed_line.y/maximum.reduce(p.sed_line.y)
         p.sed_line.x = p.sed_line.x*(1+p.z)
   if key == 'D':
      if p.current_day > -19:
         p.current_day -= 1
         p.title = 'Filter responses + %s SED (day %d)' % (p.version,
                    p.current_day)
         p.sed_line.x,p.sed_line.y = kcorr.get_SED(p.current_day, p.version)
         p.sed_line.y = p.sed_line.y/maximum.reduce(p.sed_line.y)
         p.sed_line.x = p.sed_line.x*(1+p.z)


def plot_sn(self, xrange=None, yrange=None, device=None, 
      title=None, interactive=0, single=0, dm=1, fsize=1.0, linewidth=3,
      symbols=None, colors=None, relative=0, legend=1, mask=1, label_bad=0,
      flux=0, epoch=0, **pargs):
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
   if not symbols:  symbols = default_symbols
   if not colors:  colors = default_colors
   if device is None:  device = self.device
   if fsize is None:  fsize = 1.0
   if linewidth is None:  linewidth = 3

   # See  what filters we're going to use:
   if self.filter_order is not None:
      bands = self.filter_order
   else:
      bands = list(self.data.keys())
      if not single:
         eff_wavs = []
         for filter in bands:
            eff_wavs.append(fset[filter].ave_wave)
         eff_wavs = asarray(eff_wavs)
         ids = argsort(eff_wavs)
         self.filter_order = [bands[i] for i in ids]
      else:
         mins = []
         for filter in bands:
            if getattr(self.data[filter],'model',None) is not None:
               mins.append(min(self.data[filter].model))
            else:
               mins.append(min(self.data[filter].mag))
         mins = asarray (mins)
         ids = argsort(mins)
         self.filter_order = [bands[i] for i in ids]
      bands = self.filter_order
      
   for b in bands:
      if b not in colors:  colors[b] = 1
      if b not in symbols:  symbols[b] = 4
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
      self.p = pygplot.Panel(cols, rows, device='/NULL', aspect=1.0*rows/cols)
   else:
      self.p = pygplot.Plot(device='/NULL', title=title, flipyaxis=flip,
            xrange=xrange, yrange=yrange, fsize=fsize, linewidth=linewidth,
            width=linewidth, xlabel='Epoch (days)', ylabel=ylabel, **pargs)

   i = 0
   for filter in bands:
      if not single:
         pp = pygplot.Plot(device='/NULL', aspect=1.0,  flipyaxis=flip,
            xrange=xrange, yrange=yrange, linewidth=linewidth, font=2, 
            fsize=fsize, width=linewidth, **pargs)
         # make a reference to the lightcruve we are plotting.
         pp.lc = self.data[filter]
      else:
         pp = self.p
      if not i and title is None:
         label = self.name + " " + filter
      else:
         label = filter
      if not single:
         pp.label(0.8, 0.9, label, reference='relative', color=colors[filter], 
               just=1, fsize=fsize, vjust=1)
      if mask:
         x = self.data[filter].MJD[self.data[filter].mask]
         if not flux:
            y = self.data[filter].mag[self.data[filter].mask] + \
                  single*i*dm - relative*rel_off
            ey = self.data[filter].e_mag[self.data[filter].mask]
         else:
            y = self.data[filter].flux[self.data[filter].mask]* \
                  power(10, -0.4*(single*i*dm-relative*rel_off))
            ey = self.data[filter].e_flux[self.data[filter].mask]
      else:
         x = self.data[filter].MJD
         if not flux:
            y = self.data[filter].mag + single*i*dm - relative*rel_off
            ey = self.data[filter].e_mag
         else:
            y = self.data[filter].flux* \
                  power(10, -0.4*(single*i*dm-relative*rel_off))
            ey = self.data[filter].e_flux

      # Make the PANIC data stand out
      #if filter=='Yc' or filter == 'Jc':
      #   size = fsize + 2
      #else:
      #   size = fsize + 1
      if self.Tmax is None:  
         Tmax = 0
      else:
         Tmax = self.Tmax
      pp.point(x-Tmax*epoch, y, symbol=symbols[filter], fill=1,
            fillcolor=colors[filter], size=fsize+1, label=filter+'+'+str(i*dm))
      pp.error(x-Tmax*epoch, y, dy1=ey, length=0, order=1000)
      if label_bad:
         gids = equal(self.data[filter].mask, 0)
         if any(gids):
            x = self.data[filter].MJD[gids] - Tmax*epoch
            if not flux:
               y = self.data[filter].mag[gids] + \
                     single*i*dm - relative*rel_off
            else:
               y = self.data[filter].flux[gids]* \
                     power(10, -0.4*(single*i*dm-relative*rel_off))
            pp.point(x, y, symbol=5, color='red', size=4.0)

      # Now check to see if there is a model to plot:
      if self.model.Tmax is not None and filter in self.model._fbands:
         t = arange(-10,70,1.0) + self.Tmax
         mag,err,gids = self.model(filter, t)
         if not flux:
            y = mag + i*single*dm - relative*rel_off
         else:
            zp = fset[filter].zp
            y = power(10, -0.4*(mag - zp + i*single*dm - \
                                relative*rel_off))
         pp.line(compress(gids,t-self.Tmax*epoch), compress(gids,y), 
               color=colors[filter])
      elif self.data[filter].tck is not None:
         tck = self.data[filter].tck 
         t = arange(tck[0][0], tck[0][-1], 1.0)
         mag,gids = self.data[filter].eval(t, t_tol=-1)
         if not flux:
            y = mag + i*single*dm - relative*rel_off
         else:
            zp = fset[filter].zp
            y = power(10, -0.4*(mag - zp + i*single*dm- \
                                relative*rel_off))
         if self.Tmax is None:
            Tmax = 0
         else:
            Tmax = self.Tmax
         pp.line(compress(gids,t-Tmax*epoch), compress(gids,y),
                  color=colors[filter])
      i = i + 1

      if not single:  
         self.p.add(pp)
      else:
         if legend:  self.p.legend(position='lr', fsize=1.5, bbox=(0,0))
         self.p.label(0.05,0.05,self.name,reference='relative')
   
   self.p.plot()

   if self.xrange is None and not single:
      xmin = self.p.plots[0].xmin
      xmax = self.p.plots[0].xmax
      for plot in self.p.plots[1:]:
         xmin = min(plot.xmin, xmin)
         xmax = max(plot.xmax, xmax)
      for plot in self.p.plots:  plot.xrange = [xmin,xmax]

   if self.yrange is None and not single:
      ymin = self.p.plots[0].ymin
      ymax = self.p.plots[0].ymax
      for plot in self.p.plots[1:]:
         if not flux:
            ymin = max(plot.ymin, ymin)
            ymax = min(plot.ymax, ymax)
         else:
            ymin = min(plot.ymin, ymin)
            ymax = max(plot.ymax, ymax)
      for plot in self.p.plots:  plot.yrange = [ymin,ymax]

   self.p.close()
   self.p.device=device
   self.p.plot()
   if interactive:  self.p.interact()
   self.p.close()
   return(self.p)

def mask_data(self):
   p = self.plot(single=0, label_bad=1, epoch=1, mask=0)
   for b in self.data:
      # a place to store red x's
      self.data[b].pts = {}
   p.plot()
   p.interact(bindings={'A':bind_mask, 'u':bind_mask})
   p.close()


def bind_mask(x,y,key,inst):
   '''Call-back function used to select data when using the mask_data() 
   member function.'''
   lc = inst.lc

   # now, what are we doing?
   if key == 'A':
      # first, if all data is masked:
      if not any(lc.mask):
         return
      # find distances to valid (unmasked) data:
      dists = power(lc.mag[lc.mask] - y, 2) + power(lc.t[lc.mask] - x, 2)
      i = argmin(dists)
      # need to find index in original (masked+unmasked) data
      id = indices(lc.mask.shape)[0][lc.mask][i]
      lc.mask[id] = 0
      lc.pts[id] = inst.point([lc.t[id]],[lc.mag[id]], symbol=5, color='red', size=4.0)
   elif key == 'u':
      # first, if none are masked, nothing to do
      if all(lc.mask):
         return
      # find distances to invalid (masked) data:
      gids = logical_not(lc.mask)
      dists = power(lc.mag[gids] - y, 2) + power(lc.t[gids] - x, 2)
      i = argmin(dists)
      # need to find index in original (masked+unmasked) data
      id = indices(lc.mask.shape)[0][gids][i]
      lc.mask[id] = 1
      pid = inst.points.index(lc.pts[id])
      del inst.points[pid]
      

def plot_lira(t, t2, t_maxes, BV, eBV, BV2, tmin, tmax, c):

    p = pygplot.Plot(device='9/XS', font=2, fsize=1.5, xlabel='Epoch (Vmax)', 
          ylabel='B-V')
    p.point(t-t_maxes[0], BV, symbol=4, fill=1)
    p.error(t-t_maxes[0], BV, dy1=eBV, length=0)
    p.point(t2, BV2, symbol=4, fill=1, fillcolor='red', label='B-V (deredshifted)')
    p.line([tmin, tmax], [0.732-0.0095*(tmin-55), 0.732-0.0095*(tmax-55)], 
          color='blue', label='Lira Law')
    p.line([tmin, tmax], [c[0]+c[1]*(tmin-55), c[0]+c[1]*(tmax-55)], color='red',
          label='Fit')
    p.legend(position='ur', fsize=1)
    p.plot()
    p.close()


def xrange_select(x, y, key, inst, tol=1e-6):
   x2,y2,key = pygplot.ppgplot.pgband(4,0,x,y)
   if len(inst.parent.plots) > 1:
      inst2 = inst.parent.plots[1]
   else:
      inst2 = None
   if inst.xrange is None:
      scale = abs(inst._xmax - inst._xmin)
   else:
      scale = abs(inst.xrange[1] - inst.xrange[0])
   if abs(x - x2)/scale < tol:
      inst.xrange = None
      if inst2 is not None:  inst2.xrange = None
   else:
      inst.xrange = [x,x2]
      if inst2 is not None:  inst2.xrange = [x,x2]

def box_range_select(x,y,key,inst, tol=1.e-6):
   '''Select an x and y range simultaneously.'''
   x2,y2,key = pygplot.ppgplot.pgband(2,0,x,y)
   if len(inst.parent.plots) > 1:
      inst2 = inst.parent.plots[1]
   else:
      inst2 = None
   
   if inst.xrange is None:
      xscale = abs(inst._xmax - inst._xmin)
   else:
      xscale = abs(inst.xrange[1] - inst.xrange[0])
   if inst.yrange is None:
      yscale = abs(inst._ymax - inst._ymin)
   else:
      yscale = abs(inst.yrange[1] - inst.yrange[0])
   if abs(y - y2)/yscale < tol or abs(x - x2)/xscale < tol:
      inst.xrange = inst.yrange = None
      if inst2 is not None:  inst2.xrange = None
   else:
      inst.xrange = [x,x2]
      if inst2 is not None:  inst2.xrange = [x,x2]
      inst.yrange = [y,y2]


def add_knot(x, y, key, inst):
   '''Add a knot point to the spline.  This is called through the interactive interface).'''
   self = inst.lc_inst
   if self.model_type is None or self.model_type != 'spline':
      print("Adding knots not supported for this model")
      return
   x,y = lc_data_conv(inst, x, y)
   if x <= self.tck[0][0] or x >= self.tck[0][-1]:
      print("Error, can't add knots outside the current range.")
   old_knots = self.tck[0]
   print(old_knots)
   knots = concatenate([[x],self.tck[0]])
   knots = sort(knots)
   self.tck = (knots, self.tck[1], self.tck[2])
   print(knots)
   try:
      update_plot(self, task=-1)
   except:
      self.tck[0] = old_knots
      update_plot(self, task=-1)

def move_knot(x, y, key, inst):
   '''Move the point closest to x to a new location.'''
   self = inst.lc_inst
   if self.model_type is None or self.model_type != 'spline':
      print("Moving knots not supported for this model")
      return
   x,y = lc_data_conv(inst, x, y)
   old_knots = self.tck[0]
   deltas = absolute(self.tck[0] - x)
   id = argmin(deltas)
   x2,y2,key = pygplot.ppgplot.pgband(4,0,x,y)
   self.tck[0][id] = x2
   try:
      update_plot(self, task=-1)
   except:
      self.tck[0] = old_knots
      update_plot(self, task=-1)

def delete_knot(x, y, key, inst):
   self = inst.lc_inst
   if self.model_type is None or self.model_type != 'spline':
      print("Deleting knots not supported for this model")
      return
   x,y = lc_data_conv(inst, x, y)
   old_knots = self.tck[0]
   deltas = absolute(self.tck[0] - x)
   id = argmin(deltas)
   list = self.tck[0].tolist()
   del list[id]
   self.tck = (array(list), self.tck[1], self.tck[2])
   try:
      update_plot(self, task=-1)
   except:
      self.tck[0] = old_knots
      update_plot(self, task=-1)

def change_s(x, y, key, inst):
   self = inst.lc_inst
   if self.model_type is None or self.model_type != 'spline':
      print("Changing s not supported for this model")
      return
   if self.tck is not None:
      if key == 's':
         self.s -= sqrt(2*self.s)
         if self.s < 0:  self.s = 0
      else:
         self.s += sqrt(2*self.s)
      update_plot(self, task=0)
   else:
      if self.line.N is not None:  self.line.N = None
      if self.line.sigma is None:  self.line.sigma=3.0
      if key == 's':
         self.line.sigma *= 1.1
      else:
         self.line.sigma /= 1.1
      update_plot(self, task=-2)

def change_N(x, y, key, inst):
   self = inst.lc_inst
   if self.line.sigma is not None:  self.line.sigma = None
   if self.line.N is None:  self.line.N = 5
   if key == 'n':
      if self.line.N > 1:  self.line.N -= 1
   else:
      self.line.N += 1
   update_plot(self, task=-2)

def update_plot(self, task):
   '''Update the plot with new spline info.'''
   inst = self.mp
   if len(inst.plots) != 2:
      return
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
      inst.plots[0].labels[0].string = 'N = %d   sigma = %.1f   r-chi-sq = %.2f' % \
            (N,sigma,self.line.chisq())
   else:
      inst.plots[0].labels[0].string = 'Np = %d  Nk = %d  s = %.1f r-chi-sq = %.2f' % \
            (len(self.mag), len(self.tck[0]), self.s, self.rchisq)
   if task != -2:
      del inst.plots[0].points[1]
   if not inst.plots[0].flipyaxis:
      y = compress(self.mask, self.flux - power(10, -0.4*(m_model-self.filter.zp)))
      inst.plots[0].lines[0].y = power(10, -0.4*(m - self.filter.zp))
      inst.plots[1].points[0].y = y
      inst.plots[1].errors[0].y = y
      if task != -2:
         inst.plots[0].point(self.tck[0], power(10, -0.4*(self.eval(self.tck[0], t_tol=None)[0] - \
               self.filter.zp)), color='red', symbol=4)
   else:
      y = compress(self.mask, self.mag - m_model)
      inst.plots[0].lines[0].y = m
      inst.plots[1].points[0].y = y
      inst.plots[1].errors[0].y = y
      if task != -2:
         inst.plots[0].point(self.tck[0] - self._epoch*self.parent.Tmax, 
               self.eval(self.tck[0], t_tol=None)[0], color='red', 
               symbol=4)
  
def lc_data_conv(p, x, y):
   '''Given a plot instance p and (x,y) coordinates, convert to the real
   data (which can be different because of plotting as magnitudes, flux,
   epoch or no, etc.'''
   # first, if plotting in epochs:
   realx = x + p.epoch*p.lc_inst.parent.Tmax
   if not p.flux:
      realy = power(10, -0.4*(y - p.lc_inst.filter.zp))
   return(realx,realy)

def plot_lc(self, device='/XSERVE', interactive=0, epoch=1, flux=0, gloes=0,
          symbol=4):
   '''Plot this light-curve.  You can specify a PGPLOT device, the default is an X
   server.  If interactive=1, basic pygplot keybindings are available.  If flux=1,
   plot in flux space.  If epoch=1, plot time relative to Tmax.  If gloes=1, 
   use GLoEs to smooth the data and produce a model.  You can specify the
   symbol to plot with 'symbol'.'''
   if flux: 
      flipaxis = 0
   else:
      flipaxis = 1
   p = pygplot.Plot(flipyaxis=flipaxis, title='%s %s lightcurve' % (self.parent.name,
      self.band), xlabel='Epoch (days)', ylabel='mag', font=2, leftpad=0.5)
   p.lc_inst = self
   p.epoch = epoch
   p.flux = flux
   if self.parent.Tmax is None:
      Tmax = 0
   else:
      Tmax = self.parent.Tmax

   if gloes or self.model_type is not None or self.band in self.parent.model._fbands:
      self.mp = pygplot.Panel(device=device, llcs=[(0.0,0.2),(0.0,0.0)],
            urcs=[(1.0,1.0),(1.0,0.2)], justs=[12,13])
      p2 = pygplot.Plot(flipyaxis=flipaxis, xlabel='Epoch (days)', font=2,
            ylabel='residuals', leftpad=0.5)
      p2.lc_inst = self
   else:
      self.mp = pygplot.Panel(1,1,device=device)
   self.mp.add(p)
   if gloes or self.model_type is not None or self.band in self.parent.model._fbands:
      self.mp.add(p2)
   if flux:
      y = power(10, -0.4*(self.mag - self.filter.zp))
      ey = y*self.e_mag/1.0857
   else:
      y = self.mag
      ey = self.e_mag
   p.point(self.MJD - epoch*Tmax, y, symbol=symbol)
   p.error(self.MJD - epoch*Tmax, y, dy1=ey)

   # Order of preference:  plot models, else plot splines, else plot
   #  nothing unless we explicitly ask for gloess.
   if self.band in self.parent.model._fbands:
      t = arange(-10,70,1.0) + Tmax
      m,em,mask = self.parent.model(self.band, t)
      m_m,m_em,m_mask = self.parent.model(self.band, self.MJD)
      x = self.MJD[self.mask*m_mask]
      if flux and self.parent.model.model_in_mags:
         p.line(t[mask]-epoch*Tmax, power(10, -0.4*(m-self.filter.zp)))
         y = compress(self.mask*m_mask, 
               self.flux - power(10,-0.4*(m_m - self.filter.zp)))
         dy = self.e_flux[self.mask*m_mask]
      elif not flux and not self.parent.model.model_in_mags:
         p.line(t[mask]-epoch*Tmax, -2.5*log10(m[mask]) + self.filter.zp)
         y = compress(self.mask*m_mask, self.mag + 2.5*log10(m_m) + self.filter.zp)
         dy = self.e_mag[self.mask*m_mask]
      elif flux:
         p.line(t[mask]-epoch*Tmax,m[mask])
         y = self.flux[self.mask*m_mask] - m_m[self.mask*m_mask]
         dy = self.e_flux[self.mask*m_mask]
      else:
         p.line(t[mask]-epoch*Tmax,m[mask])
         y = self.mag[self.mask*m_mask] - m_m[self.mask*m_mask]
         dy = self.e_mag[self.mask*m_mask]

      p2.point(x - epoch*Tmax, y, symbol=symbol)
      p2.error(x- epoch*Tmax, y, dy1=dy)
      p2.line([t[0]-epoch*Tmax,t[-1]-epoch*Tmax], [0,0])

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
         p.line(compress(mask,t - epoch*Tmax), 
                compress(mask, power(10, -0.4*(m - self.filter.zp))))
         if self.tck is not None:
            p.point(self.tck[0] - epoch*Tmax, power(10, -0.4*(self.eval(self.tck[0], t_tol=None)[0] \
               - self.filter.zp)), color='red', symbol=4)
         x = compress(self.mask*m_mask, self.MJD)
         y = compress(self.mask*m_mask, self.flux - power(10, -0.4*(m_model-self.filter.zp)))
         dy = compress(self.mask*m_mask, self.e_flux)
         p2.point(x - epoch*Tmax,y, symbol=symbol)
         p2.error(x - epoch*Tmax,y, dy1=dy)
      else:
         p.line(compress(mask,t - epoch*Tmax), compress(mask,m))
         if self.tck is not None:
            p.point(self.tck[0] - epoch*Tmax, self.eval(self.tck[0], t_tol=None)[0], 
                  color='red', symbol=4)
         x = compress(self.mask*m_mask, self.MJD)
         y = compress(self.mask*m_mask, self.mag - m_model)
         p2.point(x - epoch*Tmax,y, symbol=symbol)
         dy = compress(self.mask*m_mask, self.e_mag)
         p2.error(x - epoch*Tmax,y,dy1=dy)
      if self.tck is not None:
         p2.line(array([self.tck[0][0],self.tck[0][-1]]) - epoch*Tmax, [0,0])
      else:
         p2.line(array([self.MJD[0],self.MJD[-1]]) - epoch*Tmax, [0,0])
      # useful stats:
      if self.tck is not None:
         p.label(0.5, 1.03, 'Np = %d  Nk = %d  s = %.1f r-chi-sq = %.2f' % \
            (len(self.mag), len(self.tck[0]), self.s, self.rchisq), 
            reference='relative', fsize=1, just=0.5)
      else:
         if self.line.N is not None:  N = self.line.N
         else: N = -1
         if self.line.sigma is not None:  sigma = self.line.sigma
         else: sigma = -1
         p.label(0.5, 1.03, 'N = %d   sigma = %.1f   r-chi-sq = %.2f' % \
               (N,sigma,self.line.chisq()), reference='relative', fsize=1, just=0.5)


   self.mp.plot()
   if len(self.mp.plots) == 2:
      self.mp.plots[1].xrange = [self.mp.plots[0].xmin, self.mp.plots[0].xmax]
      self.mp.replot()
   if interactive:
      self.mp.interact(bindings={'a':add_knot, 'd':delete_knot, 
            'm':move_knot, 's':change_s, 'S':change_s,
            'n':change_N, 'N':change_N, 'x':xrange_select,
            'b':box_range_select})
   self.mp.close()
   return(self.mp)
   
def plot_kcorrs(self, device='13/XW', colors=None, symbols=None):
   '''Plot the k-corrections, both mangled and un-mangled.'''
   # See  what filters we're going to use:
   bands = list(self.ks.keys())
   if len(bands) == 0:
      return

   if colors is None:  colors = default_colors
   if symbols is None:  symbols = default_symbols

   eff_wavs = []
   for filter in bands:
      eff_wavs.append(fset[filter].ave_wave)
   eff_wavs = asarray(eff_wavs)
   ids = argsort(eff_wavs)
   bands = [bands[i] for i in ids]
   minx = Inf
   maxx = -Inf
   for b in bands:
      if b not in colors:  colors[b] = 1
      if b not in symbols:  symbols[b] = 4
      minx = min(minx, self.data[b].t.min())
      maxx = max(maxx, self.data[b].t.max())
   n_plots = len(bands)
   cols = int(round(sqrt(n_plots)))
   rows = (n_plots / cols)
   if n_plots % cols:  rows += 1
   p = pygplot.Panel(1, n_plots, device=device, aspect=1.0/6*n_plots)

   for b in bands:
      x = self.data[b].MJD - self.Tmax
      days = arange(int(x[0]), int(x[-1]), 1)
      rest_days = days/(1+self.z)/self.ks_s
      k,k_m = list(map(array, kcorr.kcorr(rest_days.tolist(), self.restbands[b], b, 
         self.z, self.EBVgal, 0.0, version=self.k_version)))
      k_m = equal(k_m, 1)
      pp = pygplot.Plot(fsize=1, font=2, xlabel='Epoch (days)',
            leftpad=0.1, rotylab=1, xrange=[minx,maxx])
      pp.filter = b
      pp.inst = self
      p.add(pp)
      pp.line(days[k_m],k[k_m], color=colors[b])
      pp.point(x[self.ks_mask[b]], self.ks[b][self.ks_mask[b]], 
            symbol=symbols[b], color=colors[b], fillcolor=colors[b],
            fsize=2)
      pp.label(0.9, 0.9, b, just=1, vjust=1, reference='relative')

   p.plots[n_plots/2].ylabel = 'K-correction'
   p.plot()
   p.interact(bindings={'A':plot_mangled_SED})
   p.close()
   return(p)

def plot_mangled_SED(x, y, key, inst):

   p = pygplot.Plot(device='14/XS', font=1.5, fsize=2, xlabel='Wavelength (Angstroms)',
         ylabel='Flux')
   band = inst.filter
   self = inst.inst
   id = argmin(power(x - self.data[band].t,2) + power(y - self.ks[band],2))
   day = int(self.data[band].t[id])
   wave,flux = kcorr.get_SED(day, version='H3')
   p.line(wave,flux, label='Original SED')

   if band in self.ks_mopts:
      man_flux = mangle_spectrum.apply_mangle(wave,flux, **self.ks_mopts[band][id])
      p.line(wave, man_flux, label='Mangled SED', color='007700')
   f1 = fset[band]
   f2 = fset[self.restbands[band]]
   p.line(f1.wave, f1.resp/f1.resp.max()*flux.max(), color='red')
   p.line(f2.wave*(1+self.z), f2.resp/f2.resp.max()*flux.max(), 
         color='blue')

   p.legend(fsize=1.5)
   p.xrange = [fset[band].wave.min()*0.9, fset[band].wave.max()*1.1]
   p.plot()
   p.close()
