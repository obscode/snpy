'''This class lets you interactively fit using the fit1dcurve
classes.'''
from __future__ import print_function
import numpy as num
from snpy.myplotlib import PanelPlot
from matplotlib import pyplot as plt
try:
   from . import fit1dcurve
except:
   from snpy.utils import fit1dcurve

def find_renderer(fig):
   if hasattr(fig.canvas, 'get_renderer'):
      renderer = fig.canvas.get_renderer()
   else:
      import io
      fig.canvas.print_pdf(io.BytesIO())
      renderer = fig._cachedRenderer
   return renderer

class InteractiveFit:

   def __init__(self, interp, title=None, figsize=None, fignum=None, draw=True,
         xlabel='X', ylabel='Y', extra_bindings={}):
      '''Takes a oneDcurve instance as argument.  Sets up a plot with
      fit and residuals, then sets up bindings to interactively fit.'''

      self.interp = interp
      self.x = interp.xdata*1
      self.y = interp.ydata*1
      self.ey = interp.eydata*1
      self.mask = interp.mask.copy()
      
      self.title = title
      self.xlabel = xlabel
      self.ylabel = ylabel
      self.figsize = figsize
      self.fignum = fignum

      # A callback to use when closing the figure
      self.extra_bindings = extra_bindings

      # Some help strings
      self.bind_help = [('r','redraw the graph'),
                        ('x','flag/unflag bad data')]
      if isinstance(self.interp, fit1dcurve.Spline):
         self.bind_help += [('a','Add knot point at cursor position'),
                            ('d','Delete knot closest to cursor position'),
                            ('m','Move knot closest to cursor to new position')]
      if isinstance(self.interp, fit1dcurve.Polynomial):
         self.bind_help += [('n/N','Decrease/Increase order of polynomial'),
                            ('m','specify range over which to fit')]
      if fit1dcurve.GaussianProcess is not None and \
            isinstance(self.interp, fit1dcurve.GaussianProcess):
         self.bind_help += [('s/S', 'Decrease/Increase scale by 10%'),
                            ('a/A', 'Decrease/Increase amplitude by 10%'),
                            ('d/D', 'Decrease/Increase degree of diff. by 1')]
      self.bind_help.append(('q','close the graph'))

      self.setup_graph(draw=draw)
      self.set_bindings()

   def __getattr__(self, key):
      if 'interp' in self.__dict__:
         if key in self.interp.pars:
            return self.interp.pars[key]
      raise AttributeError("no such attribute %s" % key)

   def __setattr__(self, key, value):
      if 'interp' in self.__dict__:
         if key in self.interp.pars:
            setattr(self.interp,key,value)
            self.redraw()
            return
      self.__dict__[key] = value


   def setup_graph(self, draw=True):

      self.mp = PanelPlot(1,2, pheights=[0.2,0.6], #pwidths=[0.8],
            figsize=self.figsize, num=self.fignum) #, rect=(0,0,0.8,1))
      #self.mp.right_pad=0.20

      if self.title is None:
         self.mp.title("Fitting %s\ntype '?' for help" % str(self.interp))
      else:
         self.mp.title(self.title + "\ntype '?' for help")
      self.mp.xlabel(self.xlabel)
      self.mp.axes[0].set_ylabel("residuals")
      self.mp.axes[1].set_ylabel(self.ylabel)

      # fit plot
      self.mp.axes[1].errorbar(self.x, self.y, yerr=self.ey, capsize=0,
            fmt='o')
      if not num.all(self.mask):
         m = num.logical_not(self.mask)
         self._x1, = self.mp.axes[1].plot(self.x[m], self.y[m], 'x', 
               color='red', ms=16)
      else:
         self._x1 = None
      xmin,xmax = self.interp.domain()
      xs = num.linspace(xmin,xmax, 100)
      ys,m = self.interp(xs)
      self._mod, = self.mp.axes[1].plot(xs[m], ys[m], '-', color='black')

      # if spline, plot the knot points
      if isinstance(self.interp, fit1dcurve.Spline):
         tck = self.interp.tck
         self._knots, = self.mp.axes[1].plot(tck[0][4:-4], 
               self.interp(tck[0][4:-4])[0], 'd', color='red', zorder=1000)
      else:
         self._knots = None

      # residuals plot
      ys,m = self.interp(self.x)
      resids = self.interp.residuals(mask=False)
      self._rp,dum,self._rl = self.mp.axes[0].errorbar(self.x, 
            resids, yerr=self.ey, capsize=0, fmt='o')
      if not num.all(self.interp.mask):
         self._x2, = self.mp.axes[0].plot(self.x[num.logical_not(self.mask)], 
               self.interp.residuals(mask=False)[num.logical_not(self.mask)],
               'x', ms=16, color='red')
      else:
         self._x2 = None
      self.mp.axes[0].axhline(0, color='black')

      self.mp.set_limits()
      if not num.all(m):
         resids = resids[m]
         self.mp.axes[0].set_ylim((resids.min(), resids.max()))
      self.mp.axes[1].text(1.01, 1.0, "Parameters:", va='bottom', ha='left',
            transform=self.mp.axes[1].transAxes, color='blue')
      self.mp.axes[0].text(1.01, 1.0, "Stats:", va='bottom', ha='left',
            transform=self.mp.axes[0].transAxes, color='blue')
      pars_labels = ""
      for par in self.interp.pars:
         if par == 't': continue
         pars_labels += par + ":  \n"
      self._pars_labels = self.mp.axes[1].text(1.03, 0.95, pars_labels,
            ha='left', va='top', transform=self.mp.axes[1].transAxes,
            fontdict={'size':10})
      self._stats_labels = self.mp.axes[0].text(1.03, 0.95, 
          "RMS:\nr-chi-sq:\nDW:", va='top', ha='left', 
          transform=self.mp.axes[0].transAxes, fontdict={'size':10})
      self.mp.draw()
      self.plot_stats()
      self.plot_params()
      self.mp.draw()


   def draw(self):
      '''Just an alias to self.mp.draw.'''
      self.mp.draw()
      
   def redraw(self):
      '''Redraw the graph'''
      x1,x2 = self.mp.axes[0].get_xlim()
      # Only change those things that need changing
      self._mod.remove()
      self._rp.remove()
      self._rl[0].remove()
      if self._knots:
         self._knots.remove()

      self.interp.mask = self.mask.copy()
      self.interp.setup = False

      # Update the model
      xmin,xmax = self.interp.domain()
      xs = num.linspace(xmin,xmax, 1000)
      ys,m = self.interp(xs)
      self._mod, = self.mp.axes[1].plot(xs[m], ys[m], '-', color='black')

      # Update the knots if they exist:
      if self._knots is not None:
         tck = getattr(self.interp, 'tck', None)
         if tck:
            self._knots, = self.mp.axes[1].plot(tck[0][4:-4], 
               self.interp(tck[0][4:-4])[0], 'd', color='red', zorder=1000)

      # Update the residuals
      ys,m = self.interp(self.x)
      resids = self.interp.residuals(mask=False)
      self._rp,dum,self._rl = self.mp.axes[0].errorbar(self.x, 
            resids, yerr=self.ey, capsize=0, fmt='o', color='blue')
      #self.mp.set_limits()
      if not num.all(m):
         resids = resids[m]
         self.mp.axes[0].set_ylim((resids.min(), resids.max()))
      self.plot_stats()
      self.plot_params()
      self.redraw_x()

   def help(self):
      print("The following key bindings are in effect:")
      for key,help in self.bind_help:
         print(key+":   "+help)
      print()
      print("You can also vary the following member variables:")
      self.interp.help()

   def redraw_x(self):
      '''Redraw the little red X's if needed'''
      if not num.all(self.mask):
         if self._x1 is not None:
            self._x1.set_data((self.x[num.logical_not(self.mask)], 
               self.y[num.logical_not(self.mask)]))
         else:
            self._x1, = self.mp.axes[1].plot(self.x[num.logical_not(self.mask)],
                  self.y[num.logical_not(self.mask)], 'x', color='red', ms=16)
         if self._x2 is not None:
            self._x2.set_data((self.x[num.logical_not(self.mask)], 
               self.interp.residuals(mask=False)[num.logical_not(self.mask)]))
         else:
            self._x2, = self.mp.axes[0].plot(self.x[num.logical_not(self.mask)],
               self.interp.residuals(mask=False)[num.logical_not(self.mask)], 
               'x', color='red', ms=16)
      elif self._x1 is not None:
         self._x1.remove()
         self._x1 = None
         self._x2.remove()
         self._x2 = None
      self.mp.draw()

   def _bind_x(self, event):
      '''Binding for removing points.'''
      if event.key != 'x':
         return
      x,y = event.xdata,event.ydata
      ax = event.inaxes
      if ax is self.mp.axes[0]:
         id = num.argmin(num.power(x-self.x,2) + \
               num.power(y-self.interp.residuals(mask=False),2))
      else:
         id = num.argmin(num.power(x-self.x,2) + num.power(y-self.y,2))
      self.mask[id] = num.logical_not(self.mask[id])
      self.redraw_x()

   def _bind_r(self, event):
      '''Binding for re-fitting and drawing.'''
      if event.key != 'r':
         return
      if 'r' in self.extra_bindings:
         self.extra_bindings['r'](event)
      self.redraw()

   def _bind_q(self, event):
      '''Binding for quitting the session.'''
      if event.key == 'q':
         if 'q' in self.extra_bindings:
            self.extra_bindings['q'](event)
         self.mp.close()

   def _bind_help(self, event):
      '''Binding for help'''
      if event.key == '?':
         self.help()

   def plot_stats(self):
      '''plots the statistics of the fit.'''
      label = " %.2f\n %.2f\n %.2f" % \
            (self.interp.rms(), self.interp.rchisquare(), self.interp.DW())
      id = getattr(self, '_statsid', None)
      ax = self.mp.axes[0]
      if id is None:
         bbox = self._stats_labels.get_window_extent(find_renderer(self.mp.fig))
         #bbox = bbox.inverse_transformed(ax.transAxes)
         bbox = bbox.transformed(ax.transAxes.inverted())
         self._statsid = ax.text(bbox.x1,0.95,label, va='top', ha='left',
              transform=ax.transAxes, fontdict={'size':10})
      else:
         id.set_text(label)

   def plot_params(self):
      label = ""
      for par in self.interp.pars:
         if par == 't':  continue
         this_lab = str(self.interp.pars[par])
         if len(this_lab) > 6:
            this_lab = "{:.3f}".format(self.interp.pars[par])
         label += " "+this_lab+"\n"
      id = getattr(self, '_parsid', None)
      ax = self.mp.axes[1]
      if id is None:
         bbox = self._pars_labels.get_window_extent(find_renderer(self.mp.fig))
         #bbox = bbox.inverse_transformed(ax.transAxes)
         bbox = bbox.transformed(ax.transAxes.inverted())
         self._parsid = ax.text(bbox.x1, 0.95, label, va='top', ha='left',
               transform=ax.transAxes, fontdict={'size':10})
      else:
         id.set_text(label)


   def _move_knot_point(self, event):
      id = getattr(self, '_knotid', None)
      if id is None:  return
      x,y = event.xdata,event.ydata
      if x < self.interp.tck[0][0]:
         x = self.interp.tck[0][0]*1.0001
      if x > self.interp.tck[0][-1]:
         x = self.interp.tck[0][-1]*0.9999
      xdata = self._knots.get_xdata()
      ydata = self._knots.get_ydata()
      xdata[id] = x
      ydata[id] = self.interp(x)[0] 
      self._knots.set_xdata(xdata)
      self._knots.set_ydata(ydata)
      plt.draw()

   def _bind_knot_keys(self, event):
      '''Binding adding (a) and deleting (d) knot points.'''
      tck = getattr(self.interp, 'tck', None)
      if tck is None:  return
      ax = event.inaxes
      if ax is not self.mp.axes[1]:  return
      x,y = event.xdata,event.ydata
      if event.key == 'd':
         self.interp.delete_knot(x)
         self.redraw()
      elif event.key == 'a':
         self.interp.add_knot(x)
         self.redraw()
      elif event.key == 'm':
         mpressed = getattr(self, '_mpressed', False)
         if mpressed:
            self._mpressed = False
            knotold = tck[0][4:-4][self._knotid]
            self.interp.move_knot(knotold, x)
            self.mp.fig.canvas.mpl_disconnect(self._knot_move_id)
            self.redraw()
         else:
            self._mpressed = True
            id = num.argmin(num.absolute(x - tck[0][4:-4]))
            self._knotid = id
            self._knot_move_id = self.mp.fig.canvas.mpl_connect(
                  'motion_notify_event', self._move_knot_point)

   def _move_axvline(self, event):
      self._axvline1.set_xdata([event.xdata,event.xdata])
      plt.draw()

   def _bind_poly_keys(self, event):
      if event.key == 'n':
         if self.interp.n > 0:
            self.interp.n -= 1
            self.redraw()
      elif event.key == 'N':
         self.interp.n += 1
         self.redraw()
      elif event.key == 'm':
         if getattr(self,'_range_pending',None) is None:
            self._range_pending = True
            self._rx0 = event.xdata
            self._axvline0 = event.inaxes.axvline(self._rx0)
            self._axvline1 = event.inaxes.axvline(self._rx0)
            self._range_pending = self.mp.fig.canvas.mpl_connect(
                  'motion_notify_event', self._move_axvline)
         else:
            self._rx1 = event.xdata
            self.mp.fig.canvas.mpl_disconnect(self._range_pending)
            if absolute(self._rx0-self._rx1) == 0:
               self.interp.xmin = None
               self.interp.xmax = None
            else:
               self.interp.xmin = min(self._rx0, self._rx1)
               self.interp.xmax = max(self._rx0, self._rx1)
            self.interp.setup = False
            self._range_pending = None
            self._axvline0.remove()
            self._axvline1.remove()
            self.redraw()

   def _bind_gp_keys(self, event):
      '''Bindings when the Gaussian Process is in effect.'''
      if event.key == 'A':
         self.interp.amp += self._amp0*0.1
      elif event.key == 'a':
         if self.interp.amp > self._amp0*0.1:
            self.interp.amp -= self._amp0*0.1
      elif event.key == 's':
         if self.interp.scale > self._scale0*0.1:
            self.interp.scale -= self._scale0*0.1
      elif event.key == 'S':
         self.interp.scale += self._scale0*0.1
      elif event.key == 'd':
         if self.interp.diff_degree > 1:
            self.interp.diff_degree -= 1
      elif event.key == 'D':
         self.interp.diff_degree += 1
      else:
         return

      self.interp.setup = False
      self.redraw()

   def set_bindings(self):
      '''Sets the bindings to the figure canvas.'''
      # First clear out the default ones
      ids = list(self.mp.fig.canvas.callbacks.callbacks['key_press_event'].keys())
      for id in ids:
         self.mp.fig.canvas.mpl_disconnect(id)
      self.mp.fig.canvas.mpl_connect('key_press_event',self._bind_x)
      self.mp.fig.canvas.mpl_connect('key_press_event',self._bind_r)
      self.mp.fig.canvas.mpl_connect('key_press_event',self._bind_q)
      self.mp.fig.canvas.mpl_connect('key_press_event',self._bind_help)
      if isinstance(self.interp, fit1dcurve.Spline):
         self.mp.fig.canvas.mpl_connect('key_press_event', self._bind_knot_keys)
      if isinstance(self.interp, fit1dcurve.Polynomial):
         self.mp.fig.canvas.mpl_connect('key_press_event', self._bind_poly_keys)
      if fit1dcurve.GaussianProcess is not None and \
            isinstance(self.interp, fit1dcurve.GaussianProcess):
         self._scale0 = self.interp.scale
         self._amp0 = self.interp.amp
         self._diff_degree0 = self.interp.diff_degree
         self.mp.fig.canvas.mpl_connect('key_press_event', self._bind_gp_keys)


