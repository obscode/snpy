'''This class lets you interactively fit using the fit1dcurve
classes.'''
import numpy as num
import myplotlib
from myplotlib import PanelPlot

class InteractiveFit:

   def __init__(self, interp, title=None, figsize=(12,8)):
      '''Takes a oneDcurve instance as argument.  Sets up a plot with
      fit and residuals, then sets up bindings to interactively fit.'''

      self.interp = interp
      self.x = interp.xdata*1
      self.y = interp.ydata*1
      self.ey = interp.eydata*1
      self.mask = interp.mask.copy()
      
      self.title = title
      self.figsize = figsize

      self.setup_graph()
      self.set_bindings()

   def __getattr__(self, key):
      if 'interp' in self.__dict__:
         if key in self.interp.pars:
            return self.interp.pars[key]
      raise AttributeError, "no such attribute %s" % key

   def __setattr__(self, key, value):
      if 'interp' in self.__dict__:
         if key in self.interp.pars:
            setattr(self.interp,key,value)
            self.redraw()
            return
      self.__dict__[key] = value


   def setup_graph(self):

      self.mp = PanelPlot(1,2, pheights=[0.2,0.8], pwidths=[0.8],
            figsize=self.figsize)

      if self.title is None:
         self.mp.title("Fitting %s\ntype '?' for help" % str(self.interp))
      else:
         self.mp.title(self.title + "\ntype '?' for help")
      self.mp.xlabel("X")
      self.mp.axes[0].set_ylabel("residuals")
      self.mp.axes[1].set_ylabel("Y")

      # fit plot
      self.mp.axes[1].errorbar(self.x, self.y, yerr=self.ey, capsize=0,
            fmt='o')
      if not num.alltrue(self.mask):
         self._x1, = self.mp.axes[1].plot(self.x[-self.mask], self.y[-self.mask], 'x', 
               color='red', ms=16)
      else:
         self._x1 = None
      xmin,xmax = self.interp.domain()
      xs = num.linspace(xmin,xmax, 100)
      ys,m = self.interp(xs)
      self._mod, = self.mp.axes[1].plot(xs[m], ys[m], '-', color='black')

      # residuals plot
      ys,m = self.interp(self.x)
      resids = self.interp.residuals(mask=False)
      self._rp,dum,self._rl = self.mp.axes[0].errorbar(self.x, 
            resids, yerr=self.ey, capsize=0, fmt='o')
      if not num.alltrue(self.interp.mask):
         self._x2, = self.mp.axes[0].plot(self.x[-self.mask], 
               self.interp.residuals(mask=False)[-self.mask], 'x', ms=16, color='red')
      else:
         self._x2 = None
      self.mp.axes[0].axhline(0, color='black')

      self.mp.set_limits()
      if not num.alltrue(m):
         resids = resids[m]
         self.mp.axes[0].set_ylim((resids.min(), resids.max()))
      self.mp.draw()

   def redraw(self):
      '''Redraw the graph'''
      # Only change those things that need changing
      self._mod.remove()
      self._rp.remove()
      self._rl[0].remove()

      self.interp.mask = self.mask.copy()
      self.interp.setup = False

      # Update the model
      xmin,xmax = self.interp.domain()
      xs = num.linspace(xmin,xmax, 100)
      ys,m = self.interp(xs)
      self._mod, = self.mp.axes[1].plot(xs[m], ys[m], '-', color='black')

      # Update the residuals
      ys,m = self.interp(self.x)
      resids = self.interp.residuals(mask=False)
      self._rp,dum,self._rl = self.mp.axes[0].errorbar(self.x, 
            resids, yerr=self.ey, capsize=0, fmt='o', color='blue')
      self.mp.set_limits()
      if not num.alltrue(m):
         resids = resids[m]
         self.mp.axes[0].set_ylim((resids.min(), resids.max()))
      self.redraw_x()

   def help(self):
      self.interp.help()

   def redraw_x(self):
      '''Redraw the little red X's if needed'''
      if not num.alltrue(self.mask):
         if self._x1 is not None:
            self._x1.set_data((self.x[-self.mask], self.y[-self.mask]))
         else:
            self._x1, = self.mp.axes[1].plot(self.x[-self.mask], 
                  self.y[-self.mask], 'x', color='red', ms=16)
         if self._x2 is not None:
            self._x2.set_data((self.x[-self.mask], self.interp.residuals(mask=False)[-self.mask]))
         else:
            self._x2, = self.mp.axes[0].plot(self.x[-self.mask], 
                  self.interp.residuals(mask=False)[-self.mask], 'x', color='red', ms=16)
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
         id = num.argmin(num.power(x-self.x,2) + num.power(y-self.interp.residuals(mask=False),2))
      else:
         id = num.argmin(num.power(x-self.x,2) + num.power(y-self.y,2))
      self.mask[id] = -self.mask[id]
      self.redraw_x()

   def _bind_r(self, event):
      '''Binding for re-fitting and drawing.'''
      if event.key != 'r':
         return
      self.redraw()

   def _bind_q(self, event):
      '''Binding for quitting the session.'''
      if event.key == 'q':
         self.mp.close()

   def _bind_help(self, event):
      '''Binding for help'''
      if event.key == '?':
         self.help()

   def set_bindings(self):
      '''Sets the bindings to the figure canvas.'''
      self.mp.fig.canvas.mpl_connect('key_press_event',self._bind_x)
      self.mp.fig.canvas.mpl_connect('key_press_event',self._bind_r)
      self.mp.fig.canvas.mpl_connect('key_press_event',self._bind_q)
      self.mp.fig.canvas.mpl_connect('key_press_event',self._bind_help)

