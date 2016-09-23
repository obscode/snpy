'''
myplotlib:  a wrapper around the matplotlib plotting package.

GOALS: - to make matplotlib more intuitive (at least to me) and more like my own
         wrapper around PGPLOT:  pygplot
       - Make objects that mirror the matplotlib objects, but only the most
         useful subset (again, useful to me).
       - get rid of the get_* and set_* member functions and replace them with
         pythonic behaviour.
       - Make some basic classes for doing interactive graphics.

'''

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.text as text
from numpy import *
from matplotlib import ticker
from matplotlib import gridspec

from matplotlib import rcParams
rcParams['font.size'] = 18
rcParams['font.family'] = 'serif'
rcParams['xtick.major.size'] = 8
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.major.size'] = 8
rcParams['ytick.labelsize'] = 'large'
rcParams['xtick.minor.size'] = 4
rcParams['ytick.minor.size'] = 4
rcParams['xtick.major.pad'] = 10
rcParams['ytick.major.pad'] = 10
rcParams['legend.numpoints'] = 1

def scale_fonts(factor):
   rcParams['font.size'] *= factor
   #rcParams['font.family'] *= factor
   rcParams['xtick.major.size'] *=factor
   rcParams['ytick.major.size'] *=factor
   rcParams['xtick.minor.size'] *=factor
   rcParams['ytick.minor.size'] *=factor
   rcParams['xtick.major.pad'] *=factor
   rcParams['ytick.major.pad'] *=factor
   #rcParams['axes.labelsize'] *= factor

def line_bbox(line):
   '''extract a bounding box from a line2D object.'''
   bbox = transforms.Bbox.unit()
   bbox.update_from_data_xy(line.get_xydata())
   return(bbox)

def img_bbox(img):
   '''extract a bounding box form an image object.'''
   bbox = transforms.Bbox.unit()
   x0,x1,y0,y1 = img.get_extent()
   bbox.set_points(array([[x0,y0],[x1,y1]]))
   return(bbox)

def axis_bbox(axis):
   '''Given an axis, find the bounding box for all the data therein.'''
   bboxs = [line_bbox(line) \
         for line in axis.lines if line.get_transform() is axis.transData]
   bboxs += [patch.clipbox.inverse_transformed(axis.transData) \
              for patch in axis.patches]
   bboxs += [col.get_datalim(axis.transData) for col in axis.collections]
   bboxs += [img_bbox(img) for img in axis.images]
   if len(bboxs) >= 1:
      return transforms.Bbox.union(bboxs)
   else:
      return None

class SimplePlot:
   '''A simple plotting area that can either be used on its own or
   embedded in a more complex structure.'''

   def __init__(self, fig=None, position=None, nsubx=5, nsuby=5, **kwargs):

      self.nsubx = nsubx
      self.nsuby = nsuby
      if fig is None:
         self.fig = plt.figure(**kwargs)
         self.fig.clear()
      else:
         self.fig = fig

      self.axis = self.fig.add_subplot(111)
      self.axes = [self.axis]

   def __getattr__(self, key):
      '''Map axis attributes to this class.'''
      if hasattr(self.axis, key):
         return getattr(self.axis, key)
      else:
         raise AttributeError

   def xlabel(self, label, **args):
      return self.axis.set_xlabel(label, **args)

   def ylabel(self, label, **args):
      return self.axis.set_ylabel(label, **args)

   def title(self, label, **args):
      return self.axis.set_title(label, **args)

   def set_minor_ticks(self):
      '''Call this to setup minor ticks, if so requested.'''
      x_major = self.axis.xaxis.get_majorticklocs()
      xmajor_int = x_major[1] - x_major[0]
      xminor_int = xmajor_int/self.nsubx
      self.axis.xaxis.set_minor_locator(ticker.MultipleLocator(xminor_int)) 

      y_major = self.axis.yaxis.get_majorticklocs()
      ymajor_int = y_major[1] - y_major[0]
      yminor_int = ymajor_int/self.nsuby
      self.axis.yaxis.set_minor_locator(ticker.MultipleLocator(yminor_int)) 

   def set_fontsize(self, abs=16, rel=None):
      '''Set the fontsize of all text components in the plot.  If abs is
      specified, use it as an absolute value (e.g., 16pt).  If rel is
      specified, use it as a relative change from the current value of
      each element (e.g., 1.1 will increase all by 10%).'''

      # first, find all the text objects:
      objs = self.axis.findobj(text.Text)
      if rel is None:
         # Use absolute value
         for o in objs:
            o.set_fontsize(abs)
      else:
         for o in objs:
            old = o.get_fontsize()
            o.set_fontsize(old*rel)

   def set_grid_linewidth(self, lw):
      '''Set the linewidth of the window lines and tick marks.'''
      for obj in self.axis.findobj(matplotlib.spines.Spine):
         obj.set_linewidth(lw)
      for obj in self.axis.get_xticklines() + self.axis.get_yticklines():
         obj.set_edgelinewidth(lw)

   def set_data_linewidth(self, lw):
      for obj in self.axis.findobj(matplotlib.lines.Line2D):
         obj.set_linewidth(lw)

   def get_xlabels(self):
      '''Get instances of the axis labels, being sure to omit any
      nonsensical labels (don't know where they're coming from.'''
      xlabs = self.axis.get_xticklabels()
      xlabs.append(self.axis.get_xaxis().offsetText)
      if self.axis.get_xlabel() != '':
         xlabs.append(self.axis.get_xaxis().get_label())
      return xlabs

   def get_ylabels(self):
      '''Get instances of the axis labels, being sure to omit the first
      and last if there is an offset.'''
      ylabs = self.axis.get_yticklabels()
      ylabs.append(self.axis.get_yaxis().offsetText)
      if self.axis.get_ylabel() != '':
         ylabs.append(self.axis.get_yaxis().get_label())
      return ylabs

   def set_limits(self, pad=0.01, dox=True, doy=True, all_equal=False):
      '''Set the limits based on what's in the axis limits.'''
      if dox:
         bbox = axis_bbox(self.axis)
         if bbox is not None:
            x0 = bbox.x0 - bbox.width*pad
            x1 = bbox.x1 + bbox.width*pad
            if self.axis.xaxis_inverted():
               self.axis.set_xlim((x1,x0))
            else:
               self.axis.set_xlim((x0,x1))
      if doy:
         bbox = axis_bbox(self.axis)
         if bbox is not None:
            y0 = bbox.y0 - bbox.height*pad
            y1 = bbox.y1 + bbox.height*pad
            if self.axis.yaxis_inverted():
               self.axis.set_ylim((y1,y0))
            else:
               self.axis.set_ylim((y0,y1))

   def draw(self, hide_corner_labels=1):
      '''Draw the panel and everything in it.'''
      self.set_minor_ticks()
      plt.draw()
      # Now that everything's been rendered, let's clean up shop:
      # Make room for everything:
      self.fig.tight_layout()
      plt.draw()

   def close(self):
      plt.close(self.fig)

class MultiPlot(object):
   '''A multi-plot with NXM panels each containing a SimplePlot.'''

   def __init__(self, N, M, fig=None, pwidths=None, pheights=None, 
         nsubx=5, nsuby=5, **kwargs):
      if fig is None:
         self.fig = plt.figure(**kwargs)
         self.fig.clear()
      else:
         self.fig = fig
      self.fig.clear()
      self.N = N
      self.M = M
      self.num = N*M
      
      self.nsubx = 5
      self.nsuby = 5

      # Note gridspec counts grids from the *top*. It also gives dimensions
      #  as (columns X rows), whereas myplotlib is (rows X columns)
      self.gs = gridspec.GridSpec(self.M, self.N, width_ratios=pwidths,
            height_ratios=pheights[::-1])
      self.axes = [self.fig.add_subplot(g) for g in self.gs]
      self._xlabel = None
      self._ylabel = None
      self._title = None

   def __getitem__(self, idx):
      # Implement self[idx], where idx can be an integer or 2-tuple of
      # integers. Negative indices are supported. Slices are not
      if type(idx) is int:
         # interpreted as the plain axis object
         return self.axes[idx]
      elif type(idx) is tuple:
         if len(idx) == 2:
            if type(idx[0]) is int and type(idx[1]) is int:
               i,j = idx
               if i < 0: i += self.N
               if j < 0: j += self.M
               return self.axes[self.idx(i,j)]
      raise TypeError, "index must be an integer or 2-tuple of integers"


   def title(self, string, **kws):
      '''Set a title for the Multi plots.  kws can be any
      arguments recognized by figure.text()'''

      t_ax = self.N/2    # always right, even when odd/even
      idx = self.idx(t_ax,self.M-1)

      # This will give us space using tight_layout()
      self._title = self.axes[idx].set_title(string, **kws)
      return self._title

   def xlabel(self, string, **kws):
      '''Set an x-label for the Panel plots.  kws can be any
      arguments recognized by figure.text()'''

      l_ax = self.N/2    # always right, even when odd/even
      self._xlabel = self.axes[l_ax].set_xlabel(string, **kws)
      return self._xlabel

   def ylabel(self, string, **kws):
      '''Set an x-label for the Panel plots.  kws can be any
      arguments recognized by figure.text()'''

      j_ax = self.M/2    # j coordinate of left-hand middle row
      self._ylabel = self.axes[self.idx(0,j_ax)].set_ylabel(string, **kws)
      return self._ylabel

   def set_minor_ticks(self):
      '''Call this to setup minor ticks, if so requested.'''
      for ax in self.axes:
         x_major = ax.xaxis.get_majorticklocs()
         xmajor_int = x_major[1] - x_major[0]
         xminor_int = xmajor_int/self.nsubx
         ax.xaxis.set_minor_locator(ticker.MultipleLocator(xminor_int)) 
      
         y_major = ax.yaxis.get_majorticklocs()
         ymajor_int = y_major[1] - y_major[0]
         yminor_int = ymajor_int/self.nsuby
         ax.yaxis.set_minor_locator(ticker.MultipleLocator(yminor_int)) 

   def set_fontsize(self, abs=16, rel=None):
      '''Set the fontsize of all text components in the plot.  If abs is
      specified, use it as an absolute value (e.g., 16pt).  If rel is
      specified, use it as a relative change from the current value of
      each element (e.g., 1.1 will increase all by 10%).'''

      # first, find all the text objects:
      objs = self.fig.findobj(text.Text)
      if rel is None:
         # Use absolute value
         for o in objs:
            o.set_fontsize(abs)
      else:
         for o in objs:
            old = o.get_fontsize()
            o.set_fontsize(old*rel)

   def set_grid_linewidth(self, lw):
      '''Set the linewidth of the window lines and tick marks.'''
      for obj in self.fig.findobj(matplotlib.spines.Spine):
         obj.set_linewidth(lw)
      for p in self.axes:
         for obj in p.get_xticklines() + p.get_yticklines():
            obj.set_edgelinewidth(lw)

   def set_data_linewidth(self, lw):
      for obj in self.fig.findobj(matplotlib.lines.Line2D):
         obj.set_linewidth(lw)

   def ij(self, i):
      return (i%self.N, i/self.N)

   def idx(self, i, j):
      return j*self.N + i

   def set_limits(self, pad=0.01, dox=True, doy=True, all_equal=False):
      '''Go through the rows and columns and set the x and y limits
      to fit the data.'''
      if all_equal:
         # Set all panel to have the same bounds
         bboxs = [axis_bbox(ax) \
               for ax in self.axes if axis_bbox(ax) is not None]
         bbox = transforms.Bbox.union(bboxs)
         x0 = bbox.x0 - bbox.width*pad
         x1 = bbox.x1 + bbox.width*pad
         y0 = bbox.y0 - bbox.height*pad
         y1 = bbox.y1 + bbox.height*pad
         if dox:
            for i in range(self.N):
               if self.axes[i].xaxis_inverted():
                  self.axes[i].set_xlim((x1,x0))
               else:
                  self.axes[i].set_xlim((x0,x1))
         if doy:
            for j in range(self.M):
               if self.axes[j*self.N].yaxis_inverted():
                  self.axes[j*self.N].set_ylim((y1,y0))
               else:
                  self.axes[j*self.N].set_ylim((y0,y1))
      else:
         if dox:
            for i in range(self.num):
               bbox = axis_bbox(self.axes[i])
               if bbox is not None:
                  x0 = bbox.x0 - bbox.width*pad
                  x1 = bbox.x1 + bbox.width*pad
                  if self.axes[i].xaxis_inverted():
                     self.axes[i].set_xlim((x1,x0))
                  else:
                     self.axes[i].set_xlim((x0,x1))
         if doy:
            for j in range(self.num):
               bbox = axis_bbox(self.axes[j])
               if bbox is not None:
                  y0 = bbox.y0 - bbox.height*pad
                  y1 = bbox.y1 + bbox.height*pad
                  if self.axes[j].yaxis_inverted():
                     self.axes[j].set_ylim((y1,y0))
                  else:
                     self.axes[j].set_ylim((y0,y1))

   def get_xticklabels(self, k):
      '''Get instances of the axis labels, being sure to omit any
      nonsensical labels (don't know where they're coming from.'''
      xlabs = self.axes[k].get_xticklabels()
      xlabs.append(self.axes[k].get_xaxis().offsetText)
      if self.axes[k].get_xlabel() != '':
         xlabs.append(self.axes[k].get_xaxis().get_label())
      fbbox = self.fig.patch.get_window_extent(self.get_renderer())
      xlabs = [lab for lab in xlabs \
            if fbbox.overlaps(lab.get_window_extent(self.get_renderer()))]
      return xlabs

   def get_yticklabels(self, k):
      '''Get instances of the axis labels, being sure to omit the first
      and last if there is an offset.'''
      ylabs = self.axes[k].get_yticklabels()
      ylabs.append(self.axes[k].get_yaxis().offsetText)
      if self.axes[k].get_ylabel() != '':
         ylabs.append(self.axes[k].get_yaxis().get_label())
      fbbox = self.fig.patch.get_window_extent(self.get_renderer())
      ylabs = [lab for lab in ylabs \
            if fbbox.overlaps(lab.get_window_extent(self.get_renderer()))]
      return ylabs

   def draw(self, hide_corner_labels=1):
      '''Draw the panel and everything in it.'''
      plt.draw()
      # Now that everything's been rendered, let's clean up shop:
      self.fig.tight_layout()

      # Now, if our figure-wide xlabel, ylable or title are defined,
      # we may need to nudge them over half a panel
      if self._title is not None and not self.N % 2:
         x,y = self._title.get_position()
         self._title.set_position((0.0, y))
      if self._xlabel is not None and not self.N % 2:
         x,y = self._xlabel.get_position()
         self._xlabel.set_position((0.0, y))
      if self._ylabel is not None and not self.M % 2:
         x,y = self._ylabel.get_position()
         self._ylabel.set_position((x, 0.0))
      plt.draw()

   def close(self):
      plt.close(self.fig)

class PanelPlot(MultiPlot):

   def __init__(self, N, M, fig=None, pwidths=None, pheights=None, 
         nsubx=5, nsuby=5, **kwargs):
      if fig is None:
         self.fig = plt.figure(**kwargs)
         self.fig.clear()
      else:
         self.fig = fig
      self.fig.clear()
      self.N = N
      self.M = M
      self.num = N*M

      self.nsubx = nsubx
      self.nsuby = nsuby

      # Note gridspec counts grids from the *top*. It also gives dimensions
      #  as (columns X rows), whereas myplotlib is (rows X columns)
      if pheights is not None:
         pheights = pheights[::-1]
      self.gs = gridspec.GridSpec(self.M, self.N, width_ratios=pwidths,
            height_ratios=pheights)
      self.gs.update(wspace=0, hspace=0)
      #self.axes = [self.fig.add_subplot(g) for g in self.gs]
      self.axes = []
      for j in range(self.M):
         for i in range(self.N):
            if i > 0:
               sharey = self.axes[-1]
            else:
               sharey = None
            if j > 0:
               sharex = self.axes[-self.N]
            else:
               sharex = None
            jj = i
            ii = self.M-j-1
            self.axes.append(self.fig.add_subplot(self.gs[ii,jj], sharex=sharex,
               sharey=sharey))
            if i > 0:
               self.axes[-1].tick_params(axis='y', labelleft=False)
            if j > 0:
               self.axes[-1].tick_params(axis='x', labelbottom=False)

      self._xlabel = None
      self._ylabel = None
      self._title = None

   def set_limits(self, pad=0.01, dox=True, doy=True, all_equal=False):
      # Here we have the complication that some axes are shared. And
      # so we need to consider the x-range and y-range over rows and
      # columns.
      if all_equal:
         # Set all panel to have the same bounds
         bboxs = [axis_bbox(ax) for ax in self.axes \
                                if axis_bbox(ax) is not None]
         bbox = transforms.Bbox.union(bboxs)
         x0 = bbox.x0 - bbox.width*pad
         x1 = bbox.x1 + bbox.width*pad
         y0 = bbox.y0 - bbox.height*pad
         y1 = bbox.y1 + bbox.height*pad
         if dox:
            for i in range(self.N):
               if self.axes[i].xaxis_inverted():
                  self.axes[i].set_xlim((x1,x0))
               else:
                  self.axes[i].set_xlim((x0,x1))
         if doy:
            for j in range(self.M):
               if self.axes[j*self.N].yaxis_inverted():
                  self.axes[j*self.N].set_ylim((y1,y0))
               else:
                  self.axes[j*self.N].set_ylim((y0,y1))
      else:
         if dox:
            for i in range(self.N):
               bboxs = [axis_bbox(self.axes[i+j*self.N]) for j in range(self.M)\
                     if axis_bbox(self.axes[i+j*self.N]) is not None]
               bbox = transforms.Bbox.union(bboxs)
               x0 = bbox.x0 - bbox.width*pad
               x1 = bbox.x1 + bbox.width*pad
               if self.axes[i].xaxis_inverted():
                  self.axes[i].set_xlim((x1,x0))
               else:
                  self.axes[i].set_xlim((x0,x1))
         if doy:
            for j in range(self.M):
               bboxs = [axis_bbox(self.axes[i+j*self.N]) for i in range(self.N)\
                     if axis_bbox(self.axes[i+j*self.N]) is not None]
               bbox = transforms.Bbox.union(bboxs)
               y0 = bbox.y0 - bbox.height*pad
               y1 = bbox.y1 + bbox.height*pad
               if self.axes[j*self.N].yaxis_inverted():
                  self.axes[j*self.N].set_ylim((y1,y0))
               else:
                  self.axes[j*self.N].set_ylim((y0,y1))

   def draw(self, hide_corner_labels=1):
      '''Draw the panel and everything in it.'''
      plt.draw()
      # Now that everything's been rendered, let's clean up shop:
      self.fig.tight_layout()

      # Now, if our figure-wide xlabel, ylable or title are defined,
      # we may need to nudge them over half a panel
      if self._title is not None and not self.N % 2:
         x,y = self._title.get_position()
         self._title.set_position((0.0, y))
      if self._xlabel is not None and not self.N % 2:
         x,y = self._xlabel.get_position()
         self._xlabel.set_position((0.0, y))
      if self._ylabel is not None and not self.M % 2:
         x,y = self._ylabel.get_position()
         self._ylabel.set_position((x, 0.0))
      plt.draw()

