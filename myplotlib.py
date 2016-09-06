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
from matplotlib.ticker import MaxNLocator,NullFormatter
from numpy import *
from matplotlib import ticker

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

def get_renderer(fig):
    if fig._cachedRenderer:
        renderer = fig._cachedRenderer
    else:
        canvas = fig.canvas

        if canvas and hasattr(canvas, "get_renderer"):
            renderer = canvas.get_renderer()
        else:
            # not sure if this can happen
            #warnings.warn("tight_layout : falling back to Agg renderer")
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            canvas = FigureCanvasAgg(fig)
            renderer = canvas.get_renderer()

    return renderer

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
   min_pad = 0.05         # a minimum padding aroudn the axis (so that eve
                          #  with no labels, we still see the grid

   def __init__(self, fig=None, position=None, nsubx=5, nsuby=5, **kwargs):
      self.left_pad = self.min_pad
      self.right_pad = self.min_pad
      self.bottom_pad = self.min_pad
      self.top_pad = self.min_pad

      self.nsubx = nsubx
      self.nsuby = nsuby
      if position is None:
         self.position = (0,0,1,1)    # Default:  take up whole figure
      else:
         if type(position) != type(()):
            raise ValueError, "position must be a 4-tuple"
         self.position = position
      if fig is None:
         self.fig = plt.figure(**kwargs)
         self.fig.clear()
      else:
         self.fig = fig

      pos = self.get_axis_position()
      self.axis = self.fig.add_axes(pos, autoscale_on=False)
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

   def get_axis_position(self):
      '''Based on the current plot position and padding, compute the axis
      position.'''
      x0 = self.position[0] + self.left_pad
      y0 = self.position[1] + self.bottom_pad
      width = self.position[2] - self.left_pad - self.right_pad
      height = self.position[3] - self.top_pad - self.bottom_pad
      return((x0, y0, width, height))

   def set_position(self, position):
      self.position = position
      self.reset_panel_positions()

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

   def get_renderer(self):
      return get_renderer(self.fig)

   def get_bbox(self, label):
      '''Get the bounding box of an object in Figure coordinates'''
      bbox = label.get_window_extent(self.get_renderer())
      bboxi = bbox.inverse_transformed(self.fig.transFigure)
      return bboxi

   def get_xlabels(self):
      '''Get instances of the axis labels, being sure to omit any
      nonsensical labels (don't know where they're coming from.'''
      xlabs = self.axis.get_xticklabels()
      xlabs.append(self.axis.get_xaxis().offsetText)
      if self.axis.get_xlabel() != '':
         xlabs.append(self.axis.get_xaxis().get_label())
      fbbox = self.fig.patch.get_window_extent(self.get_renderer())
      #xlabs = [lab for lab in xlabs if fbbox.overlaps(lab.get_window_extent())]
      return xlabs

   def get_ylabels(self):
      '''Get instances of the axis labels, being sure to omit the first
      and last if there is an offset.'''
      ylabs = self.axis.get_yticklabels()
      ylabs.append(self.axis.get_yaxis().offsetText)
      if self.axis.get_ylabel() != '':
         ylabs.append(self.axis.get_yaxis().get_label())
      fbbox = self.fig.patch.get_window_extent(self.get_renderer())
      #ylabs = [lab for lab in ylabs if fbbox.overlaps(lab.get_window_extent())]
      return ylabs

   def get_ylabels_bbox(self):
      '''Get the bounding box of the ylabels'''
      labels = self.get_ylabels()

      bboxes = []
      for label in labels:
         if label.get_visible():
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.inverse_transformed(self.fig.transFigure)
            bboxes.append(bboxi)
      # Now put them all together into self.N uber-boxes
      bbox = transforms.Bbox.union(bboxes)
      return bbox

   def get_xlabels_bbox(self):
      '''Get the bounding box for all the x-axes in the grid'''
      labels = self.get_xlabels()

      bboxes = []
      for label in labels:
         if label.get_visible():
            #print label
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.inverse_transformed(self.fig.transFigure)
            bboxes.append(bboxi)
      # Now put them all together into self.M uber-bboxes:
      bbox = transforms.Bbox.union(bboxes)
      return bbox

   def reset_panel_positions(self):
      '''Change the axis position to match the current paddings and 
      position on the figure.'''
      pos = self.get_axis_position()
      self.axis.set_position(pos)

   def set_left_padding(self, dx):
      '''Add extra padding on the left and re-scale everything.'''
      self.left_pad = dx + self.min_pad
      self.reset_panel_positions()

   def add_left_padding(self, dx):
      self.left_pad += dx
      self.reset_panel_positions()

   def set_right_padding(self, dx):
      '''Add extra padding on the right and re-scale everything.'''
      self.right_pad = dx + self.min_pad
      self.reset_panel_positions()

   def add_right_padding(self, dx):
      self.right_pad += dx
      self.reset_panel_positions()

   def set_bottom_padding(self, dy):
      '''Add extra padding on the bottom and re-scale everything.'''
      self.bottom_pad = dy + self.min_pad
      self.reset_panel_positions()

   def add_bottom_padding(self, dy):
      self.bottom_pad += dy
      self.reset_panel_positions()

   def set_top_padding(self, dy):
      '''Add extra padding on the bottom and re-scale everything.'''
      self.top_pad = dy + self.min_pad
      self.reset_panel_positions()

   def add_top_padding(self, dy):
      self.top_padd += dy
      self.reset_panel_positions()

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
      #self.fig.canvas.draw()
      plt.draw()
      # Now that everything's been rendered, let's clean up shop:
      # Make room for everything:

      if self.axis.get_title != "":
         title_height = self.get_bbox(self.axis.title).height
      else:
         title_height = 0
      axis_label_height = self.get_xlabels_bbox().height
      axis_label_width = self.get_ylabels_bbox().width

      top_height = title_height
      bottom_height = axis_label_height
      left_width = axis_label_width

      self.set_top_padding(top_height)
      self.set_bottom_padding(bottom_height)
      self.set_left_padding(left_width)

      plt.draw()

   def close(self):
      plt.close(self.fig)



class MultiPlot:
   '''A multi-plot with NXM panels each containing a SimplePlot.'''
   min_pad = 0.05

   def __init__(self, N, M, fig=None, pwidths=None, pheights=None, **kwargs):
      self.left_pad = self.min_pad
      self.right_pad = self.min_pad
      self.bottom_pad = self.min_pad
      self.top_pad = self.min_pad

      if fig is None:
         self.fig = plt.figure(**kwargs)
         self.fig.clear()
      else:
         self.fig = fig
      self.fig.clear()
      self.N = N
      self.M = M
      self.num = N*M
      self._title = None
      self._xlabel = None
      self._ylabel = None
      if pwidths is None:
         self.rel_widths = array([1.0/N]*N)
      else:
         if len(pwidths) != self.N:
            raise AttributeError, "lenth of pwidts must equal N"
         self.rel_widths = array(pwidths)*1.0/sum(pwidths)
      self.pwidths = (1.0 - self.left_pad - self.right_pad)*self.rel_widths
      if pheights is None:
         self.rel_heights = array([1.0/M]*M)
      else:
         if len(pheights) != self.M:
            raise AttributeError, "length of pheights must equal M"
         self.rel_heights = array(pheights)*1.0/sum(pheights)
      self.pheights = (1.0 - self.bottom_pad - self.top_pad)*self.rel_heights

      # Set some initial panels with default padding
      self.axes = []
      for k in range(self.num):
         i,j = self.ij(k)
         x0 = self.left_pad + sum(self.pwidths[0:i])
         y0 = self.bottom_pad + sum(self.pheights[0:j])
         rect = (x0,y0,self.pwidths[i],self.pheights[j])
         self.axes.append(SimplePlot(self.fig, rect))

   def title(self, string, **kws):
      '''Set a title for the Multi plots.  kws can be any
      arguments recognized by figure.text()'''

      if 'horizontalalignment' not in kws and 'ha' not in kws:
         kws['horizontalalignment'] = 'center'
      if 'verticalalignment' not in kws or 'va' not in kws:
         kws['verticalalignment'] = 'top'
      #x = self.left_pad + (1.0 - self.left_pad - self.right_pad)/2
      x = 0.5
      self._title = self.fig.text(x, 0.95, string, **kws)
      return self._title

   def xlabel(self, string, **kws):
      '''Set an x-label for the Panel plots.  kws can be any
      arguments recognized by figure.text()'''

      if 'horizontalalignment' not in kws and 'ha' not in kws:
         kws['horizontalalignment'] = 'center'
      if 'verticalalignment' not in kws and 'va' not in kws:
         kws['verticalalignment'] = 'bottom'
      #x = self.left_pad + (1.0 - self.left_pad - self.right_pad)/2
      x = 0.5
      self._xlabel = self.fig.text(x, 0.05, string, **kws)
      return self._xlabel

   def ylabel(self, string, **kws):
      '''Set an x-label for the Panel plots.  kws can be any
      arguments recognized by figure.text()'''

      if 'horizontalalignment' not in kws and 'ha' not in kws:
         kws['horizontalalignment'] = 'right'
      if 'verticalalignment' not in kws and 'va' not in kws:
         kws['verticalalignment'] = 'center'
      if 'rotation' not in kws:
         kws['rotation'] = 90
      #y = self.bottom_pad + (1.0 - self.bottom_pad - self.top_pad)/2
      y = 0.5
      self._ylabel = self.fig.text(0.05, y, string, **kws)
      return self._ylabel

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

   def get_renderer(self):
      return get_renderer(self.fig)

   def get_bbox(self, label):
      '''If a title exists, how tall is it?'''
      bbox = label.get_window_extent(self.get_renderer())
      bboxi = bbox.inverse_transformed(self.fig.transFigure)
      return bboxi

   def reset_panel_positions(self):
      '''Change the axes positions to match the current paddings.
      Also, update the positions of the Panels' labels if they are
      defined.'''
      for k in range(self.num):
         i,j = self.ij(k)
         x0 = self.left_pad + sum(self.pwidths[0:i])
         y0 = self.bottom_pad + sum(self.pheights[0:j])
         rect = [x0, y0, self.pwidths[i], self.pheights[j]]
         self.axes[k].set_position(rect)

      if self._xlabel is not None:
         #x = self.left_pad + (1 - self.left_pad - self.right_pad)/2
         x = 0.5
         y = self.bottom_pad*0.25
         self._xlabel.set_position((x,y))

      if self._ylabel is not None:
         #y = self.bottom_pad + (1 - self.bottom_pad - self.top_pad)/2
         y = 0.5
         x = self.left_pad*0.25
         self._ylabel.set_position((x,y))

      if self._title is not None:
         #x = self.left_pad + (1 - self.left_pad - self.right_pad)/2
         x = 0.5
         y = 1 - self.top_pad*0.25
         self._title.set_position((x,y))

   def reset_pwidths_pheights(self):
      self.pwidths = (1.0 - self.left_pad - self.right_pad)*self.rel_widths
      self.pheights = (1.0 - self.bottom_pad - self.top_pad)*self.rel_heights

   def set_left_padding(self, dx):
      '''Add extra padding on the left and re-scale everything.'''
      self.left_pad = dx + self.min_pad
      self.reset_pwidths_pheights()
      self.reset_panel_positions()

   def set_right_padding(self, dx):
      '''Add extra padding on the right and re-scale everything.'''
      self.right_pad = dx + self.min_pad
      self.reset_pwidths_pheights()
      self.reset_panel_positions()

   def set_bottom_padding(self, dy):
      '''Add extra padding on the bottom and re-scale everything.'''
      self.bottom_pad = dy + self.min_pad
      self.reset_pwidths_pheights()
      self.reset_panel_positions()

   def set_top_padding(self, dy):
      '''Add extra padding on the bottom and re-scale everything.'''
      self.top_pad = dy + self.min_pad
      self.reset_pwidths_pheights()
      self.reset_panel_positions()

   def set_limits(self, pad=0.01, dox=True, doy=True, all_equal=False):
      '''Go through the rows and columns and set the x and y limits
      to fit the data.'''
      if all_equal:
         # Set all panel to have the same bounds
         bboxs = [axis_bbox(ax) for ax in self.axes if axis_bbox(ax) is not None]
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
      #self.fig.canvas.draw()
      plt.draw()
      # Now that everything's been rendered, let's clean up shop:
      # Make room for everything:
      if self._title is not None:
         title_height = self.get_bbox(self._title).height
      else:
         title_height = 0
      if self._xlabel is not None:
         xlabel_height = self.get_bbox(self._xlabel).height
      else:
         xlabel_height = 0
      if self._ylabel is not None:
         ylabel_width = self.get_bbox(self._ylabel).width
      else:
         ylabel_width = 0

      top_height = title_height*1.0
      bottom_height = xlabel_height
      left_width = ylabel_width

      self.set_top_padding(top_height)
      self.set_bottom_padding(bottom_height)
      self.set_left_padding(left_width)

      plt.draw()

   def close(self):
      plt.close(self.fig)

class PanelPlot:
   '''A multi-plot with NXM panels, all squeezed together.'''

   #left_pad = 0.05
   #right_pad = 0.05
   #bottom_pad = 0.05
   #top_pad = 0.05
   min_pad = 0.05

   def __init__(self, N, M, fig=None, pwidths=None, pheights=None, xlock=1,
         ylock=1, nxmax=None, nymax=None, nsubx=5, nsuby=5, prunex='both',
         pruney='both', left_pad=0, right_pad=0, bottom_pad=0,
         top_pad=0, **kwargs):
      '''Create a panelplot, where each panel is snuggled up to its 
      neighbors and they effectively share interior axes.
      Arguments:
         N,M      :  make an NXM panel plots (N = number of columns,
                     M = number of rows)
         fig      :  Supply a figure object, othewise, one is created.
         pwidths  :  specify list of N widths (in device coordinates).  Default
                     is to make them equal.
         pheights :  specify list of M heights (in device coordinates). Default
                     is to make them equal.
         [xy]lock :  lock the x ranges and y ranges to for them to be equal
                     on shared internal axes
         n[xy]max:   maximum number of major ticks (and labels).  Good if
                     things get crowded
         nsub[xy]:   number of minor ticks to place between major ticks
         prune[xy]:  ('lower','upper','both')  prune the lower/upper/both tick
                     labels in the panels.  That way, they don't overlap at
                     the boundaries.
         
         Any other arguments are passed to figure()'''

      if fig is None:
         self.fig = plt.figure(**kwargs)
         self.fig.clear()
      else:
         self.fig = fig
      self.fig.clear()
      self.N = N
      self.M = M
      self.left_pad = self.min_pad
      self.right_pad = self.min_pad
      self.bottom_pad = self.min_pad
      self.top_pad = self.min_pad
      self.extra_left_padding = left_pad
      self.extra_right_padding = right_pad
      self.extra_top_padding = top_pad
      self.extra_bottom_padding = bottom_pad
      self.nxmax = nxmax
      self.nymax = nymax
      self.nsubx = nsubx
      self.nsuby = nsuby
      self.num = N*M
      self.axes = []
      self._title = None
      self._xlabel = None
      self._ylabel = None
      if pwidths is None:
         self.rel_widths = array([1.0/N]*N)
      else:
         if len(pwidths) != self.N:
            raise AttributeError, "lenth of pwidts must equal N"
         self.rel_widths = array(pwidths)*1.0/sum(pwidths)
      self.pwidths = (1.0 - self.left_pad - self.right_pad)*self.rel_widths
      if pheights is None:
         self.rel_heights = array([1.0/M]*M)
      else:
         if len(pheights) != self.M:
            raise AttributeError, "length of pheights must equal M"
         self.rel_heights = array(pheights)*1.0/sum(pheights)
      self.pheights = (1.0 - self.bottom_pad - self.top_pad)*self.rel_heights

      # Set some initial panels with default padding
      for k in range(self.num):
         i,j = self.ij(k)
         x0 = self.left_pad + sum(self.pwidths[0:i])
         y0 = self.bottom_pad + sum(self.pheights[0:j])
         rect = [x0,y0,self.pwidths[i],self.pheights[j]]
         opts = {}
         if i != 0 and ylock:
            # not on the left-hand edge
            opts['sharey'] = self.axes[j*N]
         if j != 0 and xlock:
            # not on the lower edge
            opts['sharex'] = self.axes[i]
         self.axes.append(self.fig.add_axes(rect, autoscale_on=False,
            **opts))
         if self.nxmax is not None:
            self.axes[-1].xaxis.set_major_locator(
                  MaxNLocator(self.nxmax, prune=prunex))
         if self.nymax is not None:
            self.axes[-1].yaxis.set_major_locator(
                  MaxNLocator(self.nymax, prune=pruney))


   def title(self, string, **kws):
      '''Set a title for the Panel plots.  kws can be any
      arguments recognized by figure.text()'''

      if 'horizontalalignment' not in kws:
         kws['horizontalalignment'] = 'center'
      if 'verticalalignment' not in kws:
         kws['verticalalignment'] = 'top'
      x = self.left_pad + (1.0 - self.left_pad - self.right_pad)/2
      self._title = self.fig.text(x, 0.95, string, **kws)
      return self._title

   def xlabel(self, string, **kws):
      '''Set an x-label for the Panel plots.  kws can be any
      arguments recognized by figure.text()'''

      if 'horizontalalignment' not in kws and 'ha' not in kws:
         kws['horizontalalignment'] = 'center'
      if 'verticalalignment' not in kws and 'va' not in kws:
         kws['verticalalignment'] = 'top'
      x = self.left_pad + (1.0 - self.left_pad - self.right_pad)/2
      self._xlabel = self.fig.text(x, 0.05, string, **kws)
      return self._xlabel

   def ylabel(self, string, labelpad=0,**kws):
      '''Set a y-label for the Panel plots.  kws can be any
      arguments recognized by figure.text()'''

      if 'horizontalalignment' not in kws and 'ha' not in kws:
         kws['horizontalalignment'] = 'right'
      if 'verticalalignment' not in kws and 'va' not in kws:
         kws['verticalalignment'] = 'center'
      if 'rotation' not in kws:
         kws['rotation'] = 90
      y = self.bottom_pad + (1.0 - self.bottom_pad - self.top_pad)/2
      self._ylabel = self.fig.text(self.left_pad - labelpad, y, string, **kws)
      return self._ylabel

   def ij(self, i):
      return (i%self.N, i/self.N)

   def set_minor_ticks(self):
      '''Call this to setup minor ticks, if so requested.'''
      for ax in self.axes:
         ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(self.nsubx))
         #x_major = ax.xaxis.get_majorticklocs()
         #xmajor_int = x_major[1] - x_major[0]
         #xminor_int = xmajor_int/self.nsubx
         #ax.xaxis.set_minor_locator(ticker.MultipleLocator(xminor_int)) 
       
         ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(self.nsuby))
         #y_major = ax.yaxis.get_majorticklocs()
         #ymajor_int = y_major[1] - y_major[0]
         #yminor_int = ymajor_int/self.nsuby
         #ax.yaxis.set_minor_locator(ticker.MultipleLocator(yminor_int)) 

   def get_renderer(self):
      return get_renderer(self.fig)

   def get_bbox(self, label):
      '''If a title exists, how tall is it?'''
      bbox = label.get_window_extent(self.get_renderer())
      bboxi = bbox.inverse_transformed(self.fig.transFigure)
      return bboxi

   def get_ylabels_bbox(self):
      '''Get the bounding box for all the y-axes in the grid'''
      labels = []
      for j in range(self.M):
         labels += self.get_yticklabels(j*self.N)
         if self.axes[j*self.N].yaxis.get_label().get_text() != '': 
            labels.append(self.axes[j*self.N].yaxis.get_label())
      bboxes = []
      for label in labels:
         if label.get_visible():
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.inverse_transformed(self.fig.transFigure)
            bboxes.append(bboxi)
      # Now put them all together into an uber-bbox:
      bbox = transforms.Bbox.union(bboxes)
      return bbox

   def get_xlabels_bbox(self):
      '''Get the bounding box for all the x-axes in the grid'''
      labels = []
      for i in range(self.N):
         labels += self.get_xticklabels(i)
         if self.axes[i].xaxis.get_label().get_text() != '': 
            labels.append(self.axes[i].xaxis.get_label())
      bboxes = []
      for label in labels:
         if label.get_visible():
            bbox = label.get_window_extent(self.get_renderer())
            bboxi = bbox.inverse_transformed(self.fig.transFigure)
            bboxes.append(bboxi)
      # Now put them all together into an uber-bbox:
      bbox = transforms.Bbox.union(bboxes)
      return bbox

   def reset_panel_positions(self):
      '''Change the axes positions to match the current paddings.
      Also, update the positions of the Panels' labels if they are
      defined.'''
      for k in range(self.num):
         i,j = self.ij(k)
         rect = [self.left_pad + sum(self.pwidths[0:i]),
               self.bottom_pad + sum(self.pheights[0:j]),
                 self.pwidths[i], self.pheights[j]]
         self.axes[k].set_position(rect)
      if self._xlabel is not None:
         x = self.left_pad + (1 - self.left_pad - self.right_pad)/2
         y = self.bottom_pad*0.4
         self._xlabel.set_position((x,y))

      if self._ylabel is not None:
         y = self.bottom_pad + (1 - self.bottom_pad - self.top_pad)/2
         x = self.left_pad*0.4
         self._ylabel.set_position((x,y))

      if self._title is not None:
         x = self.left_pad + (1 - self.left_pad - self.right_pad)/2
         y = 1 - self.top_pad*0.4
         self._title.set_position((x,y))

   def set_left_padding(self, dx):
      '''Add extra padding on the left and re-scale everything.'''
      self.left_pad = dx + self.min_pad + self.extra_left_padding
      self.pwidths = (1.0 - self.left_pad - self.right_pad)*self.rel_widths
      self.reset_panel_positions()

   def set_right_padding(self, dx):
      '''Add extra padding on the right and re-scale everything.'''
      self.right_pad = dx + self.min_pad + self.extra_right_padding
      self.pwidths = (1.0 - self.left_pad - self.right_pad)*self.rel_widths
      self.reset_panel_positions()

   def set_bottom_padding(self, dy):
      '''Add extra padding on the bottom and re-scale everything.'''
      self.bottom_pad = dy + self.min_pad + self.extra_bottom_padding
      self.pheights = (1.0 - self.bottom_pad - self.top_pad)*self.rel_heights
      self.reset_panel_positions()

   def set_top_padding(self, dy):
      '''Add extra padding on the bottom and re-scale everything.'''
      self.top_pad = dy + self.min_pad + self.extra_top_padding
      self.pheights = (1.0 - self.bottom_pad - self.top_pad)*self.rel_heights
      self.reset_panel_positions()

   def set_limits(self, pad=0.01, dox=True, doy=True, all_equal=False):
      '''Go through the rows and columns and set the x and y limits
      to fit the data.'''
      if all_equal:
         # Set all panel to have the same bounds
         bboxs = [axis_bbox(ax) for ax in self.axes if axis_bbox(ax) is not None]
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

   def get_xticklabels(self, k):
      '''Get instances of the axis labels, being sure to omit any
      nonsensical labels (don't know where they're coming from.'''
      xlabs = self.axes[k].get_xticklabels()
      xlabs.append(self.axes[k].get_xaxis().offsetText)
      fbbox = self.fig.patch.get_window_extent(self.get_renderer())
      xlabs = [lab for lab in xlabs \
            if fbbox.overlaps(lab.get_window_extent(self.get_renderer()))]
      return xlabs

   def get_yticklabels(self, k):
      '''Get instances of the axis labels, being sure to omit the first
      and last if there is an offset.'''
      ylabs = self.axes[k].get_yticklabels()
      ylabs.append(self.axes[k].get_yaxis().offsetText)
      fbbox = self.fig.patch.get_window_extent(self.get_renderer())
      ylabs = [lab for lab in ylabs \
            if fbbox.overlaps(lab.get_window_extent(self.get_renderer()))]
      return ylabs

   def draw(self, hide_corner_labels=1):
      '''Draw the panel and everything in it.'''
      #self.fig.canvas.draw()
      plt.draw()
      self.set_minor_ticks()
      # Now that everything's been rendered, let's clean up shop:
      # 1) get rid of interiour tick labels
      for k in range(self.num):
         i,j = self.ij(k)
         if i != 0:
            # Not on the left edge
            for lab in self.axes[k].get_yticklabels():
               lab.set_visible(False)
            self.axes[k].yaxis.offsetText.set_visible(False)
         if j != 0:
            # Not on the bottom edge
            for lab in self.axes[k].get_xticklabels():
               lab.set_visible(False)
            self.axes[k].xaxis.offsetText.set_visible(False)
      #2) Get rid of any corner tick labels that overlap:
      # first, make them all visible
      #if hide_corner_labels:
      #   for i in range(1,self.N):
      #      lab1s = self.get_xticklabels(i-1)[0:-1]
      #      lab2s = self.get_xticklabels(i)[0:-1]
      #      [lab.set_visible(True) for lab in lab1s]
      #      [lab.set_visible(True) for lab in lab2s]
      #   for j in range(1,self.M):
      #      lab1s = self.get_yticklabels((j-1)*self.N)[0:-1]
      #      lab2s = self.get_yticklabels(j*self.N)[0:-1]
      #      [lab.set_visible(True) for lab in lab1s]
      #      [lab.set_visible(True) for lab in lab2s]
      #   plt.draw()
      #   for i in range(1,self.N):
      #      lab1s = self.get_xticklabels(i-1)[0:-1]
      #      bboxes1 = [lab.get_window_extent() for lab in lab1s]
      #      lab2s = self.get_xticklabels(i)[0:-1]
      #      bboxes2 = [lab.get_window_extent() for lab in lab2s]
      #      for j in range(len(bboxes1)):
      #         for k in range(len(bboxes2)):
      #            if bboxes1[j].overlaps(bboxes2[k]):
      #               lab2s[k].set_visible(False)

      #   for j in range(1,self.M):
      #      lab1s = self.get_yticklabels((j-1)*self.N)[0:-1]
      #      bboxes1 = [lab.get_window_extent() for lab in lab1s]
      #      lab2s = self.get_yticklabels(j*self.N)[0:-1]
      #      bboxes2 = [lab.get_window_extent() for lab in lab2s]
      #      for i in range(len(bboxes1)):
      #         for k in range(len(bboxes2)):
      #            if bboxes1[i].overlaps(bboxes2[k]):
      #               lab2s[k].set_visible(False)

      #3) Make room for everything:
      if self._title is not None:
         title_height = self.get_bbox(self._title).height
      else:
         title_height = 0
      if self._xlabel is not None:
         xlabel_height = self.get_bbox(self._xlabel).height
      else:
         xlabel_height = 0
      if self._ylabel is not None:
         ylabel_width = self.get_bbox(self._ylabel).width
      else:
         ylabel_width = 0

      axis_label_height = self.get_xlabels_bbox().height
      axis_label_width = self.get_ylabels_bbox().width

      top_height = title_height*1.0
      bottom_height = (axis_label_height + xlabel_height)*1.1
      left_width = (axis_label_width + ylabel_width)*1.0

      self.set_top_padding(top_height)
      self.set_bottom_padding(bottom_height)
      self.set_left_padding(left_width)

      plt.draw()

   def close(self):
      plt.close(self.fig)


