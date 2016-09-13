#!/usr/bin/env python

from Numeric import *
from pygplot import *
import spline2
import sys

data = columns('example1')
p = Plot()
p.point(data[0],data[1])

tck = spline2.spline2(data[0],data[1], w=1.0/data[2], verbose=1)
xs = arange(2001)/2000.0*(max(data[0]) - min(data[0])) + min(data[0])
ys = spline2.evalsp(xs, tck)
p.line(xs,ys, color='red', label='spline')
p.plot()

# Find the extrema:
xe,ye,si = spline2.eval_extrema(tck)
dy = (p.ymax - p.ymin)/50.0
gids = less(si, 0)
p.point(compress(gids, xe), compress(gids, ye)+dy, color='green', size=2,
      symbol=854, label='maxima')
gids = greater(si, 0)
p.point(compress(gids, xe), compress(gids, ye)-dy, color='green', size=2,
      symbol=852, label='minima')

# Find the roots:
value = 40000.0
roots = spline2.eval_x(value, tck)
p.line([200,2000],[value,value], color='yellow')
p.point(roots, roots*0+value, color='yellow', symbol=4, size=2, label='roots')
      
p.legend()
p.replot()
p.close()

data2 = columns('example2')
tck = spline2.spline2(data[0], data[1], w=1.0/data[2], verbose=1)
