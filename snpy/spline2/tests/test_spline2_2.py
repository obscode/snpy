#!/usr/bin/env python

from numpy import *
from pygplot import *
import spline2
import sys
from numpy import random as RA
import scipy

xs = arange(0, 2*pi, 0.1)
ys = sin(xs)
ys = ys + RA.normal(0, 0.1, size=ys.shape[0])

mp = Panel(2,1)
p = Plot()
p.point(xs,ys, symbol=4)
p.error(xs,ys,dy1=0*xs+0.1)
mp.add(p)

p = Plot()
p.point(xs,ys, symbol=4)
p.error(xs,ys,dy1=0*xs+0.01)
mp.add(p)

tck = spline2.spline2(xs,ys, w=xs*0+1.0/0.1, verbose=1)
xxs = arange(2001)/2000.0*(max(xs) - min(xs)) + min(xs)
yys = spline2.evalsp(xxs, tck)
mp.plots[0].line(xxs,yys, color='red', label='spline2')

tck = scipy.interpolate.splrep(xs, ys, w=xs*0+1.0/0.1, task=0, s=len(xs))
yys = scipy.interpolate.splev(xxs, tck)
mp.plots[0].line(xxs,yys, color='blue', label='scipy')

tck = spline2.spline2(xs,ys, w=xs*0+1.0/0.01, verbose=1)
xxs = arange(2001)/2000.0*(max(xs) - min(xs)) + min(xs)
yys = spline2.evalsp(xxs, tck)
mp.plots[1].line(xxs,yys, color='red', label='spline2')

tck = scipy.interpolate.splrep(xs, ys, w=xs*0+1.0/0.01, task=0, s=len(xs))
yys = scipy.interpolate.splev(xxs, tck)
mp.plots[1].line(xxs,yys, color='blue', label='scipy')

mp.plots[0].legend()
mp.plots[1].legend()

mp.plots[0].line(xxs, sin(xxs), color='green', label='real')
mp.plots[1].line(xxs, sin(xxs), color='green', label='real')
mp.plot()
mp.close()
mp.device = 'comp.ps/CPS'
mp.plot()
mp.close()
