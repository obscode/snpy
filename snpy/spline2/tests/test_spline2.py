from __future__ import print_function
import pytest
import snpy.spline2 as spline2
import sys
from numpy import array,arange,sin,pi,absolute,all

#data = columns('example1')
#p = Plot()
#p.point(data[0],data[1])

@pytest.fixture
def spdata():
   x = arange(0,2*pi,0.1)
   y = sin(x)
   dy = y + 0.01
   tck = spline2.spline2(x, y, w=1.0/dy)
   return (x,y,tck)


def test_spline2(spdata):
   x,y,tck = spdata
   yeval = spline2.evalsp(x, tck)
   assert all(absolute(yeval-y) < 0.01)

def test_extrema(spdata):
   xe,ye,si = spline2.eval_extrema(tck)
   assert len(xe) == 2

def test_root(spdata):
   roots = spline2.eval_x(0, tck)
   print(roots)
   assert len(xe) == 3

