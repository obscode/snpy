#!/usr/bin/env python
from snpy import dm15temp
import os,sys
from numpy import *
#from scipy.interpolate import bisplrep
#from scipy.interpolate import bisplev
#from pygplot import *
import pickle
import scipy.interpolate as inter

base = os.path.dirname(dm15temp.__file__)
tck_file = os.path.join(base, 'tck.pickle')

temp = dm15temp.template()

dm15_low = 0.83
dm15_high = 1.93
t_low = -10.0
t_high = 70

bands = ['B','V','R','I']
dm15s = arange(31)/30.0*(dm15_high-dm15_low) + dm15_low
ts = arange(81)/80.0*(t_high-t_low) + t_low
#mp = Panel(3,3, device='tcl.ps/VCPS')

dm15_mat = []
t_mat = []
z_mat = {}
ez_mat = {}

for band in bands:
   z_mat[band] = []
   ez_mat[band] = []
   #p[band] = Plot(font=2, title=band, xlabel='Epoch', ylabel='m', flipyaxis=1)
   #mp.add(p[band])

for i in range(len(dm15s)):
   sys.stdout.write('%.2f ' % dm15s[i])
   sys.stdout.flush()
   temp.mktemplate(dm15s[i])
   dm15_mat.append(ts*0 + dm15s[i])
   t_mat.append(ts)
   for band in bands:
      z,ez,mask = temp.eval(band, ts)
      z_mat[band].append(z)
      ez_mat[band].append(ez)

x = ravel(array(t_mat))
y = ravel(array(dm15_mat))
if os.path.isfile(tck_file):
   f = open(tck_file)
   tck = pickle.load(f)
   f.close()
else:
   tck = {}
smooths = {}
for band in bands:  smooths[band] = 10
for band in bands:
   z = ravel(array(z_mat[band]))
   ez = ravel(array(ez_mat[band]))
   print "fitting band ",band
   print inter.__file__
   try:
      print x.shape
      print y.shape
      print z.shape
      #tck[band] = inter.bisplrep(x, y, z, w=1.0/ez, kx=1, ky=1, s=0.0*len(x))
      tck[band] = inter.SmoothBivariateSpline(x, y, z, w=1.0/ez, s=smooths[band]*len(x))
   except:
      print x,y,z,ez,smooths[band]*len(x)
      print type(x), x.shape
      print type(y), y.shape
      print type(z), z.shape
      print type(ez), ez.shape
      sys.exit(1)
   tck["e_"+band] = inter.SmoothBivariateSpline(x, y, ez, kx=1, ky=1, s=0.0*len(x))
f = open(tck_file, 'w')
pickle.dump(tck, f)
f.close()
