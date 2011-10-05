#!/usr/bin/env python
from snpy import CSPtemp
import os,sys
from numpy import *
from pygplot import *
import pickle
import scipy.interpolate as inter

base = os.path.dirname(CSPtemp.__file__)
tck_file = "bs_tck.pickle"

CSPtemp.use_gloes=1
temp = CSPtemp.template()

dm15_low = 0.6
dm15_high = 2.0
t_low = -10.0
t_high = 70

bands = ['u','B','g','V','r','i','Y','J','H']
dm15s = arange(31)/30.0*(dm15_high-dm15_low) + dm15_low
ts = arange(81)/80.0*(t_high-t_low) + t_low

data0 = {}
list0 = {}
SNe = {}
for band in bands:
   data0[band],list0[band] = CSPtemp.dm15tempc.get_data(band)
   list0[band] = array(list0[band])
   SNe[band] = {}.fromkeys(list0[band]).keys()
   SNe[band].sort()

tck = {}.fromkeys(bands)
for b in bands:  
   tck[b] = []
   tck['e_'+b] = []

for j in range(50):
   print "iteration ", j
   dm15_mat = []
   t_mat = []
   z_mat = {}
   ez_mat = {}
   for band in bands:
      z_mat[band] = []
      ez_mat[band] = []
   for i in range(len(dm15s)):
      sys.stdout.write('%.2f ' % dm15s[i])
      sys.stdout.flush()
      temp.mktemplate(dm15s[i])
      dm15_mat.append(ts*0 + dm15s[i])
      t_mat.append(ts)
      for band in bands:
         N = len(SNe[band])
         randints = random.randint(0, N, size=(N,))
         data1 = []
         sn_list = [SNe[band][id] for id in randints]
         masks = [list0[band] == sn for sn in sn_list]
         data1 = concatenate([compress(m, data0[band], axis=1) for m in masks],
               axis=1)
         list1 = concatenate([list0[band][m] for m in masks])
         list1 = list1.tolist()
         CSPtemp.dm15tempc.set_data(band, data1, list1)

         z,ez,mask = temp.eval(band, ts, mag=0)
         z_mat[band].append(z)
         ez_mat[band].append(ez)
   
   x = ravel(array(t_mat))
   y = ravel(array(dm15_mat))
   smooths = {}
   for band in bands:  
      smooths[band] = 0.1
      z = ravel(array(z_mat[band]))
      ez = ravel(array(ez_mat[band]))
      print "fitting band ",band
      print inter.__file__
      try:
         tck[band].append(inter.bisplrep(x, y, z, w=1.0/ez, s=smooths[band]*len(x)))
      except:
         continue
      tck["e_"+band].append(inter.bisplrep(x, y, ez, kx=1, ky=1, s=0.0*len(x)))
f = open(tck_file, 'w')
pickle.dump(tck, f)
f.close()
