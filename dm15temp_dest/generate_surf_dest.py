import dm15temp_dest
import os
from Numeric import *
from scipy.interpolate import bisplev
from scipy.interpolate import bisplrep
from pygplot import *
import pickle

temp = dm15temp_dest.template()
dm15temp_dest.use_fits=1

dm15_low = 0.6
dm15_high = 2.0
t_low = -10.0
t_high = 70

bands = ['d%d'%i for i in range(15)]
bands = bands + ['dw%d'%i for i in range(10)]
dm15s = arange(21)/20.0*(dm15_high-dm15_low) + dm15_low
ts = arange(81)/80.0*(t_high-t_low) + t_low
mp = Panel(2,3, device='tcl.ps/VCPS')

dm15_mat = []
t_mat = []
z_mat = {}
ez_mat = {}
p = {}

for band in bands:
   z_mat[band] = []
   ez_mat[band] = []
   p[band] = Plot(font=2, title=band, xlabel='Epoch', ylabel='m', flipyaxis=1)
   mp.add(p[band])

for i in range(len(dm15s)):
   sys.stdout.write('%.2f ' % dm15s[i])
   sys.stdout.flush()
   temp.mktemplate(dm15s[i])
   dm15_mat.append(ts*0 + dm15s[i])
   t_mat.append(ts)
   for band in bands:
      z,ez,mask = temp.eval(band, ts, mag=0)
      z_mat[band].append(z)
      ez_mat[band].append(ez)
      p[band].line(ts, z+5*dm15s[i])

x = ravel(array(t_mat))
y = ravel(array(dm15_mat))
if os.path.isfile('tck_destiny.pickle'):
   f = open('tck_destiny.pickle')
   tck = pickle.load(f)
   f.close()
else:
   tck = {}
for band in bands:
   z = ravel(array(z_mat[band]))
   ez = sqrt(ravel(array(ez_mat[band]))**2 + 0.05**2)
   print "fitting band ",band
   tck[band] = bisplrep(x, y, z, w=1.0/ez, s=0.5*len(x))
   tck["e_"+band] = bisplrep(x, y, ez, kx=1, ky=1, s=0.0*len(x))
f = open('tck_destiny.pickle', 'w')
pickle.dump(tck, f)
f.close()

ts = arange(321)/320.0*(t_high-t_low) + t_low
for i in range(len(dm15s)):
   for band in bands:
      z = bisplev(ts, dm15s[i], tck[band])
      p[band].line(ts, z[:,0]+5*dm15s[i], color='red')

mp.plot()
mp.close()
