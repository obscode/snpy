'''This script reads in a trace file and then re-does the mean surface
and variances, but this time bootstrapping the Sne that are loaded.'''
import matplotlib
matplotlib.use("Agg")
from pymc import *
from pymc.gp import *
from pymc.gp.cov_funs import matern
from getSNdata2 import getdata
from pylab import *
import pymc.database
from FITS import qdump
from pylab import *
import sys
import pickle
from numpy import *

f = open(sys.argv[1])
trs = pickle.load(f)
filter = sys.argv[2]
base = 'boot_'
m = 101

x,y,v,ev,mesh,xoff,xscale,yoff,yscale,foff = getdata(filter, 
      use_B_tmax=False)
vv = power(ev,2)

# ids into the data that show where SNe switch
ids = concatenate([nonzero(diff(y))[0],array([len(y)-1])])
NSN = ids.shape[0]
low = 0
sids = zeros(x.shape, dtype=int16)
for i in range(len(ids)):
   high = ids[i]
   sids[low:high+1] = i
   low = high+1

scales = trs['scale']
amps = trs['amp']
dds = trs['diff_degree']
if 'scalerat' in trs:
   scalerats = trs['scalerat']
   scalerat = median(scalerats)
else:
   scalerat = 5
vars = trs['var']

scale = median(scales)
amp = median(amps)
dd = median(dds)
var = median(vars)

mn = 0

def constant(x, val):
   return zeros(x.shape[:-1], dtype=float) + val

n = 500
xplot = linspace(-0.1,1.1, m)
yplot = linspace(-0.1,1.1, m)
dplot = dstack(meshgrid(xplot, yplot))

Msurf = zeros((n,m,m))

for i in range(n):

   print i
   M = Mean(constant, val=mn)
   C = Covariance(matern.aniso_euclidean_s, diff_degree=dd, amp=amp, 
         scale=scale, scalerat=scalerat)
   this_ids = array([], dtype=int64)
   for id in random.randint(NSN, size=(NSN,)):
      this_ids = concatenate([this_ids,nonzero(equal(id,sids))[0]])

   this_v = v[this_ids]
   this_var = power(ev[this_ids],2) + var
   this_mesh = take(mesh, this_ids, axis=0)

   observe(M, C, obs_mesh=this_mesh, obs_vals=this_v, obs_V=this_var)

   Msurf_i,Vsurf_i = point_eval(M, C, dplot)
   #maxs = Msurf_i.max(axis=1)
   #Msurf_i = Msurf_i/maxs[:,newaxis]

   Msurf[i,:,:] = Msurf_i


#Vsurf = E2surf - Msurf**2
med = median(Msurf, axis=0)
mad = 1.49*(median(absolute(Msurf - med[newaxis,:,:]),axis=0))

x = x*xscale + xoff
y = y*yscale + yoff
x0 = xplot[0]*xscale + xoff
x1 = xplot[-1]*xscale + xoff
y0 = yplot[0]*yscale + yoff
y1 = yplot[-1]*yscale + yoff

figure()
imshow(mad, extent=[x0,x1,y0,y1],
   interpolation='nearest', aspect='auto', origin='lower')
plot(x,y,'r.',markersize=4)
axis([x0,x1,y0,y1])
title('Posterior predictive standard deviation surface')
xlabel('phase (days)')
ylabel('stretch')
colorbar()
savefig('boot_var_%s.pdf' % (filter))
      
N = mad.shape[0]
extras = [['CRPIX1',1],['CRVAL1',x0],['CDELT1',(x1-x0)/N],
          ['CRPIX2',1],['CRVAL2',y0],['CDELT2',(y1-y0)/N]]
      
#qdump('%s%s_mean2.fits' % (base,filter), Msurf, extras=extras)
qdump('boot_%s_std2.fits' % (filter), mad, extras=extras)
      
