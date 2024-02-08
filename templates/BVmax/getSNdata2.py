from numpy import *
from snpy import *
import os,sys
from astropy.io import ascii

data = ascii.read('BV_max.dat')
names = list(data['col1'])
BVmax = array(data['col8'])


def getdata(filt, use_B_tmax=True, data_location='.'):
   flist = [suff+'_GP.snpy' for suff in names]

   ss = []
   for f in flist:
      fil = os.path.join(data_location,f)
      if not os.path.isfile(fil):
         print "%s not found" % (fil)
         continue
      ss.append(get_sn(fil))
   
   sts = []
   ts = []
   fluxes = []
   efluxes = []
   SNnames = []
   
   for i in range(len(ss)):
      s = ss[i]
      if filt not in s.data: 
         continue
      if s.B.Tmax is None or s.B.Tmax < 10:
         continue
      try:
         Tmax,Mmax,eMmax,rb = s.get_max(filt)
      except:
         continue
      
      this_flux = power(10, -0.4*(s.data[filt].mag - Mmax))
      this_eflux = s.data[filt].e_mag*this_flux/1.087
      this_flux = concatenate([this_flux, [1.0]])
      this_eflux = concatenate([this_eflux, [this_eflux.min()/10]])
      if use_B_tmax:
         this_t = (s.data[filt].MJD - s.B.Tmax)/(1 + s.z)
      else:
         this_t = (s.data[filt].MJD - Tmax)/(1 + s.z)
      this_t = concatenate([this_t, [0.0]])
      sids = argsort(this_t)
      this_flux = this_flux[sids]
      this_eflux = this_eflux[sids]
      this_t = this_t[sids]
      # get rid of repeated points:
      dts = this_t[1:] - this_t[:-1]
      bids = concatenate([less(absolute(dts), 2e-1),[False]])
      print bids.shape
      print this_t.shape
      while any(bids):
         for id in nonzero(bids)[0]:
            print "combining point %d and %d for %s" % (id+1,id,s.name)
            this_flux[id+1] = (this_flux[id+1]+this_flux[id])/2
            this_eflux[id+1] = (this_eflux[id+1]+this_eflux[id])/2
         this_flux = this_flux[~bids]
         this_eflux = this_eflux[~bids]
         this_t = this_t[~bids]
         dts = this_t[1:] - this_t[:-1]
         bids = concatenate([less(absolute(dts), 2e-1),[False]])
      ts = concatenate([ts,this_t])
      fluxes = concatenate([fluxes,this_flux])
      efluxes = concatenate([efluxes,this_eflux])
      sts = concatenate([sts, this_t*0+BVmax[i]/30.])
      SNnames = SNnames + [s.name]*this_t.shape[0]
   
   sts = array(sts)
   ts = array(ts)
   fluxes = array(fluxes)
   efluxes = array(efluxes)
   
   gids = greater_equal(ts, -10)*less_equal(ts,70)
   ts = ts[gids]
   sts = sts[gids]
   fluxes = fluxes[gids]
   efluxes = efluxes[gids]
   SNnames = [SNnames[i] for i in range(len(gids)) if gids[i]]
   
   # re-scale to 0,1 in both dimentions
   xoff = ts.min()
   xscale = (ts.max() - ts.min())
   x = (ts - xoff)/xscale
   
   yoff = sts.min()
   yscale = (sts.max() - sts.min())
   y = (sts - yoff)/yscale
    
   mesh = vstack((x,y)).T
   mean_flux = mean(fluxes)
   v = fluxes-mean_flux
   return x,y,v,efluxes,mesh,xoff,xscale,yoff,yscale,mean_flux,SNnames
