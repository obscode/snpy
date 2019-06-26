import matplotlib
matplotlib.use("Agg")
from pymc import *
from pymc.gp import *
from pymc.gp.cov_funs import matern
from mydistance import aniso_euclidean_s
from getSNdata2 import getdata
from pylab import *
#from FITS import qdump
from astropy.io import fits
from pylab import *
import sys
import pickle
from numpy import *

matern.add_distance_metric('aniso_euclidean_s', 'mydistance', with_x=False)

class SNmodel(MCMC):

   def __init__(self, filter, scalelims=(1,10), use_B_tmax=True, 
         min_t_scale=1., sbetalims=(0.5,10,0.6), diff_degree=None, 
         data_location='.', **args):

      x,y,v,ev,mesh,xoff,xscale,yoff,yscale,foff,names = \
            getdata(filter, use_B_tmax=use_B_tmax, data_location=data_location)
      vv = power(ev,2)
      self.filter = filter
      # ============
      # = The mean =
      # ============
      
      # Vague prior for the overall mean
      mn = Uninformative('mn', value=0)
      
      #def constant(x, val):
      #    return zeros(x.shape[:-1],dtype=float) + val

      def constant(x, val):
         return zeros(x.shape[:-1], dtype=float) + val
      
      @deterministic
      def M(mn=mn):
         return Mean(constant, val=mn)
      
      # ==================
      # = The covariance =
      # ==================
      
      # Informative prior for the degree of differentiability
      if diff_degree is None:
         diff_degree = Uniform('diff_degree',2.,10.,value=3)
      
      # Vague priors for the amplitude and scale
      amp = Exponential('amp',7.e-5, value=np.std(v))
      #scale = Exponential('scale',4e-3, value=x.max()/2)
      #scale = Uniform('scale', 10.0/xscale, 5., value=x.max()/2)
      salpha = 1
      if type(sbetalims) is type ((1,)):
         sbeta = Uniform('sbeta', sbetalims[0], sbetalims[1], 
               value=sbetalims[2])
      else:
         sbeta = sbetalims
      print salpha, sbeta, x.max()/2
      scale = InverseGamma('scale', salpha, sbeta, value=x.max()/2)
      if type(scalelims) is type ((1,)):
         scalerat = Uniform('scalerat', scalelims[0], scalelims[1], 
               value=random.uniform(scalelims[0], scalelims[1], size=1))
      else:
         scalerat = array([scalelims])
      @deterministic
      def C(amp=amp, scale=scale, diff_degree=diff_degree, sr=scalerat):
          return FullRankCovariance(matern.aniso_euclidean_s, 
                diff_degree=diff_degree, amp=amp, scale=scale,
                sr=sr)
      
      # ===================
      # = The GP submodel =
      # ===================
      lc_v = GPSubmodel('lc_v', M, C, mesh)
      
      # ============
      # = The data =
      # ============
      
      valpha = 1
      vbeta = Uniform('vbeta', 1e-6, 1, value=0.0005)
      # Vague prior for the observation variance
      var = Exponential('var',5e-9, value=np.std(v)**2)

      @deterministic
      def tau(vv=vv, var=var):
         return power(vv + var,-1)
      
      # This is the observation of the elevation variable 'v', plus 
      #   normally-distributed error.
      d = Normal('d',lc_v.f_eval,tau,value=v,observed=True)
      MCMC.__init__(self, locals(), **args)
      self.use_step_method(GPEvaluationGibbs, lc_v, var, d)

   def dump_traces(self, outfile):
      f = open(outfile, 'w')
      d = {}
      d['amp'] = self.amp.trace()
      try:
         d['sbeta'] = self.sbeta.trace()
      except:
         pass
      d['scale'] = self.scale.trace()
      try:
         d['scalerat'] = self.scalerat.trace()
      except:
         pass
      try:
         d['diff_degree'] = self.diff_degree.trace()
      except:
         pass
      d['var'] = self.var.trace()
      pickle.dump(d, f)
      f.close()

   def make_mean_surf(self,base='./',m=101):
      n = len(self.trace('var')[:])
      x = self.x
      y = self.y
      xplot = linspace(-0.1,1.1, m)
      yplot = linspace(-0.1,1.1, m)
      dplot = dstack(meshgrid(xplot, yplot))

      Msurf = zeros(dplot.shape[:2])
      E2surf = zeros(dplot.shape[:2])

      for i in xrange(n):
         self.remember(0,i)
         Msurf_i,Vsurf_i = point_eval(self.lc_v.M_obs.value,
               self.lc_v.C_obs.value, dplot)
         Msurf += Msurf_i/n
         E2surf += (Vsurf_i + Msurf_i**2)/n
      Vsurf = E2surf - Msurf**2
      SDsurf = sqrt(Vsurf)
      Msurf += self.foff

      x = self.x*self.xscale + self.xoff
      y = self.y*self.yscale + self.yoff
      x0 = xplot[0]*self.xscale + self.xoff
      x1 = xplot[-1]*self.xscale + self.xoff
      y0 = yplot[0]*self.yscale + self.yoff
      y1 = yplot[-1]*self.yscale + self.yoff

      # Plot mean and standard deviation surfaces
      close('all')
      imshow(Msurf, extent=[x0,x1,y0,y1],
            interpolation='nearest', aspect='auto', origin='lower')
      plot(x,y,'r.',markersize=4)
      axis([x0,x1,y0,y1])
      title('Posterior predictive mean surface')
      xlabel('phase (days)')
      ylabel('stretch')
      colorbar()
      savefig('%slc2_surf_mean_%s.pdf' % (base,self.filter))
      
      figure()
      imshow(SDsurf, extent=[x0,x1,y0,y1],
            interpolation='nearest', aspect='auto', origin='lower')
      plot(x,y,'r.',markersize=4)
      axis([x0,x1,y0,y1])
      title('Posterior predictive standard deviation surface')
      xlabel('phase (days)')
      ylabel('stretch')
      colorbar()
      savefig('%slc2_surf_var_%s.pdf' % (base,self.filter))
      
      N = Msurf.shape[0]
      hdu = fits.PrimaryHDU(Msurf)
      hdu.header['CRPIX1'] = 1
      hdu.header['CRPIX2'] = 1
      hdu.header['CRVAL1'] = x0
      hdu.header['CRVAL2'] = y0
      hdu.header['CDELT1'] = (x1-x0)/N
      hdu.header['CDELT2'] = (y1-y0)/N
      fts = fits.HDUList([hdu])
      fts.writeto('%s%s_mean2.fits' % (base,self.filter))

      hdu = fits.PrimaryHDU(SDsurf)
      hdu.header['CRPIX1'] = 1
      hdu.header['CRPIX2'] = 1
      hdu.header['CRVAL1'] = x0
      hdu.header['CRVAL2'] = y0
      hdu.header['CDELT1'] = (x1-x0)/N
      hdu.header['CDELT2'] = (y1-y0)/N
      fts = fits.HDUList([hdu])
      fts.writeto('%s%s_mean2.fits' % (base,self.filter))

      #extras = [['CRPIX1',1],['CRVAL1',x0],['CDELT1',(x1-x0)/N],
      #          ['CRPIX2',1],['CRVAL2',y0],['CDELT2',(y1-y0)/N]]
      
      #qdump('%s%s_mean2.fits' % (base,self.filter), Msurf, extras=extras)
      #qdump('%s%s_std2.fits' % (base,self.filter), SDsurf, extras=extras)
      
if __name__ == "__main__":
   #from basics import rdarg
   import argparse
   import string

   parser = argparse.ArgumentParser(description="Generate a LC templates")
   parser.add_argument('filter', help="filter to build a template")
   parser.add_argument('-N', help="Number of MCMC interations", 
         type=int, default=100000)
   parser.add_argument('-burn', help="Number of burn-in iterations",
         type=int, default=0)
   parser.add_argument('-thin', help="Thinning factor", type=int,
         default=1)
   parser.add_argument('-scalelim', help="Lower and upper limit of scale param",
         nargs=2, type=float, default=[1.0,10.0])
   parser.add_argument('-scale', help="Fix scale parameter to this value",
         nargs=2, type=float, default=None)
   parser.add_argument('-sbetalim', help="Lower and upper limit of sbeta param",
         nargs=2, type=float, default=[0.5,10.0])
   parser.add_argument('-sbeta', help="Starting value for sbeta param",
         type=float, default=0.6)
   parser.add_argument('-dd', help="Fix dd parameter to this value",
         type=float, default=None)
   parser.add_argument('-Bmax', help="Use time of B maximum as t=0",
         action='store_true')
   parser.add_argument('-min_t_scale', help="Minimum time scale parameter",
         type=float, default=5.0)
   parser.add_argument('-datadir', help="Location of SNooPy fit files",
         default='.')
   parser.add_argument('-out', help="Base name of output files", 
         default = './')
   #argv,scalelims = rdarg(argv, '-scalelim', None, "1:10")
   #argv,sbetalims = rdarg(argv, '-sbeta', None, "0.5:10:0.6")
   #argv,diff_degree = rdarg(argv, '-dd', float, None)
   #argv,use_B_tmax = rdarg(argv, '-Bmax', int, 0, single=1)
   #argv,min_t_scale = rdarg(argv, '-min_t_scale', float, 5.0)
   #argv,output = rdarg(argv, '-out', None, './')
   #filter = argv[0]
   args = parser.parse_args()

   if args.scale is not None:
      scalelims = args.scale
   else:
      scalelims = tuple(args.scalelim)
   print scalelims
   #if len(scalelims) == 2:
   #   scalelims = tuple(map(float,scalelims))
   #else:
   #   scalelims = float(scalelims[0])

   sbetalims = tuple(list(args.sbetalim) + [args.sbeta])

   sampler = SNmodel(args.filter, scalelims=scalelims, use_B_tmax=args.Bmax,
         min_t_scale=args.min_t_scale, sbetalims=sbetalims, 
         diff_degree=args.dd, data_location=args.datadir)
   sampler.sample(args.N,burn=args.burn,thin=args.thin)
   sampler.make_mean_surf(base=args.out)
   sampler.dump_traces('%straces_%s.pickle' % (args.out,sampler.filter
      ))
