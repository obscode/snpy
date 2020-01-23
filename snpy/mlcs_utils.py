'''A module of convenience functions to deal with MLCS2k2<-->SNooPy'''
import os
from numpy import array
from glob import glob


snpy_to_mlcs = {}
for f in ['B','V','u','g','r','i','Y','J','H']:
   snpy_to_mlcs[f] = f+"_CSP"
   snpy_to_mlcs[f+'2'] = f+"_CSP2"
snpy_to_mlcs['V0'] = 'V0_CSP'
snpy_to_mlcs['V1'] = 'V1_CSP'
snpy_to_mlcs.update({
   'U_UVOT':'uu_Swift',
   'B_UVOT':'bb_Swift',
   'V_UVOT':'vv_Swift'})
for f in ['U','B','V','R','I']:
   snpy_to_mlcs[f+"s"] = f+"_Astier"
   snpy_to_mlcs[f+"kait"] = f+"_Astier"

for f in ['u','g','r','i','z']:
   snpy_to_mlcs[f+"_s"] = f+"_SDSS"
   snpy_to_mlcs[f+"_40"] = f+"_AB"

# Stock MLCSk2 config
snpy_to_mlcs0 = snpy_to_mlcs.copy()
for filt in ['u2','g2','r2','i2','B2','V2']:
   snpy_to_mlcs0[filt] = filt+"_CSP"

mlcs_to_snpy = {}
for key in snpy_to_mlcs:
   filt = snpy_to_mlcs[key]
   mlcs_to_snpy[filt] = key

mlcs_to_snpy0 = {}
for key in snpy_to_mlcs0:
   filt = snpy_to_mlcs0[key]
   mlcs_to_snpy[filt] = key

def parse_results(loc):
   '''Take the file output from MLCS and parse the results into python
   data structures.
   
   Args:
      loc (str):  Location (folder) for the output files
      
   Returns:
      (results,CovMat,stats):
      
      results:  dictionary of  dictionaries of 2-tuples, keyed by interation
                and parameter name equal to the best-fit value and error.
      stats:  dictionary of dictionaries of statistics, keyed by iteration
              and statistic'''


   if not os.path.isdir(loc):
      raise IOError('Folder not found: {}'.format(loc))
   fits = glob(os.path.join(loc, '*.iter*.fit.out'))
   if len(fits) == 0:
      # No fit results, probably MLCS died
      return None,None
   fits.sort()
   iters = [os.path.basename(fit).split('.')[1] for fit in fits]
   base = os.path.basename(fits[0]).split('.')[0]

   tags = ['mu0','av0','del','t0','Rv','m0v','p0']
   results = {}
   stats = {}
   for it in iters:
      fit = os.path.join(loc,"{}.{}.fit.out".format(base,it))
      with open(fit, 'r') as fin: lines = fin.readlines()
      lines = [line.strip().split() for line in lines if line[0] != ';']
      results[it] = {}
      for i,tag in enumerate(tags):
         results[it][tag] = list(map(float,lines[-1][2*i+1:2*i+3]))

   with open(os.path.join(loc, base+".mcs.log"), 'r') as fin:
      lines = fin.readlines()
   lines = [line.strip().split() for line in lines if line[0] != ';']
   for line in lines:
      it,rchisq,ndof,chisq0,chisq,prisq,_,_,_,_,_,_,_,typ = line
      it = "iter{}".format(it)
      stats[it] = dict(rchisq=float(rchisq), ndof=int(ndof), 
            chisq=float(chisq),type=typ)
   return results,stats

# Convert from mlcs to snpy parameter names
mlcspar_to_snpy = {'del':'del',
                   'av0':'av0',
                   'mu0':'DM',
                   't0':'Tmax',
                   'Rv':'Rv',
                   'm0v':'Vmax',
                   'p0':'p0',
                   }

snpypar_to_mlcs = {}
for key in mlcspar_to_snpy:
   snpypar_to_mlcs[mlcspar_to_snpy[key]] = key

def parse_lc(loc, stock=True, trans=None):
   if not os.path.isdir(loc):
      raise IOError('Folder not found: {}'.format(loc))

   if trans is None:
      if stock:
         trans = mlcs_to_snpy0
      else:
         trans = mlcs_to_snpy

   lcs = {}
   files = glob(os.path.join(loc, '*.obsframe.dat'))
   for fi in files:
      with open(fi) as fin: lines = fin.readlines()
      it = os.path.basename(fi).split('.')[1]
      lines = [line.strip().split() for line in lines if line[0] != '#']
      lcs[it] = {}
      for line in lines:
         filt,t,mag,emag,kcorr,ekcorr,galext,egalext = line
         filt = trans[filt]
         if filt not in lcs[it]:
            lcs[it][filt] = [[float(t),float(mag),float(emag),float(kcorr),
                          float(ekcorr),float(galext),float(egalext)]]
         else:
            lcs[it][filt].append([float(t),float(mag),float(emag),float(kcorr),
                          float(ekcorr),float(galext),float(egalext)])
      for filt in lcs[it]:
         lcs[it][filt] = array(lcs[it][filt]).T
   return lcs
