'''A module of convenience functions to deal with SALT2 <--> SNooPy.'''
import os
from numpy import array


snpy_to_salt = {
   # CSPI
   'B':('SWOPE2','B','VEGA-CSP'),
   'V':('SWOPE2','V','VEGA-CSP'),
   'V0':('SWOPE2','V0','VEGA-CSP'),
   'V1':('SWOPE2','V1','VEGA-CSP'),
   'u':('SWOPE2','u','BD17-CSP'),
   'g':('SWOPE2','g','BD17-CSP'),
   'r':('SWOPE2','r','BD17-CSP'),
   'i':('SWOPE2','i','BD17-CSP'),
   # CSPII
   'B2':('SWOPE3','B','VEGA-CSP'),
   'V2':('SWOPE3','V','VEGA-CSP'),
   'u2':('SWOPE3','u','BD17-CSP'),
   'g2':('SWOPE3','g','BD17-CSP'),
   'r2':('SWOPE3','r','BD17-CSP'),
   'i2':('SWOPE3','i','BD17-CSP'),
   #Swift UVOT
   'U_UVOT':('UVOT','U','VEGA-CSP'),
   'B_UVOT':('UVOT','B','VEGA-CSP'),
   'V_UVOT':('UVOT','V','VEGA-CSP'),
   # Standard
   'Us':('STANDARD','U','VEGA2'),
   'Bs':('STANDARD','B','VEGA2'),
   'Vs':('STANDARD','V','VEGA2'),
   'Rs':('STANDARD','R','VEGA2'),
   'Is':('STANDARD','I','VEGA2'),
   # KAIT is also Standard
   'Ukait':('STANDARD','U','VEGA2'),
   'Bkait':('STANDARD','B','VEGA2'),
   'Vkait':('STANDARD','V','VEGA2'),
   'Rkait':('STANDARD','R','VEGA2'),
   'Ikait':('STANDARD','I','VEGA2'),
   # SDSS
   'u_s':('SDSS','u','AB_B12'),
   'g_s':('SDSS','g','AB_B12'),
   'r_s':('SDSS','r','AB_B12'),
   'i_s':('SDSS','i','AB_B12'),
   'z_s':('SDSS','z','AB_B12'),
   # KEPLERCAM
   'Bk':('KEPLERCAM','B','VEGA2'),
   'Vk':('KEPLERCAM','V','VEGA2'),
   'uk':('KEPLERCAM','u','VEGA2'),
   'rk':('KEPLERCAM','r','VEGA2'),
   'ik':('KEPLERCAM','i','VEGA2'),
   'Bk1':('KEPLERCAM','B','VEGA2'),
   'Vk1':('KEPLERCAM','V','VEGA2'),
   'uk1':('KEPLERCAM','u','VEGA2'),
   'rk1':('KEPLERCAM','r','VEGA2'),
   'ik1':('KEPLERCAM','i','VEGA2'),
   'Bk2':('KEPLERCAM','B','VEGA2'),
   'Vk2':('KEPLERCAM','V','VEGA2'),
   'uk2':('KEPLERCAM','u','VEGA2'),
   'rk2':('KEPLERCAM','r','VEGA2'),
   'ik2':('KEPLERCAM','i','VEGA2'),
   # 4-Shooter
   'U4sh':('4SHOOTER2','U','VEGA2'),
   'B4sh':('4SHOOTER2','B','VEGA2'),
   'V4sh':('4SHOOTER2','V','VEGA2'),
   'R4sh':('4SHOOTER2','R','VEGA2'),
   'I4sh':('4SHOOTER2','I','VEGA2'),
   #USNO
   'u_40':('USNO','u','VEGA2'),
   'g_40':('USNO','g','VEGA2'),
   'r_40':('USNO','r','VEGA2'),
   'i_40':('USNO','i','VEGA2'),
   'z_40':('USNO','z','VEGA2'),
   # PANSTARRS
   'ps1_u':('PS','u','PS-AB'),
   'ps1_g':('PS','g','PS-AB'),
   'ps1_r':('PS','r','PS-AB'),
   'ps1_i':('PS','i','PS-AB'),
   'ps1_z':('PS','z','PS-AB'),
   # ANDICAM
   'BANDI':('ANDICAM','B','VEGA-CSP'),
   'VANDI':('ANDICAM','V','VEGA-CSP'),
   'RANDI':('ANDICAM','R','VEGA-CSP'),
   'IANDI':('ANDICAM','I','VEGA-CSP'),
}   
   
# Stock SALT2 config
snpy_to_salt0 = snpy_to_salt.copy()
for filt in ['u','g','r','i','B']:
   snpy_to_salt0[filt] = ('SWOPE2',filt,'VEGA2')
for filt in ['u2','g2','r2','i2','B2','V2']:
   snpy_to_salt0[filt] = ('SWOPE2',filt[0],'VEGA2')
snpy_to_salt0['V0'] = ('SWOPE2','V', 'VEGA2')
snpy_to_salt0['V'] = ('SWOPE2','V2', 'VEGA2')
snpy_to_salt0['V1'] = ('SWOPE2','V1', 'VEGA2')
snpy_to_salt0['V2'] = ('SWOPE2','V', 'VEGA2')

salt_to_snpy = {}
for key in snpy_to_salt:
   inst,filt,syst = snpy_to_salt[key]
   newkey = "{}::{} {}".format(inst,filt,syst)
   salt_to_snpy[newkey] = key

salt_to_snpy0 = {}
for key in snpy_to_salt0:
   inst,filt,syst = snpy_to_salt0[key]
   newkey = "{}::{} {}".format(inst,filt,syst)
   salt_to_snpy0[newkey] = key

def parse_results(infile):
   '''Take the file output from SALT2 and parse the results into python
   data structures.
   
   Args:
      infile (str):  The file to parse
      
   Returns:
      (results,CovMat,stats):
      
      results:  a dictionary of 2-tuples, keyed by parameter and equal to 
                the best-fit value and error
      CovMat:   a dictionary keyed by [par1][par2], giving the covariance
                matrix element between par1 and par2.
      stats:    dictionary of SALT2 statistics'''

   if not os.path.isfile(infile):
      raise IOError('File not found: {}'.format(infile))

   with open(infile, 'r') as fin: lines = fin.readlines()
   results = {'Salt2Model':{},
              'UBVR':{}}
   CovMat = {'Salt2Model':{},
             'UBVR':{}}
   stats = {}

   current_model = None
   current_params = []
   for line in lines:
      if line.find('BEGIN_OF_FITPARAMS') == 0:
         current_model = line.strip().split()[1]
         current_params = []
         continue
      if line.find('END_OF_FITPARAMS') == 0:
         current_model = None
         continue
      if current_model is not None:
         fs = line.strip().split()
         if fs[0].find('Cov') == 0:
            # Covariance matrix result
            for par1 in current_params:
               for par2 in current_params:
                  if par1 in fs[0] and par2 in fs[0]:
                     C = float(fs[1])
                     if par1 in CovMat[current_model]:
                        CovMat[current_model][par1][par2] = C
                     else:
                        CovMat[current_model][par1] = {par2:C}
                     if par2 in CovMat[current_model]:
                        CovMat[current_model][par2][par1] = C
                     else:
                        CovMat[current_model][par2] = {par1:C}
         else:
            # parameter result
            current_params.append(fs[0])
            if fs[-1] == 'F':
               results[current_model][fs[0]] = (float(fs[1]),-1)
            else:
               results[current_model][fs[0]] = (float(fs[1]),float(fs[2]))
      else:
         if line[0] == '@':
            fs = line.strip().split()
            stats[fs[0][1:]] = float(fs[1])
   return results,CovMat,stats

# Convert from salt to snpy parameter names
saltpar_to_snpy = {'X0':'X0',
                   'X1':'X1',
                   'Color':'Color',
                   'RestFrameMag_0_U':'Umax',
                   'RestFrameMag_0_B':'Bmax',
                   'RestFrameMag_0_V':'Vmax',
                   'RestFrameMag_0_R':'Rmax',
                   'DayMax':'Tmax'
                   }

snpypar_to_salt = {}
for key in saltpar_to_snpy:
   snpypar_to_salt[saltpar_to_snpy[key]] = key

def parse_lc(infile, stock=True, trans=None):
   if not os.path.isfile(infile):
      raise IOError('File not found: {}'.format(infile))

   with open(infile, 'r') as fin: lines = fin.readlines()
   lines = [line for line in lines if line[0] != '@' and line[0] != '#']
   lcs = {}
   if trans is None:
      if stock:
         trans = salt_to_snpy0
      else:
         trans = salt_to_snpy
   for line in lines:
      MJD,flux,eflux,zp,filt,syst = line.strip().split()
      filt = trans[filt+" "+syst]
      if filt not in lcs:
         lcs[filt] = [[float(MJD), float(flux),float(eflux), float(zp)]]
      else:
         lcs[filt].append([float(MJD), float(flux),float(eflux), float(zp)])
   for filt in lcs:
      lcs[filt] = array(lcs[filt]).T

   return lcs
      
