import pytest
from snpy import get_sn
import numpy as num
import os

# Test if SALT2 is available
havesalt = pytest.mark.skipif('SALTPATH' not in os.environ,
      reason="SALTPATH environment variable not set (SALT2 not installed?)")

havemlcs = pytest.mark.skipif('MLCS2K2_BASEDIR' not in os.environ,
    reason="MLCS2K2_BASEDIR environment variable not set (MLCS not installed?")

@pytest.fixture
def snobj():
   import snpy
   snobj = snpy.get_sn('SN2006ax.txt')
   return snobj

# Test fitting all the models. First, here are the results we expect

results = {
      'EBV_model':
          {'DM': 34.159,
           'EBVhost': -0.003,
           'Tmax': 827.671,
           'dm15': 1.086},
      'EBV_model2':
           {'DM': 34.254,
            'EBVhost': 0.0091,
            'Tmax': 827.564,
            'st': 0.962},
      'max_model':
           {'Bmax': 15.001,
            'Hmax': 15.991,
            'Jmax': 15.642,
            'Tmax': 827.547,
            'Vmax': 15.068,
            'Ymax': 15.757,
            'gmax': 14.957,
            'imax': 15.695,
            'rmax': 15.123,
            'st': 0.969,
            'umax': 15.399},
      'color_model':
           {'Bmax': 14.834,
            'EBVhost': 0.040,
            'Rv': 3.1,
            'Tmax': 827.565,
            'st': 0.962},
}
SALT_results = {
            'X0':0.018,
            'X1':0.121,
            'Color':-0.105,
            'Tmax':827.981,
            'Bmax':15.011,
}

MLCS_results = {
      'del':-0.107,
      'av0':0.099,
      'DM':34.548,
      'Tmax':827.31,
      'Vmax':15.084
}

@pytest.mark.parametrize("model,result", list(results.items()))
def test_model(snobj, model, result):
   if model == 'EBV_model':  
      stype='dm15'
   else:
      stype = 'st'
     
   snobj.replot=False
   snobj.choose_model(model, stype=stype)
   if model == 'color_model':
      snobj.fit(dokcorr=False, Rv=3.1)
   else:
      snobj.fit(dokcorr=False)

   res = [round(snobj.model.parameters[key],3) == round(result[key],3) \
         for key in result]
   assert num.all(res)

@havesalt
def test_SALT(snobj):
   snobj.replot=False
   snobj.choose_model('SALT_model', workdir='temp')
   snobj.fit()

   res = [round(snobj.model.parameters[key],3) == round(SALT_results[key],3) \
         for key in SALT_results]
   assert num.all(res)

@havemlcs
def test_MLCS(snobj):
   snobj.replot=False
   snobj.choose_model('MLCS_model')
   snobj.fit()

   res = [round(snobj.model.parameters[key],3) == round(MLCS_results[key],3) \
         for key in MLCS_results]
   assert num.all(res)
