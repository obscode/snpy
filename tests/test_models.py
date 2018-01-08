import pytest
from snpy import get_sn
import numpy as num

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
            'st': 0.962}}

@pytest.mark.parametrize("model,result", results.items())
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
   assert num.alltrue(res)
