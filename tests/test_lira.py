import pytest
import matplotlib
import numpy as np
matplotlib.use("Agg")
try:
   import astropy
   haveastropy = pytest.mark.skipif(False, reason="astropy is installed")
except:
   haveastropy = pytest.mark.skipif(True, reason="astropy not installed")

@pytest.fixture
def snobj():
   import snpy
   s = snpy.get_sn('SN2006ax.txt')
   s.choose_model('max_model', stype='st')
   s.fit(['B','V'])
   return s

def test_lira_noparams(snobj):
   EBV,eEBV,slope,eslope = snobj.lira('B','V')
   assert (round(EBV,3) == round(0.0471892119,3) and
           round(eEBV,3) == round(0.004766264,3))

def test_lira_interp_model(snobj, interpolate=True):
   EBV,eEBV,slope,eslope = snobj.lira('B','V')
   assert (round(EBV,3) == round(0.047188195845567507,3) and
           round(eEBV,3) == round(0.00476626486707499,3))

