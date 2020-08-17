import pytest
from snpy import get_sn
from snpy.utils import fit1dcurve

allresults = {'chebyshev': 828.730,
           'hermite': 828.730,
           'hermiteE': 828.730,
           'hyperspline': 827.426,
           'laguerre': 828.730,
           'polynomial': 828.730,
           'spline': 827.290,
           'spline2': 827.426}

if fit1dcurve.gp == 'pymc':
   allresults['gp'] = 827.172
elif fit1dcurve.gp == 'sklearn':
   allresults['gp'] = 827.202

@pytest.fixture
def snobj():
   import snpy
   snobj = snpy.get_sn('SN2006ax.txt')
   return snobj

funcs = fit1dcurve.functions
results = {}
for func in funcs:
   if func in allresults:
      results[func] = allresults[func]

@pytest.mark.parametrize("func,Tmax", list(results.items()))
def test_interp(snobj, func, Tmax):
   snobj.B.template(method=func)
   assert round(Tmax,3) == round(snobj.B.Tmax,3)
