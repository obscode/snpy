import pytest
from snpy import get_sn
from snpy.utils import fit1dcurve
from snpy.version import __version__

@pytest.fixture
def snobj():
   import snpy
   snobj = snpy.get_sn('SN2006ax.txt')
   for f in snobj.data:
      snobj.data[f].template(method='spline2')
   return snobj

# Test computing of bolometric light-curve

def test_bolo_direct(snobj):
   x,y,f,lims = snobj.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='direct', interpolate='spline', SED='H3', verbose=True,
      EBVhost=0, Rv=3.1)
   checksum = sum(y)*1e-44
   assert round(checksum, 2) == round(2.24,2)

def test_bolo_direct_model(snobj):
   
   snobj.choose_model('color_model', stype='st')
   snobj.fit(dokcorr=False, Rv=3.1)
   x,y,f,lims = snobj.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='direct', interpolate='model', SED='H3')
   checksum = sum(y)*1.e-44
   assert round(checksum,2) == round(2.56, 2)

def test_bolo_SED(snobj):
   from snpy import fset
   l0 = fset['u'].wave.min()
   l1 = fset['H'].wave.max()
   x,y,f,lims = snobj.bolometric(['u','B','V','r','i','Y','J','H'], 
      method='SED', interpolate='spline', refband='B', lam1=l0,
      lam2=l1, SED='H3', EBVhost=0, Rv=3.1)
   checksum = sum(y)*1.e-44
   assert round(checksum,2) == round(2.71, 2)

