import pytest
import os
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
   return snpy.get_sn('SN2006ax.txt')

def test_load_object(snobj):
   from snpy import sn
   assert isinstance(snobj,sn)

def test_object_equality(snobj):
   from snpy import get_sn
   sn2 = get_sn('SN2006ax.txt')
   assert sn2 == snobj

def test_sn_inequality(snobj):
   from snpy import get_sn
   sn2 = get_sn('SN2006ax.txt')
   sn2.B.MJD[0] += 0.1
   assert sn2 != snobj

def test_sn_zhelio(snobj):
   assert round(snobj.z,3) == round(0.016727,3)

def test_sn_zcmb(snobj):
   assert round(snobj.zcmb,7) == round(0.0179567,7)

@haveastropy
def test_sn_distmod(snobj):
   assert round(snobj.get_distmod(),3) == round(34.331, 3)

def test_sn_EBVgal(snobj):
   assert round(snobj.EBVgal,3) == round(0.041,3)

def test_get_mag_table(snobj):
   tab = snobj.get_mag_table()
   expected_keys = ['MJD']
   for f in snobj.data:
      expected_keys.append(f)
      expected_keys.append('e_'+f)
   assert set(tab.keys()) == set(expected_keys)

def test_get_mag_table_out(snobj):
   snobj.get_mag_table(outfile='table.dat')
   assert os.path.isfile('table.dat')

def test_make_template_spline(snobj):
   from snpy.utils.fit1dcurve import oneDcurve
   snobj.data['B'].template(method='spline')
   assert isinstance(snobj.B.interp, oneDcurve)

def test_make_template_spline(snobj):
   from snpy.utils.fit1dcurve import oneDcurve
   snobj.data['B'].template(method='spline2')
   assert isinstance(snobj.B.interp, oneDcurve)

def test_get_color(snobj):
   res = snobj.get_color('B','V', interp=False)
   assert len(res)==4

def test_dump_lc(snobj):
   import os
   snobj.dump_lc()
   assert os.path.isfile(snobj.name+"_lc_B_data.dat")
   os.system('rm SN2006ax*.dat')

def test_plot_sn(snobj):
   from snpy.myplotlib import PanelPlot
   ret = snobj.plot()
   assert isinstance(ret, PanelPlot)

def test_plot_filter(snobj):
   from snpy.myplotlib import PanelPlot
   ret = snobj.B.plot()
   assert isinstance(ret, PanelPlot)

def test_plot_filters(snobj):
   from matplotlib.figure import Figure
   ret = snobj.plot_filters()
   assert isinstance(ret, Figure)

def test_plot_color(snobj):
   from matplotlib.figure import Figure
   ret = snobj.plot_color('B','V')
   assert isinstance(ret, Figure)

def test_compute_w(snobj):
   w = snobj.compute_w('B','B','V')
   assert len(w) == 4

def test_mask_emag(snobj):
   snobj.mask_emag(0.01)
   checksum = sum([sum(snobj.data[f].mask) for f in snobj.data])
   assert checksum == 115

def test_mask_epoch(snobj):
   snobj.Tmax = 827.24
   snobj.mask_epoch(-10, 10)
   checksum = sum([sum(snobj.data[f].mask) for f in snobj.data])
   assert checksum == 99

def test_mask_SNR(snobj):
   snobj.mask_SNR(100)
   checksum = sum([sum(snobj.data[f].mask) for f in snobj.data])
   assert checksum == 115

def test_closest_band(snobj):
   assert snobj.closest_band('B', ['Bs','Vs','Rs','Is']) == 'Bs'

def test_save_load(snobj):
   import snpy
   snobj.save('temp.snpy')
   t = snpy.get_sn('temp.snpy')
   assert snobj == t


