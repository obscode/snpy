import unittest
import snpy
import os

class TestSNooPy(unittest.TestCase):

   def setUp(self):
      # First, we make sure we have the test data
      self.base = os.path.dirname(snpy.__file__)
      self.data = os.path.join(self.base, 'tests')
      pass

   def test_loadTxt(self):
      s = snpy.get_sn(os.path.join(self.data, 'SN2006ax.txt'))
      self.assertEqual(s.name, "SN2006ax")

if __name__ == '__main__':
   unittest.main()

