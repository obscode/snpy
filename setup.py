#!/usr/bin/env python
import os,sys, string
try:
   import matplotlib
   have_mpl = 1
except:
   have_mpl = 0

def configuration(parent_package='', top_path=None):
   config = Configuration('snpy', parent_package, top_path)
   config.add_subpackage('dm15temp')
   config.add_subpackage('CSPtemp')
   config.add_subpackage('filters')
   config.add_subpackage('spline2')
   #config.add_subpackage('dust_getval')
   config.add_subpackage('tspack')
   config.add_subpackage('utils')
   config.add_subpackage('sqlmod')
   config.add_subpackage('SwiftTemp')
   config.add_scripts(['bin/snpy'])
   config.add_scripts(['bin/update-snpy'])
   config.add_data_dir('typeIa')
   config.add_data_dir('tests')
   config.add_data_files('ipythonrc-SN')
   config.add_data_files('st_calibration2.dat')
   config.add_data_files('dm15_calibration2.dat')
   config.add_data_files('BmX_Uniform_u_all.pickle')
   config.add_data_files('color_priors.pickle')
   return config


if __name__ == '__main__':
   # First, let's see if the user has any extra options
   extra_options = string.join(sys.argv[2:])
   try:
      from numpy.distutils.core import setup
   except ImportError:
      print "You don't have numpy installed.  Would you like me to try to"
      print "install it?  (y/n)"
      answer = ""
      while answer != 'n' and answer != 'y':
         answer = raw_input()
      if answer == 'n':
         print "Okay, but you'll need to install Numpy before installing SNPY"
         sys.exit(1)
      
      # So, at this point, we're going to try to install numpy:
      flag,mesg = get_package('numpy',
            'http://sourceforge.net/projects/numpy/files/NumPy/1.2.1/numpy-1.2.1.tar.gz/download',
            options=extra_options)
      if flag != 1:
         print mesg
         sys.exit(1)
   from numpy.distutils.core import setup
   from numpy.distutils.misc_util import Configuration
   from numpy.distutils.system_info import get_info

   setup(version='2.0b',
         author='Chris Burns (Carnegie Observatories)',
         author_email='cburns@obs.carnegiescience.edu',
         **configuration(top_path='').todict())
