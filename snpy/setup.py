#!/usr/bin/env python

import sys
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
   config = Configuration('snpy', parent_package, top_path)
   config.add_subpackage('dm15temp')
   config.add_subpackage('CSPtemp')
   config.add_subpackage('filters')
   config.add_subpackage('spline2')
   #config.add_subpackage('tspack')
   config.add_subpackage('utils')
   config.add_subpackage('sqlmod')
   config.add_subpackage('SwiftTemp')
   #config.add_scripts(['bin/snpy'])
   #config.add_scripts(['bin/update-snpy'])
   config.add_data_dir('typeIa')
   #config.add_data_files('ipythonrc-SN')
   config.add_data_files('st_calibration2.dat')
   config.add_data_files('dm15_calibration2.dat')
   #config.add_data_files('BmX_Uniform_u_all.pickle')
   config.add_data_files('color_priors.pickle')
   return config


if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(**configuration(top_path='').todict())
