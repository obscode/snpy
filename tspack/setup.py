#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration
   config = Configuration('tspack', parent_package, top_path)
   config.add_extension('tspack',
                        sources=['tspack.pyf','tspack.f'])
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(description='Python wrapper to the Tension spline by Robert J. Renka',
      author='Chris Burns',
      author_email='cburns@ociw.edu',
      **configuration(top_path='').todict())
