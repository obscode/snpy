#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration
   config = Configuration('spline2', parent_package, top_path)
   config.add_extension('spline2c',
                        sources=['spline2_nonint.c','spline2_wrap.c']
                        )
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup( version="1.0",
      description="Python interface spline2",
      author="Chris Burns",
      author_email="cburns@ociw.edu",
      **configuration(top_path='').todict())
