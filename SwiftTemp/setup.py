import sys
import os
from glob import glob


def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration

   config = Configuration('SwiftTemp', parent_package, top_path)
   config.add_data_files('SNIa_m2_template_fe.dat','SNIa_w1_template_fe.dat',
         'SNIa_w2_template_fe.dat')
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(version="0.1",
        description="Python Interface to Swift template generator",
        author="Chris Burns",
        author_email="cburns@ociw.edu",
        url="http://www.ociw.edu/~cburns",
        **configuration(top_path='').todict())

