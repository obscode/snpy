#!/usr/bin/env python
def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration

   config = Configuration('utils', parent_package, top_path)
   config.add_data_files('Ia_R_splines.pickle')
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(version="1.0",
        description="Bunch of misc. utilities.", 
        author="Chris Burns",
        author_email="cburns@ociw.edu",
        **configuration(top_path='').todict())
