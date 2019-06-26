#!/usr/bin/env python
#from distutils.core import setup, Extension
import os
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from numpy.distutils.system_info import get_info

# see if we are dealing with a known os:
if os.name != "posix":
   raise Exception, "Sorry, only POSIX operating systems are supported"

lapack_info = get_info('lapack_opt', 1)
config = Configuration('mydistance', parent_package=None, top_path=None)

config.add_extension(name='aniso_euclidean_s',
      sources = ['src/mydistances.f'], extra_info=lapack_info)

setup(version="1.0",
      description="Python interface custom anisotropic distance",
      author="Chris Burns",
      author_email="cburns@carnegiescience.edu",
      packages=['mydistance'],
      **(config.todict()))
