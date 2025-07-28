#!/usr/bin/env python

# numpy.distutils is now defunct, trying to use the modern version of
# setuptools. So setup.py is super-simple and all meta data is in
# pyproject.toml
# EXCEPT:  to build C-extensions, seems we have to do that here...?

from setuptools import setup, Extension
import numpy
setup(
      ext_modules=[
         Extension(name="snpy.spline2.spline2c",
                   sources=["snpy/spline2/spline2_nonint.c",
                            "snpy/spline2/spline2_wrap.c"],
                   include_dirs = [numpy.get_include()],
                   )
         ],
      scripts = ["bin/snpy","bin/update-snpy"],
      )
