#!/usr/bin/env python

# This is the numpy.distutils verison of setup.py.

import os,sys, string
from pkg_resources import parse_version
try:
   import numpy
except ImportError:
   print("Error:  setup.py requires numpy")
   sys.exit(1)
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from numpy.distutils.system_info import get_info
try:
   import matplotlib
   have_mpl = 1
except:
   have_mpl = 0

numpy_min_version = '1.7'

def get_numpy_status():
   """
      Returns a dictionary containing a boolean specifying whether NumPy
      is up-to-date, along with the version string (empty string if
      not installed).
   """
   numpy_status = {}
   try:
      import numpy
      numpy_version = numpy.__version__
      numpy_status['up_to_date'] = parse_version(
         numpy_version) >= parse_version(numpy_min_version)
      numpy_status['version'] = numpy_version
   except ImportError:
      numpy_status['up_to_date'] = False
      numpy_status['version'] = ""
   return numpy_status


def configuration(parent_package='', top_path=None):
   config = Configuration(None, parent_package, top_path)
   config.add_subpackage('snpy')
   config.add_scripts(['bin/snpy'])
   config.add_scripts(['bin/update-snpy'])
   return config


def dosetup():
   status = get_numpy_status()
   if not status['up_to_date']:
      raise ImportError("You need at least version {} to run SNooPy"
                        .format(numpy_min_version))


   with open(os.path.join('snpy','version.py'), 'r') as fi:
      line = fi.readline()
      v = line.split('=')[1].strip()
      version = v[1:-1]

   setup(version=version,
         name='snpy',
         description="SNooPy:  Supernova light-curve analysis tool",
         author='Chris Burns (Carnegie Observatories)',
         author_email='cburns@carnegiescience.edu',
         url='http://csp.obs.carnegiescience.edu/data/snpy',
         license='MIT',
         classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Framework :: IPython',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7',
            'Topic :: Scientific/Engineering :: Astronomy'],
         setup_requires=['numpy','pytest-runner'],
         tests_require=['pytest'],
         install_requires=[
            'numpy',
            'scipy',
            'pymysql',
            'matplotlib',
            'ipython',
            'astropy',
            'certifi',
            'pandas',
            'scikit-learn'],
         **configuration(top_path='').todict())

if __name__ == '__main__':
   dosetup()
