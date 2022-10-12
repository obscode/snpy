import os


def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration

   config = Configuration('ATLAStemp', parent_package, top_path)
   config.add_data_dir('fits')
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(version="0.1",
        description="Python Interface to ATLAS template generator",
        author="Chris Burns",
        author_email="cburns@ociw.edu",
        url="http://www.ociw.edu/~cburns",
        **configuration(top_path='').todict())
   if tdir is not None:  os.system('rm -rf %s' % tdir)

