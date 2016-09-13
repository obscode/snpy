def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration

   config = Configuration('filters', parent_package, top_path)
   config.add_data_dir('filters')
   config.add_data_dir('standards')
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(version="1.0",
        description="A python module for dealing with (mostly SN related) filter systmes", 
        author="Chris Burns",
        author_email="cburns@ociw.edu",
        **configuration(top_path='').todict())
