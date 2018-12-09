def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration

   config = Configuration('dm15temp', parent_package, top_path)
   config.add_extension('dm15tempc',
                        sources=['dm15temp_wrap.c','dm15temp_CRB.c'],
                        libraries=['m']
                        )
   config.add_data_files('dm15temps.dat')
   config.add_data_files('tck.pickle')
   config.add_data_dir('templates')
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(version="1.0",
        description="Python Interface to J-L Prieto's template generation code", 
        author="Chris Burns",
        author_email="cburns@ociw.edu",
        **configuration(top_path='').todict())
