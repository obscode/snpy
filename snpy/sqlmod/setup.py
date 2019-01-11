def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration

   config = Configuration('sqlmod', parent_package, top_path)
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(version="1.0",
        description="A system for importing the correct SQL access functions.", 
        author="Chris Burns",
        author_email="cburns@ociw.edu",
        **configuration(top_path='').todict())
