import sys
import os
from glob import glob

def get_gsl_info(path=''):
   '''Look in the standard places for the GSL library.'''
   inc_dirs = []
   lib_dirs = []
   libraries = []
   p = os.popen('%s --libs' % (os.path.join(path, 'gsl-config')))
   libs = p.readline();  p.close()
   if libs == '':
      print "Warning.  Could not run gsl-config!"
      return None,None,None
   libs = libs.split()
   for item in libs:
      if item[0:2] == '-l':  libraries.append(item[2:])
      if item[0:2] == '-L':  lib_dirs.append(item[2:])
   p = os.popen('%s --cflags' % (os.path.join(path, 'gsl-config')))
   cflags = p.readline();  p.close()
   cflags = cflags.split()
   for item in cflags:
      if item[0:2] == '-I':
         inc_dirs.append(item[2:])

   return(inc_dirs, lib_dirs, libraries)

def get_temp_GSL():
   import urllib, tarfile, StringIO,tempfile
   global tdir
   u = urllib.urlopen('ftp://ftp.gnu.org/gnu/gsl/gsl-1.8.tar.gz')
   f = StringIO.StringIO(u.read())
   tar = tarfile.open(fileobj=f)
   tdir = tempfile.mkdtemp()
   cwd = os.pwd()
   os.chdir(tdir)
   tar.extractall()
   os.chdir(os.listdir()[0])
   tar.close()
   u.close()
   stat = os.system('./configure --prefix=%s --disable-shared' % tdir)
   if stat != 0:
      return 0
   stat = os.system('make')
   if stat != 0:
      return 0
   stat = os.system('make install')
   if stat != 0:
      return 0

   return(os.path.join(tdir,'include'), os.path.join(tdir,'lib'))


def configuration(parent_package='', top_path=None):
   from numpy.distutils.misc_util import Configuration
   global tdir
   incdirs,libdirs,libs = get_gsl_info()

   tdir = None    # if we need a temporary dir, this will be it

   if incdirs is None:
      print "I couldn't run the GNU Scientific Library (GSL) gsl-config script."
      print "Would you like me to download/install it temporarily? (y/n)"
      answer = ""
      while answer not in ['y','n']:
         answer = raw_input()
      if answer == 'n':
         print "Okay, you'll have to install it yourself, then.  When you've"
         print "done so, please make sure the gsl-config script is in your"
         print "PATH."
         sys.exit(1)
      else:
         stat = get_temp_GSL()
         if stat == 0:
            print "Sorry, the automatic download/build failed.  You're going"
            print "to have to install GSL yourself."
            if tdir is not None:  os.system('rm -rf %s' % tdir)
            sys.exit(1)
         incdirs,libdirs,libs = get_gsl_info(path=os.path.join(tdir,'bin'))

   config = Configuration('CSPtemp', parent_package, top_path)
   config.add_extension('dm15temp2c',
                       sources=['dm15temp_wrap.c','gloes.c','dm15temp.c'],
                       libraries=libs,
                       include_dirs=incdirs,
                       library_dirs=libdirs)
   config.add_data_files('templates.dat','tck.pickle','bs_error.pickle')
   config.add_data_dir('templates')
   config.add_data_dir('fits')
   config.add_scripts('generate_surf.py')
   return config

if __name__ == '__main__':
   from numpy.distutils.core import setup
   setup(version="0.1",
        description="Python Interface to CSP template generator",
        author="Chris Burns",
        author_email="cburns@ociw.edu",
        url="http://www.ociw.edu/~cburns",
        **configuration(top_path='').todict())
   if tdir is not None:  os.system('rm -rf %s' % tdir)

