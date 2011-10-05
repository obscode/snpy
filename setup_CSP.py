#!/usr/bin/env python
#from numpy.distutils.misc_util import Configuration
#from numpy.distutils.system_info import get_info
import os,sys, string
try:
   import matplotlib
   have_mpl = 1
except:
   have_mpl = 0

# Here are the packages we need to run SNPY
required_packages = [\
      ['scipy','http://sourceforge.net/projects/scipy/files/scipy/0.7.1/scipy-0.7.1.tar.gz/download'],
      ['IPython', 'http://ipython.scipy.org/dist/0.10/ipython-0.10.tar.gz'],
      ['pyfits','http://www.stsci.edu/resources/software_hardware/pyfits/pyfits-2.1.1.tar.gz']]

# These are optional for added functionality.
optional_packages = [\
      ['MySQLdb', 'http://sourceforge.net/projects/mysql-python/files/mysql-python-test/1.2.3c1/MySQL-python-1.2.3c1.tar.gz/download']]

# Now, if we have matplotlib installed, then the PGPLOT wrappers are simply optional.
#  PGPLOT is a pain to install, but some may prefer it's speed.
if have_mpl:
   optional_packages.append(['pygplot','http://www.ociw.edu/Code/python/burns-python-scripts/pygplot/pygplot-releases/pygplot-0.97.tar.gz'])
else:
   required_packages.append(['pygplot','http://www.ociw.edu/Code/python/burns-python-scripts/pygplot/pygplot-releases/pygplot-0.97.tar.gz'])

def get_package(package, url, setup=1, test=1, verbose=1, options=''):
   '''Retrieve a python package from the source url, unpack it, 
   and install it'''
   from urllib import urlopen
   import tarfile
   from StringIO import StringIO
   import tempfile

   if verbose:  print "Downloading %s..." % (package) 
   u = urlopen(url)
   t = StringIO(u.read())
   if not setup:
      f = open(os.path.basename(url),'w')
      f.write(u.read())
      return 1,"Package saved as %s" % (os.path.basename(url))

   if verbose:  print "Extracting..."
   # extract the package and run the old python setup install
   try:
      tar = tarfile.open('blah',fileobj=t)
   except:
      return 0,"Sorry, URL %s does not have a valid TAR file, package %s install failed" % (package, url)
   tmp = tempfile.mkdtemp()
   cwd = os.getcwd()   # remember where we started
   os.chdir(tmp)
   # instead of extractall(), we'll use the old way
   #tar.extractall()
   names = tar.getmembers()
   for name in names:  tar.extract(name)
   files = os.listdir('.')
   if len(files) == 1:
      os.chdir(os.path.join(tmp,files[0]))

   if verbose:  
      print "Running python setup.py install %s (this may take a while)..." %\
            (options)
   logfile = os.path.join(cwd, package+".log")
   res = os.system('%s setup.py install %s > %s 2> %s' % \
         (sys.executable, options, logfile, logfile))
   if res != 0:
      return 0,"Sorry, package %s failed to install, see %s.log" % (package,package)
   os.chdir(cwd)
   if verbose:  print "Cleaning up..."
   os.system('rm -rf %s' % tmp)
   
   if not test:
      return 1, "Package %s setup succesfully" % (package)
   
   try:
      res = __import__(package)
   except ImportError:
      return 0, "Package %s was installed, but failed to load" % (package)

   return 1, "Package %s installed and loaded properly" % (package)
   

def configuration(parent_package='', top_path=None):
   config = Configuration('snpy', parent_package, top_path)
   config.add_subpackage('dm15temp')
   config.add_subpackage('CSPtemp')
   config.add_subpackage('filters')
   config.add_subpackage('spline2')
   #config.add_subpackage('dust_getval')
   config.add_subpackage('tspack')
   config.add_subpackage('utils')
   #config.add_subpackage('sqlmod')
   config.add_scripts(['bin/snpy'])
   config.add_data_dir('typeIa')
   config.add_data_files('ipythonrc-SN')
   return config


if __name__ == '__main__':
   # First, let's see if the user has any extra options
   extra_options = string.join(sys.argv[2:])
   try:
      from numpy.distutils.core import setup
   except ImportError:
      print "You don't have numpy installed.  Would you like me to try to"
      print "install it?  (y/n)"
      answer = ""
      while answer != 'n' and answer != 'y':
         answer = raw_input()
      if answer == 'n':
         print "Okay, but you'll need to install Numpy before installing SNPY"
         sys.exit(1)
      
      # So, at this point, we're going to try to install numpy:
      flag,mesg = get_package('numpy',
            'http://sourceforge.net/projects/numpy/files/NumPy/1.2.1/numpy-1.2.1.tar.gz/download',
            options=extra_options)
      if flag != 1:
         print mesg
         sys.exit(1)
   from numpy.distutils.core import setup
   from numpy.distutils.misc_util import Configuration
   from numpy.distutils.system_info import get_info

   # Now, make sure we have all the other stuff too:
   for package,url in required_packages:
      try:
         res = __import__(package)
      except ImportError:
         if package == "pyfits":
            # FITS will also do... so see if that's installed
            try:
               res = __import__('FITS')
               continue
            except ImportError:
               pass


         print "Package %s does not appear to be installed.  Would you" % \
               (package)
         print "like me to install it now?  (y/n)"
         answer = ""
         while answer != 'y' and answer != 'n':
            answer = raw_input()
         if answer == 'n':
            print "Okay, but you'll need to install %s before installing SNPY"\
                  % (package)
            print "You can find it here:"
            print url
            sys.exit(1)

         flag,mesg = get_package(package, url, options=extra_options)
         if flag != 1:
            print mesg
            sys.exit(1)
         print mesg

   # Now, see if optional packages are there
   for package,url in optional_packages:
      try:
         res = __import__(package)
      except ImportError:
         print "Package %s is does not appear to be installed.  It is optional" %\
               (package)
         print "Would you like me to try installing it now?  (y/n)"
         answer = ""
         while answer != 'y' and answer != 'n':
            answer = raw_input()
         if answer == 'n':
            continue
         flag,mesg = get_package(package, url, options=extra_options)
         print mesg
   setup(version='0.5.0',
         author='Chris Burns (Carnegie Observatories)',
         author_email='cburns@ociw.edu',
         **configuration(top_path='').todict())
