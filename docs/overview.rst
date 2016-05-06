What Is SNooPy?
===============

SNooPy grew organically out of a need to fit NIR light-curves of type Ia
supernovae. There was no plan, so the design may look a bit weird in places.
Currently, SNooPy is a python module for fitting and analyzing the photometry of
supernovae, particularly type Ia's.  It is meant to be run interactively through
the ``ipython`` shell, but of course, you can write non-interactive scripts as
well.  SNooPy comes with a shell script that launches ipython, loads ``snpy``
modules and makes sure the plotting driver is correctly configured. Just run it
from the shell::

  % snpy

The module, :mod:`~snpy` has several classes and utilities. Most of the interaction
comes from the use of the supernova class :class:`~snpy.sn.sn`. You instantiate a supernova
by loading data from a file using :func:`~snpy.sn.get_sn` or connecting to a 
database :mod:`~snpy.sqlmod`.
The convenience function :func:`~snpy.sn.get_sn` tries to figure out which::

   In[1]: s = get_sn('some_file.txt')   # Load ASCII file
   In[2]: s = get_sn('SN2004ef')        # Query database

The supernova is now bound to the variable ``s``. It has a a lot  of
member data (``s.ra``, ``s.decl``, ``s.z``, ``s.zcmb``) that define it's 
position, and several member functions for getting stuff done.  Most imporantly,
there is the :meth:`~snpy.sn.sn.fit` member function that fits a model to the 
data::

   In[3]: s.fit(['B','V'])

   In[4]: print "The value of dm15 is",s.dm15,"+/-",s.e_dm15
   The value of dm15 is 1.34687129103 +/- 0.0113777593656


You can plot the photometric data along with any models that you have fit,
perhaps saving it to a file::

   In[5]: s.plot(outfile="MyNiceFit.pdf")

All the usual photometric transformations (S-corrections, K-corrections, 
time-dilation) are taken care of under the hood. But you might want to plot
out the K-corrections to check they are behaving::

   In[6]: s.plot_kcorr()

You can choose different models to fit (see ...), depending on what information
you want (distance, light-curve properties, extinction properties, etc)::

   In[7]: s.choose_model('max_model')
   In[8]: s.fit()

Each parameter's best-fit value can be accessed as a member variable after
the fit. Or you can set it by hand before fitting in order to start the
fitter in a different part of parameter space::

   In[9]: print "dm15 = ",s.dm15,"+/-",s.e_dm15
   dm15 =  1.34687129103 +/- 0.0113777593656

   In[10]: s.dm15 = 1.7
   
   In[11]: s.fit()
   dm15 =  1.34793759392 +/- 0.0114967483058

It's good to save all this work you've done so you can access it later. Use the
save function, which uses python's ``pickle`` module to serialize the 
:class:`~snpy.sn.sn`
instance and write it to a file.::

   In[12]: s.save('SN2004ef_maxmodel.snpy')

You can load it back in later to using the :func:`~snpy.get_sn` function::

   In[1]: s = get_sn('SN2004ef_maxmodel.snpy')

Light-Curves
------------

Each filter has its own light-curve. Each light-curve's data is stored in
another class called :class:`~snpy.lc.lc`. It has member variables that hold
the photometric data:  ``MJD``, ``mag``, ``e_mag``, ``flux``, and ``e_flux``.
The :class:`~snpy.lc.lc` instances are stored as a dictionary member variable
of the parent :class:`~snpy.sn.sn` class called ``sn.data``::

   In[8]: print s.data
   {'B': <snpy.lc.lc instance at 0x11346e1b8>, 'g': <snpy.lc.lc instance at 0x11346e200>, 'i': <snpy.lc.lc instance at 0x11346e248>, 'H': <snpy.lc.lc instance at 0x11346e290>, 'K': <snpy.lc.lc instance at 0x11346e2d8>, 'J': <snpy.lc.lc instance at 0x11346e320>, 'r': <snpy.lc.lc instance at 0x11346e368>, 'u': <snpy.lc.lc instance at 0x11346e3b0>, 'V': <snpy.lc.lc instance at 0x11346e3f8>, 'Y': <snpy.lc.lc instance at 0x11346e440>}

You can also refer to each filter by name as a member variable of the parent
:class:`~snpy.sn.sn` class::

   In[9]: print "B-band data covers MJD=", s.B.MJD.min(), "to MJD=", s.B.MJD.max()
   B-band data covers MJD= 53255.17803 to MJD= 53329.09918

You can plot just the single filter's data (along with residuals if a model
has been fit for that filter)::

   In[10]: s.B.plot()

Or maybe you don't want to fit a SNooPy template, but would prefer to fit a
spline to the light-curve and derive some useful quantities::

   In[11]: s.B.spline_fit(method='spline2')
   In[12]: print "B peaks at",s.B.Tmax," and has dm15=",s.B.dm15
   B peaks at 53263.6826124  and has dm15= 1.41942033475

Actually, :meth:`~snpy.sn.sn.spline_fit` is a bit of a misnomer, as you can 
fit with many different kinds of interpolators (spline, polynomial, Guassian 
process, etc).  See (....).

Other Useful Bits
-----------------

If you have astropy installed in your python environment, you can get
the Hubble flow value of the distance modulus (using default LambdaCDM
cosmology)::

   In[13]: s.choose_model('EBV_model')
   In[14]: s.fit()
   In[15]: print "Ia distance:",s.DM,"+/-",s.e_DM,"\nHubble distance",s.get_distmod().value
   Ia distance: 35.359938857 +/- 0.014760098834
   Hubble distance 35.4486366281

SNooPy has several spectral energy distrubutions for type Ia SNe built in.
You can access them through the :mod:`snpy.kcorr` sub-module (that's where 
they are used most heavily). This will retrieve the [Hsiao+2007]_ SED at 
maximum::

   wave,flux = kcorr.get_SED(day=0, version='H3')

SNooPy also has a module for dealing with the photometric properties of 
filters called :mod:`~snpy.filters`. Each :class:`~snpy.lc.lc` instance has an 
instance of the :class:`~snpy.filters.filter` filter it represents::

  In[16]: print s.B.filter
  B:  TAM scanned B filter for Swope at LCO + CTIO extinction

  In[17]: print "The Hsiao SED template as B=",s.B.filter.synth_mag(wave,flux),"at maximum"
  he Hsiao SED template as B= 0.0211822721773 at maximum

The :meth:`~snpy.sn.sn.spline_fit` function is actually a wrapper around a
useful unified interface to interpolating 1-D data:
:mod:`~snpy.utils.fit1dcurve`.  After calling the 
:meth:`~snpy.sn.sn.spline_fit` function, an instance of the
:class:`~snpy.utils.fit1dcurve.oneDcurve` class is available as a member 
variable of the :class:`~snpy.lc.lc` instace. It can be used to further 
analyze the data::

   In[18]: s.B.spline_fit(method='spline2')

   In[19]: t,m,c = s.B.interp.find_extrema()

   In[20]: for i in range(len(t)):
     ....:     print ["maximum ","minimum "][c[i]<0]+"found at (%.2f,%.2f)" % (t[i],m[i])

   maximum found at (53260.78,17.32)
   minimum found at (53264.45,17.37)
   maximum found at (53265.80,17.37)
   minimum found at (53275.88,18.04)
   maximum found at (53286.39,17.78)
   minimum found at (53310.60,19.21)
   maximum found at (53311.57,19.21)
   minimum found at (53319.67,19.61)
   maximum found at (53323.81,19.58)

Peruse the rest of the documentation for more in-depth explanations of the
API and all the features SNooPy has. 

Getting Help
------------

Python has an internal help system which utilizes comments at the
beginning of functions and classes (so-called docstrings). Simply
use the built-in help() function to get help on an item. Here are
some examples (output is not shown to save space)::

   In [1] help(sn)
   In [2] help(sn.fit)
   In [3] help(sn.plot)
   In [4] help(lc)

Line 1 gets help about the entire :class:`~snpy.sn.sn` class, which will list
all the functions defined therein, including internal ones that are
not meant to be used by end users (but of course are available to
be hacked, but may lack good documentation). Lines 2 and 3 get more
specific help on individual member functions. Line 4 gets help on
the :class:`~snpy.lc.lc` (light-curve) class. You can ask for help on any python
object (including variables).

