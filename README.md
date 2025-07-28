PURPOSE
-------

SNPY (or more affectionately SNooPy) is a Python package for fitting the
light-curves of TypeIa supernovae.  It is NOT an off-the-shelf program that
will do all the work for you from the command line.  It is meant to be run
interactively.  It provides many tools, but does not lock you into a specific
routine for fitting and therefore allows for a great deal of customization.
The package also contains a number of stand-alone sub-packages that can be
installed seperately.

You can find intructions for installing snpy in the INSTALL file and at
the following website:

http://csp.obs.carnegiescience.edu/data/snpy/installing_snoopy2

Documentation can be found the the docs/ directory under the source tree.  
But briefly, here is what SNOOPY can do:

- Fit lightcurve templates to data.  These templates were constructed
  from the Carnegie Supernova Project's low-z data set (as outlined in
  [Contreras et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009AAS...21442704C/abstract) and [Folatelli et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010AJ....139..120F/abstract)).  These include
  templates for the CSP uBVgriYJH filter set.  You can also use Jose-
  Louis Prieto's templates to fit Johnson/Kron/Cousins BRVI lightcurves.

- Using templates, fit for distance (and dm15, E(B-V), stretch, time of 
  maximum, etc.  This uses the calibration of [Folatelli et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010AJ....139..120F/abstract) 
  for the CSP filter set, or [Prieto et al. (2006)](https://ui.adsabs.harvard.edu/abs/2006ApJ...647..501P/abstract) for the BVRI set.

- Alternatively, fit for light-curve parameters only (dm15, stretch,
  time of maximum, maximum light, colors, etc).

- Compute k-corrections based on the spectral energy distribution
  template of [Hsiao et al. (2007)](https://ui.adsabs.harvard.edu/abs/2007ApJ...663.1187H/abstract) (+ [Lu et al. 2003](https://ui.adsabs.harvard.edu/abs/2023ApJ...948...27L/abstract) if would like to use 
  the updated NIR template).  This includes color-matching the
  SED template to the observed filters.  You can do plain K-corrections
  or cross-band K-corrections.  This is normally done as part of the 
  fitting procedure, but can be done separately.

- Fit splines, Gaussian Processes, and polynomials to light-curves to
  do your own analysis, independent of the light-curve templates.

- Built-in Lira Law to estimate E(B-V).

- Various plotting routines.

- Interact with a properly set-up SQL database.

The subpackages include:

- filters:   a package for working with filter systems.  computing
  zero-points, producing synthetic photometry, etc.

- dm15temp:  a python wrapper to Jose-Louis Prieto's light-curve
  generator

- CSPtemp:  a generator for the CSP light-curve templates.

- spline2:  a python wrapper to the Hyperspline (or spline2) algorithm,
  which uses the Durbin-Watson statistic rather than chi-square.  Good
  if you don't trust your error bars.

- tspack:  a python wrapper to the TSPACK library:  generating tension
  splines.'''

