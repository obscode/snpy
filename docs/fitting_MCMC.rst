MCMC Fitting
============

While the LM fitter can give adequate fits to most light-curves, there are
times when it would be useful to impose priors on the various parameters of
the fit. Such cases might be:

* The ``color_model`` is used on an object whose reddening is very low,
  thereby making the model insensitive to ``Rv_host``.
* The light-curve has a very poorly sampled light-curve (maybe only one
  point) and you want to use a prior on the stretch and/or time of
  maximum.
* You wish to combine the SN distance estimate with another, independent,
  estimate.
* etc.

.. _sub-MCMC-fitting
Doing the Fit
-------------

Begin fitting in the usual way using the ``fit()`` :ref:`sub-Doing-the-Fit`
function. This will provide a good starting point for the MCMC chains. Next,
use ``fitMCMC()`` to start the MCMC simulation. The calling sequence is 
much the same as ``fit()``, with some extra arguments to control the MCMC
sampling and the ability to assign priors to the model parameters.

Priors
------

When calling ``fitMCMC()``, you can assign priors to any of the model
parameters by assigning them as arguments (similar to how one kept a 
parameter fixed using ``fit()``.  There are three ways to do so:

1. Three built-in priors (Normal, Exponential, and Uniform) can be
   specified by strings. The string should begin with ``N``, ``E``,
   or ``U`` to specify Normal, Exponential, or Uniform respectively.
   Arguments for the prior then follow as a comma-separated list:
 
   * ``N,mu,sigma``: Normal distribution with mean ``mu`` and 
     standard deviation ``sigma``. Example:  ``"N,0,0.5"``.
   * ``E,tau``:  Positive exponential with scale ``tau``.
     Example: ``"E,0.7"``.
   * ``U,lower,upper``: Uniform distribution between ``lower`` and ``upper``.
     Example: ``"U,0,10"``.
2. A floating-point scalar can specified to keep the parameter fixed at
   the given value.
3. A python function that takes a single argument (the value of the parameter)
   and returns the log-probability.

Some models also have built-in priors at work. At the moment, the ``color_model`` :ref:`sub-color-model` has these kinds of priors, which are used on
the reddening law ``Rv_host``. There are currently three priors that can be 
used: ``uniform``, ``mix``, and ``bin``. These correspond to a uniform, 
Gaussian mixture model, and Normal distribution binned by color. They are
chosen by specifying the ``rvprior`` argument to ``fitMCMC()``.

Examples
--------

