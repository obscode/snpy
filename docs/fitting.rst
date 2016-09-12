Fitting Light-curves
====================

All the light-curve fitting is done through a member function of the
:class:`~snpy.sn.sn` class: :meth:`~snpy.sn.sn.fit`.You basically tell ``fit()``
which filters to fit (or all by default) and which parameters to hold
constant. This function sets up some initial book-keeping, and then
hands off the work to a :class:`~snpy.model.model` instance. The ``model``
class defines the model that will be used to fit the light-curves.
I did it this way so that adding new models (or modifying existing
ones) would be easier. SNooPy comes with four built-in models: 
:ref:`sub-ebv_model`, :ref:`sub-ebv_model2`, :ref:`sub-color_model`,
and :ref:`sub-max_model`, which
I will explain in the next sections. At any time, you can switch
between the two by using the :meth:`~snpy.sn.sn.choose_model()` member function
of the ``sn`` class.
As of SNooPy version 2, you can fit the ``EBV_model2`` and ``max_model``
using either :math:`\Delta m_{15}` or a new stretch-like parameter :math:`s_{BV}` which
is introduced in the CSP's most recent analysis paper [Burns2014]_. The choice
between the two is made when your select the model using 
:meth:`~snpy.sn.sn.choose_model` and specify the ``stype`` optional
argument to be ``'dm15'`` or ``'st'``.

With the advent of increased complexity of the models and the need
for adding priors to parameters, a new experimental version of the
``fit()`` function is now available: fitMCMC(). Unlike the original,
which uses the Levenberg-Marquardt least-squares fitter, fitMCMC uses
a Markov Chain Monte Carlo fitter called ``emcee``. The advantage
of this is that MCMC, being a Bayesian framework, allows you to specify
priors on your parameters very easily. This can be particularly useful
if the model you are trying to fit is insenstive to one or more of
its parameters. This can happen when your light-curve has very few
points (so that shape is not well defined) or you are using the ``color_model``
at low ``E(B-V)``. Before, you would simply have had to keep
those parameters fixed a reasonable values. Now you can specify simple
priors so that the uncertainties are far more realistic.

.. toctree::
   :maxdepth: 2

   models
   fitting_LM
   fitting_MCMC
   photometric_systems

.. [Burns2014] Burns et al. ApJ, 189, 32B (2014). `ADS <http://adsabs.harvard.edu/abs/2014arXiv1405.3934B>`_
