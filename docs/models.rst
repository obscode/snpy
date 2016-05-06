Models
======

By design, the specification of models is done in its own namespace: ``model``.
All ``model`` classes are subclasses of ``model.model``. In the following
sections, I detail what each model can be used for and the mathematics behind
it.

``SNooPy`` now has two different parametrizations for the shape of the
light-curves. One is the traditional :math:`\Delta m_{15}` of
[Phillips 1999]_ and the newer 'color-stretch' introduced in 
[Burns+2011]_, which more consistently deals with fast-declining events. Note
that the ``EBV_model`` will only work with :math:`\Delta m_{15}`-based
templates.

.. _sub-ebv_model:

EBV_model
---------

This light-curve model is a variation of that given by [Prieto+2006]_. It
is most useful for determining a distance using the Phillips relation
and correcting for dust extinction in the host.
We compare the observed data in :math:`N` filters to the following model:

.. math::
   :nowrap:

   \begin{eqnarray*}
      m_{X}\left(t\right) & = & T_{Y}\left(t^{\prime},\Delta m_{15}\right)+M_{Y}\left(\Delta m_{15}\right)+\mu+R_{X}E\left(B-V\right)_{gal}+\\
                                  &   & R_{Y}E\left(B-V\right)_{host}+K_{X,Y}
   \end{eqnarray*}

where :math:`m_{X}` is the observed magnitude in band :math:`X`, :math:`t`
is the observed time relative to B maximum, :math:`t^\prime` is the
de-redshifted time relative to B maximum, :math:`\Delta m_{15}` is the decline
rate parameter [Phillips1999]_, :math:`M_{Y}` is the absolute
magnitude in filter :math:`Y` in the rest-frame of the supernova,
:math:`E(B-V)_{gal}` and :math:`E(B-V)_{host}` are the reddening due to
galactic foreground and host galaxy, respectively, :math:`R_{X}` and
:math:`R_{Y}` are the total-to-selective absorptions for filters :math:`X` and
:math:`Y`, respectively, :math:`K_{XY}` is the cross-band k-correction from
rest-frame :math:`X` to observed filter :math:`Y`. Note that the k-corrections
depend on the epoch and *can* depend on the host and galaxy extinction (as
these modify the shape of the spectral template). In the fitting, one has the
choice of the spectral template of [Nugent+2002]_ or [Hsiao+2007]_.  The latter
is the default. You can choose to allow the K-corrections to vary during the
fit, keep them fixed, or not use them at all. See the parameters for the
:meth:`~snpy.sn.sn.fit` function.

The template :math:`T\left(t,\Delta m_{15}\right)` can be generated from
the code of  \citet{2006ApJ...647..501P} or \citet{Burns2011}. In
the former case, you will be fitting to rest-frame :math:`B_{s}V_{s}R_{s}I_{s}`
[#f1]_ while in the latter case, you can fit to :math:`uBVgriYJH` from the CSP
data [Contreras+2010]_. You can mix and match which template
generator you use: it is all a matter of which filter you choose for
the ``restbands`` instance variable (see :ref:`sub-Doing-the-Fit`).
In all, this model fits 4 parameters: ``Tmax``, ``dm15``,
``EBVhost``, and ``DM``.

Note that to determine the host reddening, you need to fit at least
two distinct *restframe* filters. For now, I've left the galactic
and host reddening laws (:math:`R_{X}` and :math:`R_{Y}`) as member variables
rather than parameters to be fit
[#f2]_. This could change in the future.


.. _sub-ebv_model2:

EBV_model2
----------

This is the same model as EBV_model, except that it can fit both
:math:`\Delta m_{15}`- and :math:`s_{BV}`-based templates. And instead of the calibration
presented in [Folatelli+2010]_, the calibration from the
upcoming CSP analysis paper [Burns+2011]_ is used. When using
the ``choose_model`` function to select this model, you can
specify ``stype='dm15'`` or ``stype='st'`` to select the
two different light-curve parameters.

.. _sub-max_model:

max_model
---------

An alternative to assuming some kind of relationship between the different
filters (i.e., some kind of reddening law as is done in the ``EBV_model``),
one can simply fit the maximum magnitude of each light-curve. This
is how the ``max_model`` model works. It uses the same light-curve
templates as ``EBV_model``, but instead fits the following model:

.. math::
   :nowrap:

   \[
   m_{X}\left(t\right) = T_{Y}\left(t^\prime,\Delta m_{15}\right)+m_{Y}+R_{X}E\left(B-V\right)_{gal}+K_{X,Y}
   \]

where :math:`m_{Y}`, the peak magnitude in filter :math:`Y` is now a parameter.
Of course, depending on the number of filters you fit, you will have
that many :math:`m_{Y}` parameters. Note that K-corrections and Milky-way
extinction are performed exactly as in the ``EBV_model`` case.
Therefore, for :math:`N` filters, max_model will have :math:`N+2` parameters:
``Tmax``, ``dm15``, and :math:`N` \emph{f}``max``, where
\emph{f} is the rest-band filter name (so if, for instance, you fit
an observed light-curve to restframe :math:`B`, there will be a ``Bmax``
parameter). As with the ``EBV_model2`` model, you can choose
which light-curve parameter you wish to use by specifying the ``stype``
argument to ``choose_model()``.

**NOTE:** Please be aware that any non-linear fitter will only
find a \emph{local} minimum in :math:`\chi^{2}`. It is up to you, the user,
to try different starting points in parameter space and see which
one gives the overall best-fit. After each fit, the member variable
``rchisq`` is updated with the reduced-:math:`\chi^{2}`. You can therefore
use this to discriminate between different solutions.


.. _sub-color_model:

color_model
-----------

In this model, intrinsic colors from [Burns+2014]_ are used to
infer the amout of extinction :math:`E(B-V)` as well as the shape reddening
law :math:`R_{V}`. Mathematically, the model being fit to the light-curves
is:

.. math::
   :nowrap:

   \begin{eqnarray*}
      m_{X}\left(t\right) & = & T_{Y}\left(t^{\prime},\Delta m_{15}\right)+B_{max}+(X-B)\left(s_{BV}\right)+R_{X}E\left(B-V\right)_{gal}+R_{Y}\left(R_{V}\right)E\left(B-V\right)_{host} + \\
                                  &   & K_{X,Y}
   \end{eqnarray*}

where :math:`B_{max}` is the de-reddened and K-corrected :math:`B` maximum (treated
as a free parameter) and :math:`(X-B)\left(s_{BV}\right)` is the intrinsic
:math:`X-B` color, which is a function of :math:`s_{BV}`. In [Burns+2014]_
it is modeled as a 2nd degree polynomial in :math:`\left(s_{BV}-1\right)`.
All other variables have the same meaning as in previous models. The
model has 5 free parameters: :math:`s_{BV}`, :math:`T_{max}`, :math:`B_{max}`, :math:`E(B-V)`,
and :math:`R_{V}`. Note that the distance modulus is not included in this
model.

A major complication of this model is that two parameters, :math:`R_{V}`
and :math:`E(B-V)`, appear as multiplicative factors in the same term.
At best, they will be highly covariant. At worst (for low values of
:math:`E(B-V)`), the model becomse insensitive to :math:`R_{V}`. For this reason,
it may be necessary to impose priors on :math:`R_{V}`. Two such priors
were introduced in [Burns+2014]_: a Gaussian mixture model that
applies to all SNeIa, and a binned prior, where a separate Gaussian
prior is applied to the SNIa depending on its value of :math:`E(B-V)`.
Because SNooPy uses the LM least-squares algorithm by default, there
is no natural way to incorporate these priors using the standard ``fit()``
routine. Instead, use the ``fitMCMC()`` routine, which fits the
light-curves using Markov Chain Monte Carlo and allows priors to be
specified on parameters. In this case, use the keyword argument ``rvprior='mix'``
for the Gaussian mixture model, or ``rvprior='bin'`` for the
binned prior. See section :ref:`sub-MCMC-fitting` for more details.

.. [Phillips1999] Phillips, M.M., AJ, 118, 1766 (1999)
   http://adsabs.harvard.edu/abs/1999AJ....118.1766P
.. [Nugent+2002] Nugent et al., PASP, 114, 803 (2002). 
   http://adsabs.harvard.edu/abs/2002PASP..114..803N
.. [Hsiao+2007] Hsiao et al., ApJ, 663, 1187 (2007).
   http://adsabs.harvard.edu/abs/2007ApJ...663.1187H
.. [Stritzinger+2005] Stritzinger et al., PASP, 117, 810 (2005)
   http://adsabs.harvard.edu/abs/2005PASP..117..810S
.. [Contreras+2010] Contreras et al., AJ, 139, 519 (2010).
   http://adsabs.harvard.edu/abs/2010AJ....139..519C
.. [Burns+2011] Burns et al., AJ, 141, 19B (2011).
   http://adsabs.harvard.edu/abs/2011AJ....141...19B
.. [Burns+2014] Burns et al., ApJ, 789, 32B (2014).
   http://adsabs.harvard.edu/abs/2014ApJ...789...32B

.. [#f1] The 's' subscript refers to 'standard', which is to say the Bessel
   filters from [Stritzinger+2005]_
.. [#f2] So far, data I've analyzed hasn't been good enough to distinguish
   between reddening laws.
