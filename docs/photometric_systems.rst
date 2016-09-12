.. _sub-Photometric-systems:

Photometric Systems
===================

The templates that SNooPy uses to fit light-curves were trained on the CSP 
natural photometric system (see :ref:`natural_v_standard` below). In order to
use them, therefore, one needs to convert from whatever system you measured in
to the CSP system. This is usually referred to as performing an S-correction
[Stritzinger+2005]_.

SNooPy will compute the S-corrections for you as long as you completely
specify the photometric system in which your data was obtained. This usually
means specifying:

   * the shape of the filter+camera+telescope+atmostphere transmission function
   * the zero-point of the filter:  what :math:`m=0` means, or:
   * a primary standard SED and magnitude of the standard.

In practice, you define all your filters in SNooPy's :mod:`~snpy.filters` module
using a unique ID and then load your data into SNooPy using these IDs.

Specifying your Photometric System
----------------------------------

Internally, SNooPy relates magnitudes to *photon* fluxes through the following
formula:

.. math::
   :nowrap:

   \begin{equation*}
      m_X = -2.5*\log_{10}\left[\frac{1}{ch}\int f_\lambda\left(\lambda\right)
                    S_X\left(\lambda\right) \lambda d\lambda\right] + ZP_X
   \end{equation*}
                
where :math:`m_X` is the magnitude for filter :math:`X`, :math:`f_\lambda` is
the flux in :math:`erg\cdot s^{-1} \cdot cm^{-2} \cdot \AA^{-1}`, :math:`S_X`
is the complete telescope/camera/filter/atmosphere transmission function, and
:math:`ZP_X` is the zero-point for the filter :math:`X`. The term :math:`ch`
is just the Planck constant times the speed of light. To make the units work
out, :math:`ch = 1.9864\times 10^{-8} \ erg \cdot \AA`. The zero-point
:math:`ZP_X` is therefore the magnitude of an object that produces a 
photon flux of :math:`1\ s^{-1}\cdot cm^{-1}` through the instrument.

The first thing you need to know is :math:`S_X\left(\lambda\right)`. This may
be directly measured or approximated by multiplying the efficiencies of
the various optical components of your instrument. The atmosphere is usually
the trickiest to model. However you do it, you should tabulate the 
transmission as a function of wavelength (in Angstroms) and have these in an 
ASCII file. The normalization of the transmission is unimportant as it can
simply be absorved into the zero-point.

The second step is to determine the zero-point. Be careful here when using
published zero-points. They may be defined in a manner different than how
SNooPy defines them. Always refer to the equation above to check for
consistency. The easiest way to compute the zero-point yourself is to use a
spectrophotometric standard object's :math:`f_\lambda` and its magnitude. For
example, if you were using the [Landolt2007]_ standards for :math:`B` and
:math:`V`, then you could use the 
`CALSPEC <http://www.stsci.edu/hst/observatory/crds/calspec.html>`_ 
SED for Vega. But what to use for its magnitude? If you were on the
[Landolt2007]_ system, you can trace the calibration all the way back to
[JohnsonMorgan1953]_, who define :math:`B = V = 0.03`. If using
[Smith+2002]_ standards to calibrate your ugriz data, you could use the
CALSPEC SED of :math:`BD+17^\circ 4708` and its magnitudes from that paper.

Using AB Magnitudes with SNooPy
-------------------------------

The AB system doesn't refer to any specific spectrophotometric standard, but 
rather ties everything to a hypothetical source with constant 
:math:`f_\nu = 3631 Jy`. Mathematically, a broad-band AB magnitude is defined
as:

.. math::
   :nowrap:

   \begin{equation*}
      m_{AB} = -2.5\log_{10}\left[\frac{
    h^{-1}\int f_\nu\left(\nu\right) S_X\left(\nu\right) d\left(\ln \nu\right)}{
    h^{-1} 3631 Jy \int S_X\left(\nu\right) d\left(\ln \nu\right)}\right]
   \end{equation*}

Using this definition in the first equation, you can compute a zero-point that
will let SNooPy treat the AB magnitudes just like any other:

.. math::
   :nowrap:

   \begin{equation*}
      ZP_{AB} = 16.847 + 2.5\log_{10}\left[\int \frac{S_X\left(\lambda\right)}
                                           {\lambda} d\lambda\right]
   \end{equation*}

You can either perform this calculation yourself and use it as the zero-point
for your filter in SNooPy, or simply tell SNooPy that your magnitudes are AB
magnitudes and it will do this for you behind the scenes (see next section).

Putting it Into SNooPy
----------------------

At this point, you hopefully have:

   * one ASCII file for each filter, tabulating the wavelength (in Angstroms)
      and transmission.
   * one of the following:
      * computed zero-points for each filter, or:
      * a chosen spectrophotometric standard and its
         magnitude in each filter, or:
      * are using an AB system and are happy to let SNooPy figure out the
         zero-points for you.

The first thing you need to do is locate SNooPy's filter database. If you
downloaded the source, then the filter database is located in::
   
   /path/to/source/snpy/filters/filters

If you used the snpy-bootstrap.py script, the source folder will be in the
root of the virtual environment (either the one you set up or the existing
one you installed into)::

   /path/to/venv/snpy/filters/filters

In the ``filters`` folder, there is one folder for each observatory (e.g., ``LCO``, ``HST``, etc). Either pick one of these, or make a new one. Within an observatory's folder, there is a separate forder for each instrument. Make a new one
for your system. Copy your filters functions into this folder and create a new
file called ``filters.dat``.

The ``filters.dat`` file should have one line per filter (comments
are denoted by starting the line with a ``#``). Each line must have at least
four fields separated by spaces:  ``ID``, ``filename``, ``zero-point``, and
``comment``. Any extra fields will be considered to be part of the ``comment``.
The ``ID`` field should be unique (SNooPy will alert you if it is not). The
``filename`` should simply be the name of the filter's data file. The
``comment`` field should describe the filter. The ``zero-point`` can be 
one of the following:

   1) a floating-point number representing the zero-point
   2) a string of the form ``SED=mag``, where ``SED`` is the ID of
      a spectro-photometric standard in SNooPy (run ``standards.list_SEDs()``
      in SNooPy to see a list) and ``mag`` is its magnitude in your system.
   3) the string ``"AB"``, which indicates the filter is on the AB system.

Once you have all that set up, re-install SNooPy the usual way (``python
setup.py install``) or run the ``update-snpy`` script. Now, try starting
SNooPy. If there is a problem with the way you setup the filters, you'll get an
error. If not, you can now start using your filter ID when fitting data and
SNooPy will take care of the S-corrections.

.. _natural_v_standard:

Natural vs. Standard Photometry
-------------------------------

When specifying the magnitudes of your standard through your filters, you need
to know if you are on a standard system or a natural system. What's the
difference?

Standard photometry usually refers to a published photometric system defined
by a standard set of filters and a list of standard stars with magnitudes
through those filters distributed across the sky. Two very common systems
are the [Landolt2007]_ and [Smith+2002]_ systems. [Landolt2007]_ uses the
standard UBVRI filters and is tied to a set of A0 stars, which collectively
define zero colors and has Vega setting the absolute flux. [Smith+2002]_ 
uses the Sloan filters ugriz and is nominally an AB system, which uses a
hypotheical constant :math:`f_\nu = 3631 Jy` source, but is ultimately
calibrated against three standard stars :math:`BD+17^\circ 4708`,
:math:`BD+26^\circ 2606` and :math:`BD+21^\circ 0607` with magnitudes
tabulated in [Smith+2002]_.

An observer on a different telescope will have different filter functions than
those used by either of these standard systems. As such, there will be
systematic differences in the flux levels measured from these standards that is
a function of the color of the star. We can correct for these systematics by
measuring "color terms" empirically and applying them to the standards we
measure, thereby putting our photometry on the standard system.

There are two major problems with this approach:  1) the "standard" filters are
not well defined; and 2) the color terms, being determined using standard stars,
will only work for stars, not supernovae. Supernova observers have therefore
come up with the concept of "natural photometry".

Simply put, the natural photometric system is a system based on a set of 
standard stars whose magnitudes are what would have been measured through
the instrument+telescope being used. We can do this because we can use the 
color terms in reverse and transform the [Landolt2007]_ and [Smith+2002]_ 
standard magnitudes into natural magnitudes.

Put another way:  instead of using color terms to transform our instrumental
magnitudes to match the standard system, use the color terms in reverse to 
transform the standard magnitudes to match our instrumental system.

The reason you need to know this is when assigning magnitudes to the
fundamental standards. If you are on the standard system, then simply use
the magnitudes of the fundamental standards as published. If, on the other
hand, you are on a natural system, you need to transform the published
magnitudes of the fundamental standard to your natural system using the
color terms in reverse.

For example, consider the CSP :math:`B` and :math:`V` filters. CSP
photometry is in the natural system. According to [JohnsonMorgan1953]_,
:math:`B_{Vega} = V_{Vega} = 0.03`. The :math:`B` and :math:`V`
color terms for the CSP are:

.. math::
   :nowrap:

   \begin{eqnarray*}
      B_{CSP} &=& B - 0.061 \left(B - V\right) \\
      V_{CSP} &=& V + 0.058 \left(V - i^\prime\right) \\
   \end{eqnarray*}

And so :math:`B` for Vega in the CSP natural system would be the same as in the [Landolt2007]_ system: 0.03 mag. But :math:`V` requires an :math:`i^\prime` 
magnitude for Vega. For this we have to rely on synthetic [Smith+2002]_ 
photometry:  :math:`i^\prime_{Vega} = 0.382`. Using this in the color term
equation gives :math:`V_{Vega} = 0.0096` in the CSP natural system.

Photometry in one natural system is much easier to transform to another
natural system as we only need the filter functions and zero-points.


.. [Stritzinger+2005] Stritzinger et al., PASP, 117, 810 (2005)
   http://adsabs.harvard.edu/abs/2005PASP..117..810S
.. [Landolt2007] Landolt, A., AJ, 104, 340L (1992)
   http://adsabs.harvard.edu/abs/1992AJ....104..340L
.. [Smith+2002] Smith, A. et al., AJ, 123, 2121 (2002)
   http://adsabs.harvard.edu/abs/2002AJ....123.2121S
.. [JohnsonMorgan1953] Johnson, H.L. and Moran W.W., ApJ, 117, 313 (1953)

