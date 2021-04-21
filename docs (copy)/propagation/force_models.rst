Force Models
======================

Overview
------------
The force models in SatMAD are used generally within the numerical propagation algorithms
(see :ref:`numprop-intro`).

Two-Body Force Model
---------------------

The two-body acceleration is given by:

    .. math::

        \ddot{\vec{r}} = - \dfrac{\mu}{r^3} \vec{r}

where :math:`\ddot{\vec{r}}` is the inertial acceleration, :math:`\vec{r}` is the position vector
(with :math:`r` its norm). :math:`\mu` is equal to the constant :math:`GM` (Gravitational Constant times
the Mass of the main attracting body).

The equation holds regardless of whether the orbit is elliptic, parabolic or hyperbolic.

For an object orbiting the Earth, this equation should be solved with the coordinates in the GCRS frame.
For an object orbiting the Sun, the coordinates should be in the Heliocentric inertial frame, with the
:math:`\mu` selected as that of the Sun.

Reference/API
-------------
.. automodule:: satmad.propagation.force_models
    :members:
    :undoc-members: