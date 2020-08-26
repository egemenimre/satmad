Central Body
====================
.. _central_body-intro:

Celestial bodies in SatMAD are defined with the :class:`.CelestialBody` class. A celestial body
stores its relevant values (e.g :math:`GM` constant) for various calculations such as force models
within orbit propagation. As such, each `CelestialBody` object is a central repository of
information for that body.

There are default celestial bodies such as `EARTH`, `MOON` and `SUN` for convenience.
A simple usage example is:

    >>> from satmad.core.celestial_bodies import SUN
    >>> SUN.mu
    <Quantity 1.32712442e+11 km3 / s2>
    >>> from satmad.core.celestial_bodies import EARTH
    >>> EARTH.ellipsoid.re
    <Quantity 6378137. m>

A `CelestialBody` object should include, as a minimum, `name`, `info` and `mu` (Gravitational constant). In addition,
it can have an `ellipsoid` (e.g. GRS80 for the Earth) as well as default coordinate definitions `inert_coord` and
`body_fixed_coord`. The inertial coordinate definition (`inert_coord`) is required to be able to run a propagation
around this central body. For the Earth, this is set to GCRS. The body fixed coordinate (`body_fixed_coord`) is
required for conversions between the inertial and body fixed coordinate frames (for example, to be able to compute
where the satellite sensor is pointing at on the ground). For the Earth, this is set to ITRS.


Reference/API
-------------
.. automodule:: satmad.core.central_body
    :members:
    :undoc-members: