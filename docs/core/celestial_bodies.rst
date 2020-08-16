Celestial Bodies
==================
Celestial bodies in SatMAD are defined with the :class:`.CelestialBody` class. A celestial body
stores its relevant values (e.g :math:`GM` constant) for various calculations such as force models
within orbit propagation. As such, each `CelestialBody` object is a central repository of
information for that body.

There are default celestial bodies such as `EARTH`, `MOON` and `SUN` for convenience.
A simple usage example is:

    >>> from satmad.core.celestial_bodies import EARTH
    >>> EARTH.mu
    <Quantity 398600.4418 km3 / s2>


Reference/API
-------------
.. automodule:: satmad.core.celestial_bodies
    :undoc-members:
    :members:
