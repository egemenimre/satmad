# Celestial Bodies

## Overview 

A Celestial Body is an object in space around which another object can rotate. As such, it should have certain properties to enable operations such as coordinate transformations and propagation.

This module defines the default Celestial Bodies, which are instances of the {py:class}`.CelestialBody` class. These  default celestial bodies are currently `EARTH`, `MOON` and `SUN`. The properties of these Celestial bodies are defined using certain standards and established values (e.g. IERS Technical Note No. 36 [[TCF1]](../references.md#time-and-coordinate-frames)). However, they are different from the standard values offered by Astropy.

A simple usage example is:

    >>> from satmad.core.celestial_bodies import SUN
    >>> SUN.mu
    <Quantity 1.32712442e+11 km3 / s2>
    >>> from satmad.core.celestial_bodies import EARTH
    >>> EARTH.ellipsoid.re
    <Quantity 6378137. m>

## Reference/API

```{eval-rst}
.. automodule:: satmad.core.celestial_bodies
    :undoc-members:
    :members:
```