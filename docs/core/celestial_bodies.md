# "Celestial Bodies" in Space

## Overview

A Celestial Body is an object in space around which another object can rotate. As such, it should have certain properties to enable operations such as coordinate transformations and propagation.

Celestial bodies in SatMAD are defined with the {py:class}`.CelestialBody` class. A celestial body stores its relevant values (e.g Gravitational Constant $GM$) for various calculations such as force models within orbit propagation. As such, each `CelestialBody` object is a central repository of information for that body.

The `celestial_bodies` module defines the default Celestial Bodies, which are instances of the {py:class}`.CelestialBody` class. These  default celestial bodies are currently `EARTH`, `MOON` and `SUN`. The properties of these Celestial bodies are defined using certain standards and established values (e.g. IERS Technical Note No. 36 [[TCF1]](../references.md#time-and-coordinate-frames)). However, they are different from the standard values offered by Astropy. A simple usage example is:

    >>> from satmad.core.celestial_bodies import SUN
    >>> SUN.mu
    <Quantity 1.32712442e+11 km3 / s2>
    >>> from satmad.core.celestial_bodies import EARTH
    >>> EARTH.ellipsoid.re
    <Quantity 6378137. m>

A `CelestialBody` object should include, as a minimum, `name`, `info` and `mu` (Gravitational Constant). In addition, it can have an `ellipsoid` (e.g. GRS80 for the Earth) as well as default coordinate definitions `inert_coord` and `body_fixed_coord`. The inertial coordinate definition (`inert_coord`) is required to be able to run a propagation around this central body. For the Earth, this is set to GCRS. The body fixed coordinate (`body_fixed_coord`) is required for conversions between the inertial and body fixed coordinate frames (for example, to be able to compute where the satellite sensor is pointing at on the ground). For the Earth, this is set to ITRS.

Apart from the properties that can be queried, the `CelestialBody` object offers the {py:meth}`.get_coord_list` method. This method provides positions (and velocities) of planets and other celestial bodies. For example, the following code snippet computes the position of the Sun in ICRS frame every hour for a day:

    >>> from astropy.time import Time
    >>> from astropy import units as u
    >>> time_list = Time("2020-01-01T11:30:00", scale="utc") + range(0, 24, 1) * u.hour
    >>> pv_list = SUN.get_coord_list(time_list, velocity=True, ephemeris="builtin")
    >>> pv_list_gcrs = pv_list.transform_to("gcrs")

This is just a convenience method to wrap the {py:meth}`~astropy.coordinates.get_body_barycentric_posvel` method and return a SkyCoord object. Then the results can be converted to another coordinate system (in this example, to GCRS). Velocity output, or the ephemeris type (e.g. `builtin` or `jpl`) can also be specified.

## Reference/API

```{eval-rst}
.. automodule:: satmad.core.celestial_bodies_lib
    :undoc-members:
    :members:
```

```{eval-rst}
.. automodule:: satmad.core.celestial_body
    :members:
    :undoc-members:
```
