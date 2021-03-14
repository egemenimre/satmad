Coordinate Systems and Frames
==================================

In SatMAD, frames and coordinate systems as well as conversions between them are handled through
`Astropy <https://docs.astropy.org/en/latest/coordinates/index.html>`_. This section introduces the
additional frames defined by SatMAD.

Inertial Frames of Celestial Bodies
-----------------------------------
In addition to Earth, it is possible to define a Celestial Body in space through the
:class:`.CelestialBody` class. For each such Celestial Body, it is then possible to realise an inertial
or "Celestial Reference System" (CRS). This enables the user to define coordinates in this local
inertial coordinate system and then run an orbit propagation around it. For example, the following
would define a "Sun CRS" (equivalent to Heliocentric Celestial Reference System) and a "Moon CRS"
by simply subclassing :class:`.CelestialBodyCRS`.

.. code-block:: python
    :linenos:

    from satmad.coordinates.frames import CelestialBodyCRS
    from satmad.core.celestial_bodies import MOON, SUN

    class SunCRS(CelestialBodyCRS):
        body = SUN

    class MoonCRS(CelestialBodyCRS):
        body = MOON

While Sun and Moon are pre-defined, it is possible to define planets like Jupiter or moons like Io or Ceres
by simply defining them as an instance of the :class:`.CelestialBody` class.

Earth Based Additional Frames (J2000)
----------------------------------------------

The built-in frames offered by `Astropy <https://docs.astropy.org/en/latest/coordinates/index.html>`_
do not include some frames that are used in satellite applications. To bridge this gap, this package
offers Mean Pole and Equinox at J2000.0 Reference System (:class:`.J2000`).

:class:`.J2000` coordinate frame is similar to GCRS but rotated by a constant frame bias
[TCF1]_:

.. math:: \vec{r}_{J2000} = B \times \vec{r}_{GCRS}

This rotation is applicable only to the equinox based approach, and is only an approximation.
The difference between GCRS and J2000 is less than 1m for the Low Earth Orbit, therefore these two
can be used interchangeably with a small error.

The :class:`.J2000` class is similar to (and compatible with) the `Astropy Built-in Frames
<https://docs.astropy.org/en/latest/coordinates/index.html#built-in-frame-classes>`_.

.. code-block:: python
    :linenos:

    from astropy import units as u
    from astropy.coordinates import CartesianRepresentation, CartesianDifferential
    from satmad.coordinates.j2000 import J2000

    v_j2000 = CartesianDifferential([-4.7432196000, 0.7905366000, 5.5337561900], unit=u.km / u.s)
    r_j2000 = CartesianRepresentation([5102.50960000, 6123.01152000, 6378.13630000], unit=u.km)
    rv_j2000 = J2000(r_j2000.with_differentials(v_j2000), obstime=time, representation_type="cartesian", differential_type="cartesian")


Reference / API
---------------

.. automodule:: satmad.coordinates.frames
    :members:
    :undoc-members:






