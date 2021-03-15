Ground Location
==================

A Ground Location is simply a location on a planet (or, more specifically, on the ellipsoid of a
:class:`.CelestialBody`). The ground locations are usually defined with their geodetic positions (latitude, longitude
and height on the ellipsoid) or with the cartesian coordinates fixed on the body fixed frame (e.g. ITRS, the Earth
Fixed Coordinate System).

The :class:`.GroundLocation` class is a generalisation of the `EarthLocation` class in astropy and the code is mostly
based on it. However, unlike the `EarthLocation` class which is bounded on the Earth, the `GroundLocation` can be
anywhere, therefore an ellipsoid (e.g. Earth or Moon) should be explicitly specified. Otherwise "Earth GRS80 Ellipsoid"
is assumed as default.

A :class:`.GroundLocation` object can be initialised using a cartesian position vector (with length units) or with
latitude, longitude and altitude / height values as seen in the examples below. The :class:`.GeodeticLocation` class
ensures the latitude and longitude angles wrap correctly.

.. code-block:: python
    :linenos:

    # Ground location on the Moon, defined with lat/lon
    gnd_loc_moon = GroundLocation(10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=MOON_ELLIPSOID_IAUWG2015)

    # Ground location on the Earth with default ellipsoid, defined with GeodeticLocation
    gnd_loc_earth = GroundLocation(GeodeticLocation(Longitude(10 * u.deg), Latitude(15 * u.deg), 150 * u.m))

    # Ground location on the Earth with default ellipsoid, defined with cartesian position
    gnd_loc_cart = GroundLocation.from_geocentric(5000 * u.km, 3600 * u.km, 1500 * u.km)


Once initialised, the :class:`.GroundLocation` object can yield its geodetic position with the `geodetic` parameter or
`to_geodetic()` method. The result is a :class:`.GeodeticLocation` class.

Similarly, it can yield its geocentric cartesian position with the `geocentric` parameter or
`to_geocentric()` method. This geocentric position is equal to the position in the Central Body Fixed frame,
for example, ITRS for the Earth. The result is a three-element tuple of cartesian position.

Reference/API
-------------
.. automodule:: satmad.core.ground_location
    :undoc-members:
    :members: