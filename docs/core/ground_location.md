# Ground Location

## Overview

A Ground Location is simply a location on a planet (or, more specifically, on the ellipsoid of a {py:class}`.CelestialBody`). The ground locations are usually defined with their geodetic positions (latitude, longitude and height on the ellipsoid) or with the cartesian coordinates fixed on the body fixed frame (e.g. ITRS, the Earth Fixed Coordinate System).

The {py:class}`.GroundLocation` class is a generalisation of the `EarthLocation` class in astropy and the code is mostly based on it. However, unlike the `EarthLocation` class which is bounded on the Earth, the `GroundLocation` can be anywhere, therefore an ellipsoid (e.g. Earth or Moon) should be explicitly specified. Otherwise, "Earth GRS80 Ellipsoid" is assumed as default.

A {py:class}`.GroundLocation` object can be initialised using a cartesian position vector (with length units) or with latitude, longitude and altitude / height values as seen in the examples below. The {py:class}`.GeodeticLocation` class ensures the latitude and longitude angles wrap correctly.

```python
from astropy import units as u
from astropy.coordinates import Latitude, Longitude

from satmad.core.celestial_bodies_lib import MOON_ELLIPSOID_IAUWG2015
from satmad.core.ground_location import GeodeticLocation, GroundLocation

# Ground location on the Moon, defined with lat/lon
gnd_loc_moon = GroundLocation(10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=MOON_ELLIPSOID_IAUWG2015)

# Ground location on the Earth with default ellipsoid, defined with GeodeticLocation
gnd_loc_earth = GroundLocation(GeodeticLocation(Longitude(10 * u.deg), Latitude(15 * u.deg), 150 * u.m))

# Ground location on the Earth with default ellipsoid, defined with cartesian position
gnd_loc_cart = GroundLocation.from_geocentric(5000 * u.km, 3600 * u.km, 1500 * u.km)
```

Once initialised, the {py:class}`.GroundLocation` object can yield its geodetic position with the `geodetic` parameter or `to_geodetic()` method. The result is a {py:class}`.GeodeticLocation` class.

Similarly, it can yield its geocentric cartesian position with the `geocentric` parameter or `to_geocentric()` method. This geocentric position is equal to the position in the Central Body Fixed frame, for example, ITRS for the Earth. The result is a three-element tuple of cartesian position. Perhaps more useful version of this is the `to_body_fixed_coords()` method, where the Geodetic Position is converted to a geocentric coordinates object as used by astropy (such as `ITRS` class for Earth). For this, one or more time values are to be provided, as well as a Celestial Body `body_fixed_coord` to generate the resulting object. If there are no body-fixed coordinates are defined for the Celestial Body, then a `TypeError` will be raised.

```python
import numpy as np

from astropy import units as u
from astropy.time import Time, TimeDelta

from satmad.core.celestial_bodies_lib import EARTH
from satmad.core.ground_location import GroundLocation

# Generate three time values, one day apart
time = Time("2020-04-10T00:00:00", scale="utc") + np.arange(1, 4) * TimeDelta(1, format="jd")

# Generate the Ground Location with the default Earth ellipsoid
gnd_loc = GroundLocation(10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=EARTH.ellipsoid, body_fixed_frame = EARTH.body_fixed_coord_frame)

# Generate the ground locations in body fixed frame (in this instance ITRS)
gnd_loc_body_fixed = gnd_loc.to_body_fixed_coords(obstime=time)

print(gnd_loc_body_fixed)
```

## Reference/API

```{eval-rst}
.. automodule:: satmad.core.ground_location
    :undoc-members:
    :members:
```