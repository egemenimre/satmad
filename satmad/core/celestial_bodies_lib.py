# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Defines a Celestial Body, which acts as the central repository of information for
the planets or similar Central Bodies.

"""
from astropy import units as u
from astropy.coordinates import GCRS, HCRS, ITRS
from numpy import inf

from satmad.coordinates.frames import (
    CelestialBodyCRS,
    CelestialBodyFixed,
    CelestialBodyJ2000Equatorial,
    CelestialBodyTODEquatorial,
)
from satmad.core.celestial_body import CelestialBody, CelestialBodyEllipsoid

# **************** Celestial Body Coordinate System Definitions ****************


class MarsCRS(CelestialBodyCRS):
    """Mars Celestial Reference System. This is simply the ICRS shifted to the
    centre of the Mars with the velocity adjusted with respect to the Mars.
    """

    body_name = "Mars"
    ephemeris_type = "jpl"


class MarsTODEquatorial(CelestialBodyTODEquatorial):
    """
    A coordinate frame representing the True-of-Date Equatorial System of Mars.
    """

    body_name = "Mars"
    cb_crs = MarsCRS


class MarsJ2000Equatorial(CelestialBodyJ2000Equatorial):
    """
    A coordinate frame representing the Equatorial System of Mars at J2000 Epoch.
    """

    body_name = "Mars"
    cb_crs = MarsCRS


class MarsBodyFixed(CelestialBodyFixed):
    """
    A coordinate frame representing the Body Fixed System of Mars.
    """

    body_name = "Mars"
    cb_crs = MarsCRS


class MoonCRS(CelestialBodyCRS):
    """Moon Celestial Reference System. This is simply the ICRS shifted to the
    centre of the Moon with the velocity adjusted with respect to the Moon.

    This uses the ephemeris type `jpl` to ensure that Moon velocity is computed
    (`builtin` cannot compute velocity). This is critical for coordinate
    transformations that involve velocity.
    """

    body_name = "Moon"
    ephemeris_type = "jpl"


class MoonTODEquatorial(CelestialBodyTODEquatorial):
    """
    A coordinate frame representing the True-of-Date Equatorial System of Moon.
    """

    body_name = "Moon"
    cb_crs = MoonCRS


class MoonBodyFixed(CelestialBodyFixed):
    """
    A coordinate frame representing the Body Fixed System of Moon.
    """

    body_name = "Moon"
    cb_crs = MoonCRS


# **************** GM values ****************

GM_earth = 3.986004418e14 * u.m ** 3 / u.s ** 2
"""Geocentric gravitational constant
(TCG-compatible value). Also compatible with the WGS84 coordinate system.
Includes mass of the atmosphere.

(See IERS Technical Note No. 36 (2010) Table 1.1
[TCF1] in:doc:`References <../references>`)"""

GM_sun = 1.32712442099e20 * u.m ** 3 / u.s ** 2
"""Heliocentric gravitational constant
(TCB-compatible value, computed from the TDB-compatible value).

(See IERS Technical Note No. 36 (2010) Table 1.1
[TCF1] in :doc:`References <../references>`)"""

GM_moon = 4.902802711497899e12 * u.m ** 3 / u.s ** 2
r"""Lunar  gravitational constant.

The value is derived from "Moon-Earth mass ratio" times "Geocentric
gravitational constant" (:math:`0.0123000371 \times 398600.64418`).

(See IERS Technical Note No. 36 (2010) Table 1.1
[TCF1] in :doc:`References <../references>`)."""

# **************** Ellipsoid Definitions ****************

EARTH_ELLIPSOID_IERS2003 = CelestialBodyEllipsoid(
    "Earth Ellipsoid IERS 2003",
    6378136.6 * u.m,
    298.25642 * u.dimensionless_unscaled,
    # mu=3.986004418 * u.m ** 3 / u.s ** 2,
    j2=1.0826359e-3 * u.dimensionless_unscaled,
)
"""Earth Ellipsoid defined in IERS 2010 Numerical Standards.

(See IERS Technical Note No. 36 (2010) Table 1.1
[TCF1] in :doc:`References <../references>`)."""

EARTH_ELLIPSOID_GRS80 = CelestialBodyEllipsoid(
    "Earth Ellipsoid GRS80",
    6378137 * u.m,
    298.257222101 * u.dimensionless_unscaled,
    # mu=3.986005e14 * u.m ** 3 / u.s ** 2,
    j2=1.08263e-3 * u.dimensionless_unscaled,
    om=7.292115e-5 * u.rad / u.s,
)
"""Earth Ellipsoid defined in Geodetic Reference System GRS80.

(See IERS Technical Note No. 36 (2010) Table 1.2
[TCF1] in :doc:`References <../references>`)."""

EARTH_ELLIPSOID_WGS84 = CelestialBodyEllipsoid(
    "Earth Ellipsoid WGS84",
    6378137.0 * u.m,
    298.257223563 * u.dimensionless_unscaled,
    # mu=3.986004418e14 * u.m ** 3 / u.s ** 2,
    j2=1.082629821313e-3 * u.dimensionless_unscaled,
    om=7.292115e-5 * u.rad / u.s,
)
"""Earth Ellipsoid defined in World Geodetic System WGS84.

(See `WGS84 definition
<https://earth-info.nga.mil/php/download.php?file=coord-wgs84>`_)."""

MOON_ELLIPSOID_IAUWG2015 = CelestialBodyEllipsoid(
    "Moon Ellipsoid IAU WG 2015", 1737400 * u.m, inf * u.dimensionless_unscaled
)
"""Moon Ellipsoid - flattening is zero as the reference value means a perfect sphere.

(See Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements
(2015) Table 5  [TCF2] in :doc:`References <../references>`)."""

# **************** Celestial Body Definitions ****************

EARTH = CelestialBody(
    "Earth",
    "Default Earth Model.",
    GM_earth.to(u.km ** 3 / u.s ** 2),
    inert_coord=GCRS,
    body_fixed_coord=ITRS,
    ellipsoid=EARTH_ELLIPSOID_WGS84,
)
"""Default Earth Model."""

SUN = CelestialBody(
    "Sun",
    "Default Sun Model.",
    GM_sun.to(u.km ** 3 / u.s ** 2),
    inert_coord=HCRS,
    ellipsoid=CelestialBodyEllipsoid(
        "Sun Ellipsoid IAU Resolution B3 2015",
        695700 * u.km,
        inf * u.dimensionless_unscaled,
    ),
)
"""Default Sun Model.

Sun Ellipsoid is based on
"Mamajek, E.E.; Prsa, A.; Torres, G.; et, al. (2015), "IAU 2015 Resolution B3
on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties"
<https://arxiv.org/pdf/1510.07674.pdf>
"""


class MoonJ2000Equatorial(CelestialBodyJ2000Equatorial):
    """
    A coordinate frame representing the Equatorial System of Moon at J2000 Epoch.

    This corresponds to the Moon Principal Axis (PA) system. This has a fixed and small
    rotational offset with respect to the Mean Earth/rotation (ME) system.
    """

    body_name = "Moon"
    cb_crs = MoonCRS


MOON = CelestialBody(
    "Moon",
    "Default Moon Model.",
    GM_moon.to(u.km ** 3 / u.s ** 2),
    ellipsoid=MOON_ELLIPSOID_IAUWG2015,
    inert_coord=MoonJ2000Equatorial,
)
"""Default Moon Model."""
