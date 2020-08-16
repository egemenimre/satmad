# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Celestial constants to be used in SatMAD.

"""

from astropy import units as u

GM_earth = 3.986004418e14 * u.m ** 3 / u.s ** 2
"""Geocentric gravitational constant
(TCG-compatible value).

(See IERS Technical Note No. 36 (2010) Table 1.1 [TCF1]_)"""

GM_sun = 1.32712442099e20 * u.m ** 3 / u.s ** 2
"""Heliocentric gravitational constant
(TCB-compatible value, computed from the TDB-compatible value).

(See IERS Technical Note No. 36 (2010) Table 1.1  [TCF1]_)"""

GM_moon = 4.902802711497899e12 * u.m ** 3 / u.s ** 2
r"""Lunar  gravitational constant.

The value is derived from "Moon-Earth mass ratio" times "Geocentric
gravitational constant" (:math:`0.0123000371 \times 398600.64418`).

(See IERS Technical Note No. 36 (2010) Table 1.1  [TCF1]_)."""


class Ellipsoid:
    """Defines an ellipsoid, the shape of most celestial bodies.

    Parameters
    ----------
    re : Quantity
        Equatorial Radius or Semimajor Axis of the Ellipsoid. (km)
    inv_f : Quantity
        Inverse flattening (:math:`1/f`). (dimensionless)
    """

    @u.quantity_input(re=u.km, inv_f=u.dimensionless_unscaled)
    def __init__(self, re, inv_f):
        self.re = re
        self.inv_f = inv_f


EARTH_ELLIPSOID_IERS2003 = Ellipsoid(
    6378136.6 * u.m, 298.25642 * u.dimensionless_unscaled
)
"""Earth Ellipsoid defined in IERS 2010 Numerical Standards.

(See IERS Technical Note No. 36 (2010) Table 1.1  [TCF1]_)."""

EARTH_ELLIPSOID_GRS80 = Ellipsoid(
    6378137 * u.m, 298.257222101 * u.dimensionless_unscaled
)
"""Earth Ellipsoid defined in Geodetic Reference System GRS80.

(See IERS Technical Note No. 36 (2010) Table 1.2  [TCF1]_)."""

EARTH_ELLIPSOID_WGS84 = Ellipsoid(
    6378137 * u.m, 298.257223563 * u.dimensionless_unscaled
)
"""Earth Ellipsoid defined in World Geodetic System WGS84.

(See `WGS84 definition <https://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350_2.html>`_)."""
