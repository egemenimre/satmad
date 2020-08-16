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
