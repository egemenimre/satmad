# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Ground Location class tests.

"""
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import (
    ITRS,
    CartesianRepresentation,
    EarthLocation,
    Latitude,
    Longitude,
)
from astropy.time import Time, TimeDelta
from pytest import approx

from satmad.core.celestial_bodies_lib import (
    EARTH,
    EARTH_ELLIPSOID_GRS80,
    EARTH_ELLIPSOID_WGS84,
    MOON,
)
from satmad.core.ground_location import GeodeticLocation, GroundLocation


def test_get_body_fixed():
    """Tests the body fixed coords output for multiple times."""
    time = Time("2020-04-10T00:00:00", scale="utc") + np.arange(1, 4) * TimeDelta(
        1, format="jd"
    )

    gnd_loc = GroundLocation(
        10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=EARTH_ELLIPSOID_GRS80
    )

    gnd_loc_body_fixed = gnd_loc.to_body_fixed_coords(EARTH, obstime=time)

    r_itrs_true = ITRS(
        CartesianRepresentation(
            np.asarray(
                [
                    [6068714.27712043, 1070078.065267, 1640138.96300037],
                    [6068714.27712043, 1070078.065267, 1640138.96300037],
                    [6068714.27712043, 1070078.065267, 1640138.96300037],
                ]
            ).transpose(),
            unit=u.m,
        ),
        obstime=time,
        representation_type="cartesian",
        differential_type="cartesian",
    )

    # print(gnd_loc_body_fixed)

    assert str(gnd_loc_body_fixed) == str(r_itrs_true)


def test_get_body_fixed_no_coords():
    """Tests the body fixed coords output but the celestial body has no
    body fixed coordinates defined."""
    with pytest.raises(TypeError):

        gnd_loc = GroundLocation(
            10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=MOON.ellipsoid
        )

        gnd_loc.to_body_fixed_coords(MOON)


def test_copy():
    """Tests the initialisation with copy."""
    gnd_loc = GroundLocation(
        10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=EARTH_ELLIPSOID_WGS84
    )

    gnd_loc_copy = GroundLocation(gnd_loc)

    assert gnd_loc == gnd_loc_copy


def test_init_fail_1():
    """Tests the initialisation failure case with cartesian input error."""
    with pytest.raises(TypeError):
        gnd_loc = GroundLocation(
            10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=EARTH_ELLIPSOID_WGS84
        ).to_geocentric

        GroundLocation(gnd_loc)


def test_init_fail_2():
    """Tests the initialisation failure case with cartesian input with no units error."""
    with pytest.raises(TypeError):
        gnd_loc = GroundLocation(
            10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=EARTH_ELLIPSOID_WGS84
        ).to_geocentric()

        # strip unit information
        GroundLocation.from_geocentric(
            gnd_loc[0].to_value(),
            gnd_loc[1].to_value(),
            gnd_loc[2].to_value(),
        )


def test_round_trip():
    """Round trip testing from geodetic to geocentric to geodetic."""
    gd_loc_origin = GeodeticLocation(
        Longitude(10 * u.deg), Latitude(15 * u.deg), 150 * u.m
    )

    # the values are in cartesian coordinates
    gd_loc_copy = GroundLocation(gd_loc_origin).geodetic

    assert gd_loc_origin.lat.to_value(u.deg) == approx(
        gd_loc_copy.lat.to_value(u.deg), abs=1e-14
    )

    assert gd_loc_origin.lon.to_value(u.deg) == approx(
        gd_loc_copy.lon.to_value(u.deg), abs=1e-14
    )

    assert gd_loc_origin.height.to_value(u.m) == approx(
        gd_loc_copy.height.to_value(u.m), abs=1e-8
    )


def test_lat_lon_alt_init():
    """Tests the GroundLocation geodetic init against astropy EarthLocation."""
    gnd_loc_astropy = EarthLocation.from_geodetic(
        10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid="GRS80"
    )

    gnd_loc_gd = GroundLocation.from_geodetic(
        10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=EARTH_ELLIPSOID_GRS80
    )

    gnd_loc = GroundLocation(
        10 * u.deg, 15 * u.deg, 150 * u.m, ellipsoid=EARTH_ELLIPSOID_GRS80
    )

    print(gnd_loc_astropy.geocentric)
    # print(gnd_loc_gd)

    assert gnd_loc_gd == gnd_loc == gnd_loc_astropy


def test_cartesian_init():
    """Tests the GroundLocation cartesian init against astropy EarthLocation."""
    gnd_loc_astropy = EarthLocation.from_geodetic(10 * u.deg, 15 * u.deg, 1500 * u.m)

    gnd_loc_astropy_xyz = gnd_loc_astropy.to_geocentric()

    gnd_loc_gd = GroundLocation.from_geocentric(
        gnd_loc_astropy_xyz[0],
        gnd_loc_astropy_xyz[1],
        gnd_loc_astropy_xyz[2],
        ellipsoid=EARTH_ELLIPSOID_WGS84,
    )

    gnd_loc = GroundLocation(
        gnd_loc_astropy_xyz[0],
        gnd_loc_astropy_xyz[1],
        gnd_loc_astropy_xyz[2],
        ellipsoid=EARTH_ELLIPSOID_WGS84,
    )

    # print(gnd_loc_astropy)
    # print(gnd_loc_gd)

    assert gnd_loc_gd == gnd_loc == gnd_loc_astropy
