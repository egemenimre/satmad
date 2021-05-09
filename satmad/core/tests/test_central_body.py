# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Tests the Central Body and Celestial Body modules.

"""
from astropy import units as u
from astropy.coordinates import (
    get_body_barycentric_posvel,
    CartesianDifferential,
    get_body_barycentric,
)
from astropy.time import Time

from satmad.core.celestial_bodies_lib import SUN, MOON


def test_coord_list():
    """Tests the get_coord_list method in both pos_vel and pos versions."""

    time_list = Time("2020-01-01T11:30:00", scale="utc") + range(0, 500, 27) * u.hour

    pv_list = SUN.get_coord_list(time_list, velocity=True, ephemeris="builtin")

    r, v = get_body_barycentric_posvel(SUN.name, time_list, ephemeris="builtin")
    v_sun = CartesianDifferential(v.xyz)
    r_sun = r.with_differentials(v_sun)

    # noinspection PyTypeChecker
    assert all(pv_list.cartesian == r_sun)

    pv_list = MOON.get_coord_list(time_list, velocity=False, ephemeris="builtin")

    r_moon = get_body_barycentric(MOON.name, time_list, ephemeris="builtin")

    # noinspection PyTypeChecker
    assert all(pv_list.cartesian == r_moon)
