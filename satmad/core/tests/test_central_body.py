# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Tests the Central Body and Celestial Body modules.

"""
from astropy.coordinates import get_body_barycentric_posvel
from astropy.time import Time
from astropy import units as u

from satmad.core.celestial_bodies_lib import SUN


def test_coord_list():
    time_list = Time("2020-01-01T11:30:00", scale="utc") + range(0, 500, 27) * u.hour

    pv_list = SUN.get_coord_list(time_list, velocity=True, ephemeris="builtin")

    pv_list_astropy = get_body_barycentric_posvel(
        SUN.name, time_list, ephemeris="builtin"
    )

    print(pv_list)
