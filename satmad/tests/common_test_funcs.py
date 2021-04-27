# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Common test functions.

"""
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from astropy.time import Time

from satmad.coordinates.frames import init_pvt


def parse_rv_line(rv_line, coord_sys):
    """Converts a line of t, r, v string into a SkyCoord object."""

    rv_items = [elem for elem in rv_line.strip().split(" ") if elem]

    time = Time(rv_items[0], scale="utc")

    v = CartesianDifferential(
        [float(rv_items[4]), float(rv_items[5]), float(rv_items[6])],
        unit=u.km / u.s,
    )
    r = CartesianRepresentation(
        [float(rv_items[1]), float(rv_items[2]), float(rv_items[3])], unit=u.km
    )

    rv_init = init_pvt(coord_sys, time, r.with_differentials(v))

    return rv_init
