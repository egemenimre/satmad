# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Common test functions.

"""
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianDifferential,
    CartesianRepresentation,
)
from astropy.time import Time
from numpy import inf

from satmad.coordinates.frames import MoonCRS, init_pvt
from satmad.core.celestial_body import CelestialBody, CelestialBodyEllipsoid

GMAT_MOON = CelestialBody(
    "Moon",
    "GMAT Moon Model.",
    4902.8005821478 * (u.km ** 3 / u.s ** 2),
    ellipsoid=CelestialBodyEllipsoid(
        "GMAT Moon Ellipsoid", 1738.2 * u.km, inf * u.dimensionless_unscaled
    ),
    inert_coord=MoonCRS,
)

GMAT_EARTH = CelestialBody(
    "Earth",
    "GMAT Earth Model.",
    398600.4415 * (u.km ** 3 / u.s ** 2),
    inert_coord=GCRS,
    body_fixed_coord=ITRS,
    ellipsoid=CelestialBodyEllipsoid(
        "GMAT Earth Ellipsoid",
        6378.1363 * u.km,
        298.26706833298533 * u.dimensionless_unscaled,
    ),
)


def pos_err(rv_test, rv_true):
    """
    Computes positional error between two vectors defined in a frame.

    Parameters
    ----------
    rv_test : State with position vector under test
    rv_true : State with true position vector

    Returns
    -------
    Quantity
        Norm of the position difference

    """
    r_diff = (
        rv_test.cartesian.without_differentials()
        - rv_true.cartesian.without_differentials()
    )
    return r_diff.norm().to(u.mm)


def vel_err(rv_test, rv_true):
    """
    Computes velocity error between two vectors defined in a frame.

    Parameters
    ----------
    rv_test : State with velocity vector under test
    rv_true : State with true velocity vector

    Returns
    -------
    Quantity
        Norm of the velocity difference

    """
    v_diff = rv_test.velocity - rv_true.velocity
    return v_diff.norm().to(u.mm / u.s)


def parse_rv_line(rv_line, coord_sys):
    """Converts a line of t, r, v string into a SkyCoord object."""

    rv_items = [elem for elem in rv_line.strip().split(" ") if elem]

    time = Time(rv_items[0], scale="utc")

    v = CartesianDifferential(
        [float(rv_items[4]), float(rv_items[5]), float(rv_items[6])], unit=u.km / u.s,
    )
    r = CartesianRepresentation(
        [float(rv_items[1]), float(rv_items[2]), float(rv_items[3])], unit=u.km
    )

    rv_init = init_pvt(coord_sys, time, r.with_differentials(v))

    return rv_init
