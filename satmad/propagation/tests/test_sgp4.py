# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
SGP4 propagator tests.

"""

import pytest
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord,
)

from satmad.coordinates.frames import J2000
from satmad.propagation.sgp4_propagator import (
    SGP4GeneralError,
    SGP4Propagator,
    SGP4SatDecayedError,
)
from satmad.propagation.tle import TLE


@pytest.fixture
def init_tle_decay():
    """Generates the TLE with decay test setup."""

    name = "MICROSAT-R DEB"
    line1 = "1 44160U 19006AX  20178.55672017  .00609591  12600-3  17697-2 0  9996"
    line2 = "2 44160  95.2356 283.8798 0142512 327.7074  31.5529 15.78083811 64981"

    return TLE.from_tle(line1, line2, name)


@pytest.fixture
def init_tle_vgd():
    """Generates the TLE with LEO test setup."""

    name = "VANGUARD 1"
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    return TLE.from_tle(line1, line2, name)


def test_propagation(init_tle_vgd):
    # init TLE
    tle = init_tle_vgd

    # init SGP4
    sgp4 = SGP4Propagator()

    #  Propagate 3 days into future
    time = tle.epoch + 3.0 * u.day
    rv_gcrs = sgp4.propagate(tle, time)

    # print(rv_gcrs)

    # Generate truth values
    # Vallado IAU-76/FK5 - pg.234
    v_j2000_true = CartesianDifferential(
        [-2.233348094, -4.110136162, -3.157394074], unit=u.km / u.s
    )
    r_j2000_true = CartesianRepresentation(
        [-9059.9413786, 4659.6972000, 813.9588875], unit=u.km
    )

    rv_j2000_true = J2000(
        r_j2000_true.with_differentials(v_j2000_true),
        obstime=time,
        representation_type="cartesian",
        differential_type="cartesian",
    )

    rv_gcrs_true = SkyCoord(
        rv_j2000_true.transform_to(GCRS(obstime=time)),
        representation_type="cartesian",
        differential_type="cartesian",
    )

    r_diff = (
        rv_gcrs.cartesian.without_differentials()
        - rv_gcrs_true.cartesian.without_differentials()
    )

    v_diff = rv_gcrs.velocity - rv_gcrs_true.velocity

    # print(r_diff.norm().to(u.mm))
    # print(v_diff.norm().to(u.mm / u.s))

    assert r_diff.norm().to(u.mm) < 1700 * u.mm
    assert v_diff.norm().to(u.mm / u.s) < 0.71 * u.mm / u.s


def test_decay(init_tle_decay):
    """Tests satellite altitude below threshold - satellite has decayed.
        - should raise `SGP4SatDecayedError`."""
    with pytest.raises(SGP4SatDecayedError):
        # init TLE
        tle = init_tle_decay

        #  Propagate 3 days into future
        time = tle.epoch + 30.0 * u.day
        SGP4Propagator.propagate(tle, time)


def test_crashed(init_tle_decay):
    """Tests satellite altitude below earth radius - satellite has crashed.
        - should raise `SGP4GeneralError`."""
    with pytest.raises(SGP4GeneralError):
        # init TLE
        tle = init_tle_decay

        #  Propagate 3 days into future
        time = tle.epoch + 300.0 * u.day
        SGP4Propagator.propagate(tle, time)
