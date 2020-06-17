# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
TLE class tests.

"""
import numpy as np
import pytest
from astropy import units as u
from astropy.time import Time
from pytest import approx

from satmad.propagation.tle import TLE


class TLEData:
    """TLE string data with name, line1 and line2"""

    def __init__(self, line1, line2, name="No Name"):
        self.line1 = line1
        self.line2 = line2
        self.name = name


@pytest.fixture
def tle_sso():
    """Test Fixture with TLE for SSO repeat groundtrack."""
    name = "SENTINEL 2A"
    line1 = "1 40697U 15028A   20164.50828565  .00000010  00000-0  20594-4 0  9999"
    line2 = "2 40697  98.5692 238.8182 0001206  86.9662 273.1664 14.30818200259759"
    return TLEData(line1, line2, name)


@pytest.fixture
def init_tle_sso(tle_sso):
    """Generates the TLE with SSO repeat groundtrack test setup."""
    return TLE.from_tle(tle_sso.line1, tle_sso.line2, tle_sso.name)


@pytest.fixture
def tle_leo():
    """Test Fixture with TLE for LEO."""
    name = "RASAT"
    line1 = "1 37791U 11044D   20160.49218133  .00000015  00000-0  11838-4 0  9998"
    line2 = "2 37791  98.1102 253.4779 0020835 300.2201  59.6937 14.64769942470937"
    return TLEData(line1, line2, name)


@pytest.fixture
def init_tle_leo(tle_leo):
    """Generates the TLE with LEO test setup."""
    return TLE.from_tle(tle_leo.line1, tle_leo.line2, tle_leo.name)


@pytest.fixture
def check_str_leo(tle_leo):
    """Generates the TLE string with LEO test setup."""
    return tle_leo.name + "\n" + tle_leo.line1 + "\n" + tle_leo.line2 + "\n"


@pytest.fixture
def tle_geo():
    """Test Fixture with GEO TLE."""
    name = "SAMPLE GEO"
    line1 = "1 99999U 12345A   20162.50918981  .00000000  00000-0  00000-0 0 00005"
    line2 = "2 99999 000.0000 124.6202 0000000 000.0000 000.0000 01.00273791000004"
    return TLEData(line1, line2, name)


@pytest.fixture
def init_tle_geo(tle_geo):
    """Generates the GEO TLE test setup."""
    return TLE.from_tle(tle_geo.line1, tle_geo.line2, tle_geo.name)


def test_tle_init_with_lines(tle_leo, check_str_leo):
    """Test initialise from TLE."""
    tle1 = TLE.from_tle(tle_leo.line1, tle_leo.line2, tle_leo.name)

    assert str(tle1) == check_str_leo


def test_tle_with_params(init_tle_leo, check_str_leo):
    """Test initialise from parameters."""

    tle1 = init_tle_leo

    tle2 = TLE(
        tle1.epoch,
        tle1.inclination,
        tle1.raan,
        tle1.eccentricity,
        tle1.arg_perigee,
        tle1.mean_anomaly,
        tle1.mean_motion,
        tle1.bstar,
        tle1.n_dot,
        n_dotdot=tle1.n_dotdot,
        name=init_tle_leo.name,
        intl_designator=tle1.intl_designator,
        sat_num=tle1.sat_number,
        classification=tle1.classification,
        rev_nr=tle1.rev_nr,
        el_nr=tle1.el_nr,
    )

    assert str(tle2) == check_str_leo == str(tle1)


def test_init_geo(init_tle_geo):
    """Test init GEO satellite."""
    tle_geo = init_tle_geo

    epoch = Time("2020-06-10T12:13:14.000")
    longitude = 42.0 * u.deg

    tle = TLE.init_geo(
        epoch,
        longitude,
        name=tle_geo.name,
        sat_num=tle_geo.sat_number,
        intl_designator=tle_geo.intl_designator,
        rev_nr=tle_geo.rev_nr,
        el_nr=tle_geo.el_nr,
    )

    assert str(tle_geo) == str(tle)


def test_tle_init_incl_out_of_bounds(init_tle_leo):
    """Tests init with inclination input value out of bounds
    - should raise `ValueError`."""
    with pytest.raises(ValueError):
        tle1 = init_tle_leo

        TLE(
            tle1.epoch,
            2 * np.pi,
            tle1.raan,
            tle1.eccentricity,
            tle1.arg_perigee,
            tle1.mean_anomaly,
            tle1.mean_motion,
            tle1.bstar,
            tle1.n_dot,
            n_dotdot=tle1.n_dotdot,
            name=init_tle_leo.name,
            intl_designator=tle1.intl_designator,
            sat_num=tle1.sat_number,
            classification=tle1.classification,
            rev_nr=tle1.rev_nr,
            el_nr=tle1.el_nr,
        )


def test_node_rot(init_tle_sso):
    """Test orbit plane rotation rate."""
    tle = init_tle_sso

    assert tle.node_rotation_rate().to_value(u.deg / u.day) == approx(
        (0.9870658041317965 * u.deg / u.day).to_value(),
        abs=(1e-14 * u.deg / u.day).to_value(),
    )


def test_argp_rot(init_tle_sso):
    """Test argument of perigee rotation rate."""
    tle = init_tle_sso

    assert tle.argp_rotation_rate().to_value(u.deg / u.day) == approx(
        (-2.9445253809901057 * u.deg / u.day).to_value(),
        abs=(1e-14 * u.deg / u.day).to_value(),
    )


def test_incl_out_of_bounds(init_tle_leo):
    """Tests inclination setter with input value out of bounds
    - should raise `ValueError`."""
    with pytest.raises(ValueError):
        tle = init_tle_leo
        tle.inclination = 210 * u.deg


def test_raan_out_of_bounds(init_tle_leo):
    """Tests RAAN setter with input value out of bounds """
    tle = init_tle_leo
    tle.raan = 390 * u.deg
    assert tle.raan.to_value(u.deg) == approx((30 * u.deg).to_value(), abs=1e-8)


def test_getters_setters(init_tle_leo):
    """Test getters and setters of the TLE."""
    tle = init_tle_leo

    # test getters
    assert tle.sm_axis().to_value() == approx(
        (7055.953203777368 * u.km).to_value(), abs=(0.01 * u.mm).to_value(u.km)
    )
    assert tle.period().to_value(u.s) == approx(
        (5898.53720524 * u.s).to_value(), abs=(1 * u.us).to_value(u.s)
    )
    assert tle.inclination.to_value(u.deg) == approx(
        (98.1102 * u.deg).to_value(), abs=1e-6
    )
    assert tle.raan.to_value(u.deg) == approx((253.4779 * u.deg).to_value(), abs=1e-8)
    assert tle.eccentricity == approx(0.0020835, abs=1e-9)
    assert tle.mean_anomaly.to_value(u.deg) == approx(
        (59.6937 * u.deg).to_value(), abs=1e-8
    )
    assert tle.arg_perigee.to_value(u.deg) == approx(
        (300.2201 * u.deg).to_value(), abs=1e-8
    )
    assert tle.n_dot == approx(4.5451282604019e-13, rel=1e-10)

    # test setters
    tle.inclination = 30.00 * u.deg
    assert tle.inclination.to_value(u.deg) == approx(
        (30.00 * u.deg).to_value(), abs=1e-8
    )

    tle.raan = 230.00 * u.deg
    assert tle.raan.to_value(u.deg) == approx((230.00 * u.deg).to_value(), abs=1e-8)

    tle.arg_perigee = 20.00 * u.deg
    assert tle.arg_perigee.to_value(u.deg) == approx(
        (20.00 * u.deg).to_value(), abs=1e-8
    )

    tle.mean_anomaly = 20.00 * u.deg
    assert tle.mean_anomaly.to_value(u.deg) == approx(
        (20.00 * u.deg).to_value(), abs=1e-8
    )

    tle.eccentricity = 0.05
    assert tle.eccentricity == approx(0.05, abs=1e-10)

    tle.mean_motion = 1.05e-5
    assert tle.mean_motion == approx(1.05e-5, rel=1e-10)

    tle.n_dot = 1.05e-5
    assert tle.n_dot == approx(1.05e-5, rel=1e-10)

    tle.el_nr = 123
    assert tle.el_nr == 123

    tle.rev_nr = 123
    assert tle.rev_nr == 123

    tle.classification = "S"
    assert tle.classification == "S"

    tle.name = "New Name"
    assert tle.name == "New Name"
