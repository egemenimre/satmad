# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Test orbit definitions.

"""

import pytest
from astropy import units as u
from astropy.coordinates import GCRS, ITRS
from astropy.time import Time
from pytest import approx

from satmad.coordinates.tests.test_baseframe_transform import pos_err, vel_err
from satmad.core.celestial_bodies_lib import EARTH
from satmad.core.celestial_body import CelestialBody
from satmad.propagation.classical_orb_elems import OsculatingKeplerianOrbElems
from satmad.tests.common_test_funcs import parse_rv_line

_EARTH_GMAT = CelestialBody(
    "Earth",
    "Earth as depicted in GMAT.",
    3.98600441500000e14 * u.m ** 3 / u.s ** 2,
    inert_coord=GCRS,
)


@pytest.fixture
def init_elems_leo():
    """Init a LEO example."""

    time = Time("2020-01-11T11:00:00.000", scale="utc")
    central_body = EARTH

    sm_axis = 7055.95320378 * u.km
    ecc = 0.0020835 * u.dimensionless_unscaled
    incl = 1.71234602 * u.rad
    raan = 4.42402394 * u.rad
    arg_perigee = 5.23982923 * u.rad
    true_an = 1.5 * u.rad

    orb_elems = OsculatingKeplerianOrbElems(
        time, sm_axis, ecc, incl, raan, arg_perigee, true_an, central_body
    )

    return orb_elems


def _parse_elems_line(epoch, elems_line, central_body):
    """Converts a line of elems string into a OsculatingOrbitalElems object."""

    elem_items = [elem for elem in elems_line.strip().split(" ") if elem]

    orb_elems = OsculatingKeplerianOrbElems(
        epoch,
        sm_axis=float(elem_items[1]) * u.km,
        eccentricity=float(elem_items[2]),
        inclination=float(elem_items[3]) * u.deg,
        raan=float(elem_items[4]) * u.deg,
        arg_periapsis=float(elem_items[5]) * u.deg,
        true_anomaly=float(elem_items[6]) * u.deg,
        central_body=central_body,
    )

    return orb_elems


def _elems_assert(osc_elems, truth_elems, errs):
    """Tests the orbital elements with given tolerances"""
    assert osc_elems.sm_axis.to_value(u.km) == approx(
        truth_elems.sm_axis.to_value(u.km), abs=errs[0].to_value(u.km)
    )

    assert osc_elems.eccentricity == approx(truth_elems.eccentricity, abs=errs[1])

    assert osc_elems.inclination.to_value(u.deg) == approx(
        truth_elems.inclination.to_value(u.deg), abs=errs[2].to_value(u.deg)
    )

    assert osc_elems.raan.to_value(u.deg) == approx(
        truth_elems.raan.to_value(u.deg), abs=errs[3].to_value(u.deg)
    )

    assert osc_elems.arg_periapsis.to_value(u.deg) == approx(
        truth_elems.arg_periapsis.to_value(u.deg), abs=errs[4].to_value(u.deg)
    )

    assert osc_elems.true_anomaly.to_value(u.deg) == approx(
        truth_elems.true_anomaly.to_value(u.deg), abs=errs[5].to_value(u.deg)
    )


def __init_ell_incl():
    """Case 1: Initialises an Earth Elliptical Inclined orbit satellite. GMAT example.

    This case deliberately uses an ITRS input, degrading the agreement between GMAT and SatMAD."""

    rv_line = "2010-01-01T13:24:28.100   3769.438525258299  6010.414642548296   1306.35068662678  -5.810744579719358   3.604099244744345   0.9947338414920035"
    gmat_elems_line = "25198.05865856482         7191.999999999875         0.02000000000000014       12.84999999999999         306.6                      314.1899999999973          99.8870000000026         "

    rv_init = parse_rv_line(rv_line, ITRS)
    gmat_elems = _parse_elems_line(rv_init.obstime, gmat_elems_line, _EARTH_GMAT)

    return (
        rv_init,
        gmat_elems,
        [
            9 * u.mm,  # sm axis
            1e-9,  # ecc
            5e-7 * u.deg,  # inc
            5e-6 * u.deg,  # raan
            5e-6 * u.deg,  # argp
            3e-6 * u.deg,  # true an
        ],
    )


def __init_ell_equator():
    """Case 2: Initialises an Earth Elliptical Equatorial orbit satellite. GMAT example."""

    rv_line = "2010-01-01T13:24:28.100   7.213392947764267e+03   8.523654531348812e+01  -2.783146976770290e-16   5.902225938368851e-02   7.421779936019859e+00   1.595360086373873e-18"
    gmat_elems_line = "25198.05865856482         7191.999999999999         0.01999999999999993       0.0         0.0                      260.7899999999987          99.88700000000058         "

    rv_init = parse_rv_line(rv_line, GCRS)
    gmat_elems = _parse_elems_line(rv_init.obstime, gmat_elems_line, _EARTH_GMAT)

    return (
        rv_init,
        gmat_elems,
        [
            0.0001 * u.mm,  # sm axis
            1e-14,  # ecc
            1e-14 * u.deg,  # inc
            1e-14 * u.deg,  # raan
            5e-11 * u.deg,  # argp
            3e-11 * u.deg,  # true an
        ],
    )


def __init_circ_incl():
    """Case 3: Initialises an Earth Circular Inclined orbit satellite. GMAT example."""

    rv_line = "2010-01-01T13:24:28.100   7.074397805771428e+03  -1.988417626249246e+00   1.295282105138377e+03  -1.757812701818108e-01   7.378906649083491e+00   9.713860595480043e-01"
    gmat_elems_line = "25198.05865856482         7191.999999999996         2.02792751048046e-16      12.84999999999999         306.6                      0                         54.07699999999998         "

    rv_init = parse_rv_line(rv_line, GCRS)
    gmat_elems = _parse_elems_line(rv_init.obstime, gmat_elems_line, _EARTH_GMAT)

    return (
        rv_init,
        gmat_elems,
        [
            0.0001 * u.mm,  # sm axis
            2e-14,  # ecc
            3e-13 * u.deg,  # inc
            5e-12 * u.deg,  # raan
            5e-13 * u.deg,  # argp
            2e-12 * u.deg,  # true an
        ],
    )


def __init_circ_equator():
    """Case 4: Initialises an Earth Circular Equatorial orbit satellite. GMAT example."""

    rv_line = "2010-01-01T13:24:28.100   4.355971648857403e+03  -5.722794334444523e+03  -1.411143337598508e-15   5.923828927156149e+00   4.508991473634389e+00   7.270719060994101e-19"
    gmat_elems_line = "25198.05865856482         7191.999999999996         2.02792751048046e-16      0.0         0.0                      0                         307.2770000000005     "

    rv_init = parse_rv_line(rv_line, GCRS)
    gmat_elems = _parse_elems_line(rv_init.obstime, gmat_elems_line, _EARTH_GMAT)

    return (
        rv_init,
        gmat_elems,
        [
            0.0002 * u.mm,  # sm axis
            3e-14,  # ecc
            1e-14 * u.deg,  # inc
            1e-14 * u.deg,  # raan
            1e-13 * u.deg,  # argp
            1e-13 * u.deg,  # true an
        ],
    )


def __init_hyperbolic_incl():
    """Initialises an Earth Hyperbolic Inclined orbit satellite. GMAT example."""

    rv_line = "2010-01-01T13:24:28.100   3.464897862367937e+02  -9.738869896647320e-02   6.344031422145902e+01   3.588822034885639e+01   3.027151368726605e+01   1.068941616196268e+01"
    gmat_elems_line = "25198.05865856482         -7191.999999999876         1.02       12.84999999999999         306.6                      314.1899999999976          99.88700000000226         "

    rv_init = parse_rv_line(rv_line, GCRS)
    gmat_elems = _parse_elems_line(rv_init.obstime, gmat_elems_line, _EARTH_GMAT)

    return (
        rv_init,
        gmat_elems,
        [
            0.0002 * u.mm,  # sm axis
            1e-14,  # ecc
            1e-14 * u.deg,  # inc
            5e-13 * u.deg,  # raan
            1e-13 * u.deg,  # argp
            1e-13 * u.deg,  # true an
        ],
    )


@pytest.mark.parametrize(
    "init_case",
    [
        "__init_ell_incl",
        "__init_ell_equator",
        "__init_circ_incl",
        "__init_circ_equator",
        "__init_hyperbolic_incl",
    ],
)
def test_rv_to_osc_elems(init_case):
    """Tests the various orbit cases from rv to elems."""

    # find the function and run it
    rv_gmat, gmat_elems, errs = globals()[init_case]()

    osc_elems = OsculatingKeplerianOrbElems.from_cartesian(rv_gmat, _EARTH_GMAT)

    _elems_assert(osc_elems, gmat_elems, errs)


@pytest.mark.parametrize(
    "init_case",
    [
        "__init_ell_equator",
        "__init_circ_incl",
        "__init_circ_equator",
        "__init_hyperbolic_incl",
    ],
)
def test_osc_elems_to_rv(init_case):
    """Tests the various orbit cases from elems to rv."""

    # find the function and run it
    rv_gmat, gmat_elems, errs = globals()[init_case]()

    rv = gmat_elems.to_cartesian()

    r_diff = pos_err(rv, rv_gmat)
    v_diff = vel_err(rv, rv_gmat)

    print(f"r diff      :  {r_diff}")
    print(f"v diff      :  {v_diff}")

    allowable_pos_diff = 5e-05 * u.mm
    allowable_vel_diff = 5e-06 * u.mm / u.s

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0


def test_parabolic(init_elems_leo):
    """Tests eccentricity setter with input value = 1
    - should raise `ValueError`."""
    with pytest.raises(ValueError):
        orb_elems = init_elems_leo
        orb_elems.eccentricity = 1.0000000001


def test_ecc_out_of_bounds(init_elems_leo):
    """Tests eccentricity setter with input value out of bounds
    - should raise `ValueError`."""
    with pytest.raises(ValueError):
        orb_elems = init_elems_leo
        orb_elems.eccentricity = -1.2


def test_incl_out_of_bounds(init_elems_leo):
    """Tests inclination setter with input value out of bounds
    - should raise `ValueError`."""
    with pytest.raises(ValueError):
        orb_elems = init_elems_leo
        orb_elems.inclination = 210 * u.deg


def test_raan_out_of_bounds(init_elems_leo):
    """Tests RAAN setter with input value out of bounds."""
    orb_elems = init_elems_leo
    orb_elems.raan = 390 * u.deg
    assert orb_elems.raan.to_value(u.deg) == approx((30 * u.deg).to_value(), abs=1e-8)


def test_getters_setters(init_elems_leo):
    """Test getters and setters of the Orbit."""
    orb_elems = init_elems_leo

    # test getters
    assert orb_elems.sm_axis.to_value() == approx(
        (7055.95320378 * u.km).to_value(), abs=(0.01 * u.mm).to_value(u.km)
    )

    assert orb_elems.period.to_value(u.s) == approx(
        (5898.539855583115 * u.s).to_value(), abs=(1 * u.us).to_value(u.s)
    )
    assert orb_elems.inclination.to_value(u.deg) == approx(
        (1.71234602 * u.rad).to_value(u.deg), abs=1e-6
    )
    assert orb_elems.raan.to_value(u.deg) == approx(
        (4.42402394 * u.rad).to_value(u.deg), abs=1e-8
    )
    assert orb_elems.eccentricity == approx(0.0020835, abs=1e-9)
    assert orb_elems.true_anomaly.to_value(u.deg) == approx(
        (1.5 * u.rad).to_value(u.deg), abs=1e-8
    )
    assert orb_elems.arg_periapsis.to_value(u.deg) == approx(
        (5.23982923 * u.rad).to_value(u.deg), abs=1e-8
    )

    # test setters
    orb_elems.inclination = 20.00 * u.deg
    assert orb_elems.inclination.to_value(u.deg) == approx(
        (20.00 * u.deg).to_value(), abs=1e-8
    )

    orb_elems.raan = 210.00 * u.deg
    assert orb_elems.raan.to_value(u.deg) == approx(
        (210.00 * u.deg).to_value(), abs=1e-8
    )

    orb_elems.arg_periapsis = 30.00 * u.deg
    assert orb_elems.arg_periapsis.to_value(u.deg) == approx(
        (30.00 * u.deg).to_value(), abs=1e-8
    )

    orb_elems.true_anomaly = 210.00 * u.deg
    assert orb_elems.true_anomaly.to_value(u.deg) == approx(
        (210.00 * u.deg).to_value(), abs=1e-8
    )

    orb_elems.eccentricity = 0.051
    assert orb_elems.eccentricity == approx(0.051, abs=1e-10)

    orb_elems.eccentricity = 0.051 * u.dimensionless_unscaled
    assert orb_elems.eccentricity == approx(0.051, abs=1e-10)

    orb_elems.sm_axis = 7000000 * u.m
    assert orb_elems.sm_axis == 7000 * u.km
