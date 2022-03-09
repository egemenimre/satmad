# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Tests associated with the numerical propagators.

"""
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianDifferential,
    CartesianRepresentation,
)
from astropy.time import Time
from pytest import approx

from satmad.coordinates.frames import init_pvt
from satmad.coordinates.tests.test_baseframe_transform import pos_err, vel_err
from satmad.core.celestial_bodies_lib import MoonCRS
from satmad.core.celestial_body import CelestialBody, CelestialBodyEllipsoid
from satmad.propagation.numerical_propagators import ODESolverType
from satmad.propagation.tests.num_prop_analysis_engine import (
    energy_along_trajectory,
    propagation_engine,
)


@pytest.fixture
def rv_init_moon_crs():
    """Initialises a Moon low-orbit satellite."""
    time = Time("2020-01-01T11:00:00.000", scale="utc")

    # GMAT test case
    v_moon_crs = CartesianDifferential([1, -1, 0.6], unit=u.km / u.s)
    r_moon_crs = CartesianRepresentation([1000, 1000, 2000], unit=u.km)
    rv_init = init_pvt(MoonCRS, time, r_moon_crs.with_differentials(v_moon_crs))

    return rv_init


_MOONGMAT = CelestialBody(
    "Moon",
    "Default Moon Model. ",
    4.90280105600000e12 * u.m**3 / u.s**2,
    ellipsoid=CelestialBodyEllipsoid(
        "Moon Ellipsoid as defined by GMAT",
        1738.2 * u.km,
        np.inf * u.dimensionless_unscaled,
    ),
    inert_coord=MoonCRS,
)
"""This is the Moon as defined by GMAT"""


def test_num_propagator_moon(rv_init_moon_crs):
    """This tests the final position and velocity against GMAT using a pure Keplerian
    5 day propagation.

    GMAT Propagator is Dormand-Prince 853, with an Accuracy of 1e-13."""

    allowable_pos_diff = 13 * u.mm
    allowable_vel_diff = 0.0055 * u.mm / u.s

    # Set up propagation config
    stepsize = 60 * u.s
    solver_type = ODESolverType.DOP853
    rtol = 1e-13
    atol = 1e-15
    init_time_offset = 0.0 * u.day
    duration = 10.00 * u.day

    # run propagation and get trajectory
    trajectory = propagation_engine(
        rv_init_moon_crs,
        stepsize,
        solver_type,
        init_time_offset,
        duration,
        rtol,
        atol,
        central_body=_MOONGMAT,
    )
    # no interpolation, just take final element
    rv_fin = trajectory.coord_list[-1]

    # GMAT test case - final point
    time = Time("2020-01-11T11:00:00.000", scale="utc")

    v_moon_crs_gmat = CartesianDifferential(
        [6.159801708156309e-01, -1.128605958878207e00, 1.075005007078085e-02],
        unit=u.km / u.s,
    )
    r_moon_crs_gmat = CartesianRepresentation(
        [1.870952777669765e03, -1.811234513779054e02, 2.305452195816935e03],
        unit=u.km,
    )
    rv_fin_gmat = init_pvt(
        MoonCRS, time, r_moon_crs_gmat.with_differentials(v_moon_crs_gmat)
    )

    r_diff = pos_err(rv_fin, rv_fin_gmat)
    v_diff = vel_err(rv_fin, rv_fin_gmat)

    print(f"r diff      :  {r_diff}")
    print(f"v diff      :  {v_diff}")

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0


@pytest.fixture
def rv_init_leo_gcrs():
    """Initialises an Earth LEO satellite."""
    time = Time("2004-04-06T07:51:28.386009", scale="utc")

    # Vallado IAU 2000 - Table 3-6
    v_gcrs = CartesianDifferential(
        [-4.7432201610, 0.7905364950, 5.5337557240], unit=u.km / u.s
    )
    r_gcrs = CartesianRepresentation(
        [5102.50895290, 6123.01139910, 6378.13693380], unit=u.km
    )

    rv_init = init_pvt(GCRS, time, r_gcrs.with_differentials(v_gcrs))
    return rv_init


def test_num_propagator(rv_init_leo_gcrs):
    """Tests the basic numerical propagation accuracy in terms of energy."""

    # Set up propagation config
    # Start with ITRS to test coord transform on propagator entry
    rv_init = rv_init_leo_gcrs.transform_to(ITRS)
    stepsize = 60 * u.s
    solver_type = ODESolverType.DOP853
    rtol = 1e-13
    atol = 1e-15
    init_time_offset = 0.5 * u.day
    duration = 10.001 * u.day

    # run propagation and get trajectory
    trajectory = propagation_engine(
        rv_init, stepsize, solver_type, init_time_offset, duration, rtol, atol
    )
    # print(trajectory.coord_list[-1])

    # Compute the energy along the trajectory
    time_list, energy_list, energy_diff_list = energy_along_trajectory(
        trajectory.coord_list
    )

    abs_energy_diff_list = abs(np.array(energy_diff_list))

    # Max energy diff (abs value) - should be low
    energy_diff_max = abs_energy_diff_list.max()
    # Mean energy diff - should be zero to have zero mean energy
    energy_diff_mean = np.array(energy_diff_list).mean()

    # initial energy diff and final energy diff to assess energy leak
    # First and last 100 elements taken and they should also ideally be zero
    init_energy_diff_mean = np.array(energy_diff_list)[:100].mean()
    final_energy_diff_mean = np.array(energy_diff_list)[-100:].mean()

    # print(f"max energy diff: {energy_diff_max}")
    # print(f"mean energy    : {energy_diff_mean}")
    #
    # print(f"init mean energy diff : {init_energy_diff_mean}")
    # print(f"final mean energy diff: {final_energy_diff_mean}")

    assert energy_diff_max < 4.50e-11
    assert energy_diff_mean < 2.95e-12

    assert init_energy_diff_mean < 1.25e-12
    assert final_energy_diff_mean < 5.20e-12
