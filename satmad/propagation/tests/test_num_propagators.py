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
    SkyCoord,
)
from astropy.time import Time

from satmad.propagation.numerical_propagators import ODESolverType
from satmad.propagation.tests.num_prop_analysis.num_prop_analysis_engine import (
    energy_along_trajectory,
    propagation_engine,
)


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
    rv_init = SkyCoord(
        r_gcrs.with_differentials(v_gcrs),
        obstime=time,
        frame=GCRS,
        representation_type="cartesian",
        differential_type="cartesian",
    )

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
