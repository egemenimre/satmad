# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Tests related to occultations, shadows and illumination.

"""
import time

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time
from pytest import approx

from satmad.coordinates.trajectory import Trajectory
from satmad.core.celestial_bodies import EARTH, SUN
from satmad.core.occultation import compute_single_body_occultation_times
from satmad.propagation.classical_orb_elems import OsculatingKeplerianOrbElems
from satmad.propagation.numerical_propagators import NumericalPropagator
from satmad.utils.timeinterval import TimeInterval


def _init_orbit():
    """Initialises the test orbit."""
    time = Time("2020-01-01T11:30:00", scale="utc")
    central_body = EARTH

    # Initialise a near-polar orbit
    sm_axis = 7191.9 * u.km
    ecc = 0.02 * u.dimensionless_unscaled
    incl = 98.0 * u.deg
    raan = 306.6 * u.deg
    arg_perigee = 314.1 * u.deg
    true_an = 100.3 * u.deg

    orb_elems = OsculatingKeplerianOrbElems(
        time, sm_axis, ecc, incl, raan, arg_perigee, true_an, central_body
    )

    # convert to cartesian coordinates - this is the initial condition for the propagation
    init_pvt = orb_elems.to_cartesian()

    return init_pvt


def _init_trajectory(pvt0, stepsize, prop_interval):
    """Initialises the trajectory."""

    # init propagator with defaults
    # run propagation and get trajectory
    trajectory = NumericalPropagator(stepsize).gen_trajectory(pvt0, prop_interval)

    return trajectory


# def test_time_error():
#     """ """
#     with pytest.raises(ValueError):
#         # Init PVT
#         pvt0 = _init_orbit()
#
#         r_illum = SkyCoord(
#             get_body_barycentric(SUN.name, pvt0.obstime, ephemeris="jpl"),
#             obstime=pvt0.obstime - 1 * u.s,
#             frame="icrs",
#             representation_type="cartesian",
#             differential_type="cartesian",
#         )
#
#         r_occult = SkyCoord(
#             get_body_barycentric(SUN.name, pvt0.obstime, ephemeris="jpl"),
#             obstime=pvt0.obstime,
#             frame="icrs",
#             representation_type="cartesian",
#             differential_type="cartesian",
#         )
#
#         compute_occultation(
#             pvt0,
#             r_occult,
#             r_illum,
#         )
#####################
# make sure all times match
# allowable_time_diff = 1 * u.ms
# if (
#     abs(r_illum_body.obstime - rv_obj.obstime) > allowable_time_diff
#     or abs(r_occult_body.obstime - rv_obj.obstime) > allowable_time_diff
# ):
#     raise ValueError(
#         f"Occultation calculation: Position vector times do not match. "
#         f"Illum - obj diff: {(r_illum_body.obstime - rv_obj.obstime).to(u.s)}, "
#         f"Occult - obj diff: {(r_occult_body.obstime - rv_obj.obstime).to(u.s)}."
#     )


def test_occultation_intervals():
    """Tests the umbra and penumbra intervals against GMAT.

    Using a stepsize of 60 seconds gives more points to evaluate and increases the
    accuracy of the entry-exit times by a few milliseconds.
    """

    output_timer_results = True

    # init timer
    begin = time.time()

    # Init trajectory
    pvt0 = _init_orbit()

    # Set up propagation config
    stepsize = 120 * u.s
    prop_interval = TimeInterval(pvt0.obstime, 2.0 * u.day)

    # init propagator with defaults
    # run propagation and get trajectory
    trajectory = _init_trajectory(pvt0, stepsize, prop_interval)

    # begin = time.time()
    sparse_stepsize = 2 * u.hr
    sparse_time_list = (
        pvt0.obstime + np.arange(-0.5, +2.5, sparse_stepsize.to_value(u.day)) * u.day
    )
    illum_traj = Trajectory(
        SkyCoord(
            get_body_barycentric(SUN.name, sparse_time_list, ephemeris="jpl"),
            obstime=sparse_time_list,
            frame="icrs",
            representation_type="cartesian",
            differential_type="cartesian",
        ).transform_to("gcrs")
    )

    occult_traj = Trajectory(
        SkyCoord(
            get_body_barycentric(EARTH.name, sparse_time_list, ephemeris="jpl"),
            obstime=sparse_time_list,
            frame="icrs",
            representation_type="cartesian",
            differential_type="cartesian",
        ).transform_to("gcrs")
    )

    # end timer
    end = time.time()
    if output_timer_results:
        print(f"Propagation and interpolations: {end - begin} seconds")

    # init timer
    begin = time.time()

    (
        occulting_body_name,
        umbra_intervals,
        penumbra_intervals,
    ) = compute_single_body_occultation_times(
        trajectory, occult_traj, illum_traj, occulting_body=EARTH, illum_body=SUN
    )

    # # plot umbra params if required
    # from satmad.plots.basic_plots import plot_time_param
    #
    # plot_time_param(time_list, umbra_params)

    # end timer
    end = time.time()
    if output_timer_results:
        print(f"Intervals finding: {end - begin} seconds")

    # ------------------- check umbra times -------------

    # print("Start-End Intervals:")
    # print(umbra_intervals)

    # check Umbra params
    gmat_umbra_times = [
        (0, "2020-01-01T11:43:39.068", "2020-01-01T12:16:44.181"),
        (2, "2020-01-01T15:05:58.227", "2020-01-01T15:39:04.946"),
        (6, "2020-01-01T21:50:36.564", "2020-01-01T22:23:46.456"),
        (15, "2020-01-02T13:01:02.907", "2020-01-02T13:34:19.753"),
        (22, "2020-01-03T00:49:10.141", "2020-01-03T01:22:32.221"),
    ]

    start_diff = 16 * u.ms
    end_diff = 50 * u.ms

    _check_entry_exit_times(umbra_intervals, gmat_umbra_times, start_diff, end_diff)

    # ------------------- check penumbra times -------------

    # print("Start-End Intervals:")
    # print(penumbra_intervals)

    # check Penumbra params
    gmat_penumbra_times = [
        (0, "2020-01-01T11:43:27.775", "2020-01-01T11:43:39.068"),
        (6, "2020-01-01T16:46:56.547", "2020-01-01T16:47:07.809"),
        (18, "2020-01-02T02:53:54.129", "2020-01-02T02:54:05.332"),
        (36, "2020-01-02T18:04:20.596", "2020-01-02T18:04:31.713"),
        (44, "2020-01-03T00:48:59.061", "2020-01-03T00:49:10.141"),
    ]

    start_diff = 73 * u.ms
    end_diff = 16 * u.ms

    _check_entry_exit_times(
        penumbra_intervals, gmat_penumbra_times, start_diff, end_diff
    )


def _check_entry_exit_times(intervals, truth_intervals, start_diff, end_diff):
    """Checks the `intervals` against the `truth_intervals` with the given
    tolerances."""

    for i, truth_event_entry, truth_event_exit in truth_intervals:

        truth_event = TimeInterval(Time(truth_event_entry), Time(truth_event_exit))
        interval = intervals.get_interval(i)

        # print(
        #     i,
        #     (interval.start - truth_event.start).to(u.ms),
        #     (interval.end - truth_event.end).to(u.ms),
        # )

        assert (interval.start - truth_event.start).to_value(u.s) == approx(
            0.0, abs=start_diff.to_value(u.s)
        )

        assert (interval.end - truth_event.end).to_value(u.s) == approx(
            0.0, abs=end_diff.to_value(u.s)
        )
