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
from satmad.core.occultation import compute_occultation
from satmad.propagation.classical_orb_elems import OsculatingKeplerianOrbElems
from satmad.propagation.numerical_propagators import NumericalPropagator
from satmad.utils.discrete_time_events import DiscreteTimeEvents
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


def test_occultation_times():
    """Tests the umbra and penumbra times against GMAT.

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

    time_list = trajectory.coord_list.obstime

    # init timer
    begin = time.time()

    # init interpolated planet positions
    # this is 10-15% faster than the list comprehension
    occult_pos_list = occult_traj(time_list)
    illum_pos_list = illum_traj(time_list)

    occultation_results = [
        compute_occultation(
            coord,
            occult_pos_list[i],
            illum_pos_list[i],
            occulting_body=EARTH,
            illum_body=SUN,
        )
        for i, coord in enumerate(trajectory.coord_list)
    ]

    # end timer
    end = time.time()
    if output_timer_results:
        print(f"Occultation finding: {end - begin} seconds")

    umbra_params = np.asarray(
        [result[3].to_value(u.km) for result in occultation_results]
    )
    penumbra_params = np.asarray(
        [result[2].to_value(u.km) for result in occultation_results]
    )

    # # plot umbra params if required
    # from satmad.plots.basic_plots import plot_time_param
    #
    # plot_time_param(time_list, umbra_params)

    # ------------------- check umbra times -------------

    # Find intervals in data
    umbra_intervals = DiscreteTimeEvents(
        time_list, umbra_params, 0.0, neg_to_pos_is_start=False
    ).start_end_intervals

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

    # Find intervals in data
    penumbra_events = DiscreteTimeEvents(
        time_list, penumbra_params, 0.0, neg_to_pos_is_start=False
    )
    # this nominally includes penumbra and umbra times. Subtract the umbra times.
    penumbra_intervals = umbra_intervals.invert().intersect_list(
        penumbra_events.start_end_intervals
    )

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
