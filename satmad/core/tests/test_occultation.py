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
from astropy.coordinates import (
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord,
    get_body_barycentric,
    get_body_barycentric_posvel,
    get_sun,
)
from astropy.time import Time
from pytest import approx

from satmad.coordinates.frames import MoonCRS
from satmad.coordinates.trajectory import Trajectory
from satmad.core.celestial_bodies import EARTH, MOON, SUN
from satmad.core.occultation import (
    multi_body_occultation_intervals,
    occultation_intervals,
    compute_occultation,
)
from satmad.propagation.classical_orb_elems import OsculatingKeplerianOrbElems
from satmad.propagation.numerical_propagators import NumericalPropagator
from satmad.utils.timeinterval import TimeInterval


def _init_orbit_luna():
    """Initialises the test orbit in a weird Lunar orbit."""
    obstime = Time("2020-01-10T11:30:00", scale="utc")

    v_luna = CartesianDifferential(
        [-1.061195314530757e00, 1.038484649573483e-01, 3.935374917002971e-02],
        unit=u.km / u.s,
    )
    r_luna = CartesianRepresentation(
        [3.001889915965003e02, 4.477826526299235e03, 6.101271598017775e-01], unit=u.km
    )
    rv_luna = SkyCoord(
        r_luna.with_differentials(v_luna),
        obstime=obstime,
        frame=MoonCRS,
        representation_type="cartesian",
        differential_type="cartesian",
    )

    return rv_luna


def _init_orbit_leo():
    """Initialises the test orbit in LEO."""
    obstime = Time("2020-01-01T11:30:00", scale="utc")
    central_body = EARTH

    # Initialise a near-polar orbit
    sm_axis = 7191.9 * u.km
    ecc = 0.02 * u.dimensionless_unscaled
    incl = 98.0 * u.deg
    raan = 306.6 * u.deg
    arg_perigee = 314.1 * u.deg
    true_an = 100.3 * u.deg

    orb_elems = OsculatingKeplerianOrbElems(
        obstime, sm_axis, ecc, incl, raan, arg_perigee, true_an, central_body
    )

    # convert to cartesian coordinates - this is the initial condition for the propagation
    init_pvt = orb_elems.to_cartesian()

    return init_pvt


def _init_trajectory(pvt0, stepsize, prop_interval, central_body=EARTH):
    """Initialises the trajectory."""

    # init propagator with defaults
    # run propagation and get trajectory
    trajectory = NumericalPropagator(
        stepsize, atol=1e-13, central_body=central_body
    ).gen_trajectory(pvt0, prop_interval)

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


def test_multi_body_occultation_intervals():
    """Tests the umbra and penumbra intervals against GMAT, where the occulting
    body is Earth and Moon.

    Using a stepsize of 60 seconds gives more points to evaluate and increases the
    accuracy of the entry-exit times by a few milliseconds.
    """

    output_timer_results = True

    # init timer
    begin = time.time()

    # Init trajectory
    pvt0 = _init_orbit_luna()

    # Set up propagation config
    stepsize = 120 * u.s
    prop_interval = TimeInterval(pvt0.obstime, 2.0 * u.day)

    # init propagator with defaults
    # run propagation and get trajectory
    trajectory = _init_trajectory(pvt0, stepsize, prop_interval, central_body=MOON)

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

    occult_traj_earth = Trajectory(
        SkyCoord(
            get_body_barycentric(EARTH.name, sparse_time_list, ephemeris="jpl"),
            obstime=sparse_time_list,
            frame="icrs",
            representation_type="cartesian",
            differential_type="cartesian",
        ).transform_to("gcrs")
    )

    r_moon, v_moon = get_body_barycentric_posvel(
        MOON.name, sparse_time_list, ephemeris="jpl"
    )
    v_moon = CartesianDifferential(v_moon.xyz)
    r_moon = r_moon.with_differentials(v_moon)

    occult_traj_moon = Trajectory(
        SkyCoord(
            r_moon.with_differentials(v_moon),
            obstime=sparse_time_list,
            frame="icrs",
            representation_type="cartesian",
            differential_type="cartesian",
        ).transform_to(MoonCRS)
    )

    # end timer
    end = time.time()
    if output_timer_results:
        print(f"Propagation and interpolations: {end - begin} seconds")

    # init timer
    begin = time.time()

    occult_bodies = {EARTH: occult_traj_earth, MOON: occult_traj_moon}
    # occult_bodies = {MOON: occult_traj_moon}
    # occult_bodies = {EARTH: occult_traj_earth}

    # Compute occultation intervals
    output_dict = multi_body_occultation_intervals(
        trajectory, occult_bodies, illum_traj, illum_body=SUN
    )

    # end timer
    end = time.time()
    if output_timer_results:
        print(f"Intervals finding: {end - begin} seconds")

    umbra_intervals_earth, penumbra_intervals_earth = output_dict[EARTH.name]
    umbra_intervals_moon, penumbra_intervals_moon = output_dict[MOON.name]

    # Start Time (UTC)            Stop Time (UTC)               Duration (s)    Occ Body        Type        Event Number  Total Duration (s)
    # 10 Jan 2020 11:52:56.832    10 Jan 2020 12:09:45.608      1008.7751815    Luna            Umbra       1             1008.7751815
    # 10 Jan 2020 19:52:48.654    10 Jan 2020 20:10:05.113      1036.4593227    Luna            Umbra       2             1036.4593227
    # 11 Jan 2020 03:52:40.660    11 Jan 2020 04:10:24.428      1063.7679860    Luna            Umbra       3             1063.7679860
    # 11 Jan 2020 11:52:32.833    11 Jan 2020 12:10:43.557      1090.7238229    Luna            Umbra       4             1090.7238229
    # 11 Jan 2020 19:52:25.159    11 Jan 2020 20:11:02.505      1117.3465637    Luna            Umbra       5             1117.3465637
    # 12 Jan 2020 03:52:17.623    12 Jan 2020 04:11:21.276      1143.6534295    Luna            Umbra       6             1143.6534295

    # Start Time (UTC)            Stop Time (UTC)               Duration (s)    Occ Body        Type        Event Number  Total Duration (s)
    # 10 Jan 2020 11:50:49.031    10 Jan 2020 11:52:56.832      127.80161456    Luna            Penumbra    1             127.80161456
    # 10 Jan 2020 12:09:45.608    10 Jan 2020 12:11:54.934      129.32634068    Luna            Penumbra    2             129.32634068
    # 10 Jan 2020 18:41:49.529    10 Jan 2020 20:26:28.746      6279.2173116    Earth           Penumbra    3             6279.2173116
    # 10 Jan 2020 19:50:43.540    10 Jan 2020 19:52:48.654      125.11435262    Luna            Penumbra    3             6279.2173116
    # 10 Jan 2020 20:10:05.113    10 Jan 2020 20:12:11.758      126.64437286    Luna            Penumbra    3             6279.2173116
    # 11 Jan 2020 03:50:38.094    11 Jan 2020 03:52:40.660      122.56627186    Luna            Penumbra    4             122.56627186
    # 11 Jan 2020 04:10:24.428    11 Jan 2020 04:12:28.529      124.10153624    Luna            Penumbra    5             124.10153624
    # 11 Jan 2020 11:50:32.687    11 Jan 2020 11:52:32.833      120.14619472    Luna            Penumbra    6             120.14619472
    # 11 Jan 2020 12:10:43.557    11 Jan 2020 12:12:45.244      121.68665263    Luna            Penumbra    7             121.68665263
    # 11 Jan 2020 19:50:27.315    11 Jan 2020 19:52:25.159      117.84421628    Luna            Penumbra    8             117.84421628
    # 11 Jan 2020 20:11:02.505    11 Jan 2020 20:13:01.895      119.38981398    Luna            Penumbra    9             119.38981398
    # 12 Jan 2020 03:50:21.971    12 Jan 2020 03:52:17.623      115.65151513    Luna            Penumbra    10            115.65151513
    # 12 Jan 2020 04:11:21.276    12 Jan 2020 04:13:18.478      117.20219948    Luna            Penumbra    11            117.20219948

    # ------------------- check umbra times -------------

    print("Umbra Start-End Intervals - Earth:")
    print(umbra_intervals_earth)
    print("Umbra Start-End Intervals - Moon:")
    print(umbra_intervals_moon)

    # check Umbra params
    gmat_umbra_times_earth = []

    start_diff = 160 * u.ms
    end_diff = 80 * u.ms

    _check_entry_exit_times(
        umbra_intervals_earth, gmat_umbra_times_earth, start_diff, end_diff
    )

    gmat_umbra_times_moon = [
        (0, "2020-01-10T11:52:56.832", "2020-01-10T12:09:45.608"),
        (2, "2020-01-11T03:52:40.660", "2020-01-11T04:10:24.428"),
        (5, "2020-01-12T03:52:17.623", "2020-01-12T04:11:21.276"),
    ]

    start_diff = 140 * u.ms
    end_diff = 80 * u.ms

    _check_entry_exit_times(
        umbra_intervals_moon, gmat_umbra_times_moon, start_diff, end_diff
    )

    # ------------------- check penumbra times -------------

    print("Penumbra Start-End Intervals - Earth:")
    print(penumbra_intervals_earth)
    print("Penumbra Start-End Intervals - Moon:")
    print(penumbra_intervals_moon)

    # check Penumbra params
    gmat_penumbra_times_moon = [
        (0, "2020-01-10T11:50:49.031", "2020-01-10T11:52:56.832"),
        (2, "2020-01-10T19:50:43.540", "2020-01-10T19:52:48.654"),
        (9, "2020-01-11T20:11:02.505", "2020-01-11T20:13:01.895"),
    ]

    start_diff = 73 * u.ms
    end_diff = 120 * u.ms

    _check_entry_exit_times(
        penumbra_intervals_moon, gmat_penumbra_times_moon, start_diff, end_diff
    )

    gmat_penumbra_times_earth = [
        (0, "2020-01-10T18:41:49.529", "2020-01-10T20:26:28.746")
    ]

    start_diff = 45 * u.s
    end_diff = 40 * u.s

    _check_entry_exit_times(
        penumbra_intervals_earth, gmat_penumbra_times_earth, start_diff, end_diff
    )


def test_occultation_intervals():
    """Tests the umbra and penumbra intervals against GMAT, where the occulting
    body is Earth only.

    Using a stepsize of 60 seconds gives more points to evaluate and increases the
    accuracy of the entry-exit times by a few milliseconds.
    """

    output_timer_results = True

    # init timer
    begin = time.time()

    # Init trajectory
    pvt0 = _init_orbit_leo()

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

    # Compute occultation intervals
    (umbra_intervals, penumbra_intervals) = occultation_intervals(
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
