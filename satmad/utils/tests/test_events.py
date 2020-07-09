# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Tests the `discrete_time_intervals` package functionalities.

"""
import pytest
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation
from astropy.time import Time, TimeDelta

from satmad.coordinates.tests.test_trajectory import (
    generate_timelist,
    propagation_engine,
)
from satmad.propagation.tle import TLE
from satmad.utils.discrete_time_events import DiscreteTimeEvents


def prepare_timelist(init_time, duration, steps):
    """Prepares time list setup."""

    dt, time_list = generate_timelist(init_time, duration, steps)

    return {"dt": dt, "time_list": time_list}


def generate_leo_trajectory(time_list):
    """Generates the LEO trajectory."""

    # Init TLE
    line1 = "1 37791U 11044D   15230.18696540  .00000191  00000-0  45093-4 0  9993"
    line2 = "2 37791 098.1618 317.9053 0023092 067.3069 293.0602 14.64400558213855"

    trajectory = propagation_engine(line1, line2, time_list)["traj"]

    return trajectory


def generate_geo_42e_trajectory(time_list):
    """Generates a GEO (42 deg East) trajectory."""

    # Init TLE
    tle = TLE.init_geo(Time("2015-10-04T00:00:00.000", scale="utc"), 42 * u.deg)
    line1, line2 = tle.export_tle()

    trajectory = propagation_engine(line1, line2, time_list)["traj"]

    return trajectory


def generate_geo_140w_trajectory(time_list):
    """Generates a GEO (140 deg West) trajectory."""

    # Init TLE
    tle = TLE.init_geo(Time("2015-10-04T00:00:00.000", scale="utc"), -140 * u.deg)
    line1, line2 = tle.export_tle()

    trajectory = propagation_engine(line1, line2, time_list)["traj"]

    return trajectory


def event_engine(trajectory_generator, init_time, duration, steps):
    """Generates events list to be tested."""
    time_list = prepare_timelist(init_time, duration, steps)["time_list"]
    trajectory = trajectory_generator(time_list)

    # Init location - default Ellipsoid is WGS84
    ground_loc = EarthLocation(lat=41.015137, lon=28.979530, height=100 * u.m)

    # generate the alt-az list
    sat_alt_az_list = trajectory.coord_list.transform_to(
        AltAz(location=ground_loc, obstime=time_list)
    )

    # Find intervals in data
    events = DiscreteTimeEvents(time_list, sat_alt_az_list.alt.deg, 5)

    return events


@pytest.mark.parametrize(
    "steps, time_error, ang_error_1, ang_error_2",
    [
        (4000, 0.13 * u.ms, 9.0e-7 * u.deg, 0.02 * u.deg),
        (2000, 2.0 * u.ms, 1.6e-4 * u.deg, 0.12 * u.deg),
        (1000, 160 * u.ms, 2.1e-3 * u.deg, 3.80 * u.deg),
    ],
)
def test_basic_pass_find(steps, time_error, ang_error_1, ang_error_2):
    """Tests the basic event finding functionality."""
    # test init time
    init_time = Time("2015-10-04T00:00:00.000", scale="utc")

    duration = TimeDelta(1.0, format="jd")

    events = event_engine(generate_leo_trajectory, init_time, duration, steps)

    # print(events.start_end_intervals)
    # print(events.max_min_table)

    # Check max events
    # ----------------
    # Case 1: "2015-10-04T11:51:44.967501459", alt: 6.914915588508499 deg
    diff = (events.max_min_table[2]["value"] - 6.914915588508499) * u.deg

    # print(f"diff : {diff} --- (stepsize: {(duration/steps).to(u.s)})")
    # print(
    #     f'{events.max_min_table[2]["time"].datetime64} **** 2015-10-04T11:51:44.967501459'
    # )

    assert pytest.approx(diff.to_value(u.deg), abs=ang_error_1.to_value(u.deg)) == 0.0

    # Case 2: "2015-10-04T10:15:28.951945200", alt: 63.38630365664624 deg
    diff = (events.max_min_table[1]["value"] - 63.38630365664624) * u.deg

    # print(f"diff : {diff} --- (stepsize: {(duration/steps).to(u.s)})")
    # print(
    #     f'{events.max_min_table[1]["time"].datetime64} **** 2015-10-04T10:15:28.951945200'
    # )
    # time error:
    # - exact match for 21 sec timestep
    # - 2.5 sec off for 43 sec timestep
    # - 6.2 sec off for 86 sec timestep

    assert pytest.approx(diff.to_value(u.deg), abs=ang_error_2.to_value(u.deg)) == 0.0

    # Check start / end intervals
    # ---------------------------
    diff = (
        events.start_end_intervals.get_interval(2).start
        - Time("2015-10-04T11:49:36.742030377")
    ).to(u.ms)

    # print(f"diff : {diff} --- (stepsize: {(duration/steps).to(u.s)})")

    assert pytest.approx(diff.to_value(u.ms), abs=time_error.to_value(u.ms)) == 0.0

    diff = (
        events.start_end_intervals.get_interval(4).end
        - Time("2015-10-04T21:23:19.417927488")
    ).to(u.ms)

    # print(f"diff : {diff} --- (stepsize: {(duration/steps).to(u.s)})")

    assert pytest.approx(diff.to_value(u.ms), abs=time_error.to_value(u.ms)) == 0.0


def test_borderline():
    """Tests the borderline cases."""
    # test init time
    init_time = Time("2015-10-04T08:36:00.000", scale="utc")

    duration = Time("2015-10-04T21:20:00.000", scale="utc") - init_time

    steps = 4000

    events = event_engine(generate_leo_trajectory, init_time, duration, steps)

    print(f"search: {events.search_interval}")
    print(events.start_end_intervals)

    assert (
        str(events.start_end_intervals.get_interval(0).start)
        == "2015-10-04T08:36:00.000"
    )

    assert (
        str(events.start_end_intervals.get_interval(4).end) == "2015-10-04T21:19:48.540"
    )


def test_visibility():
    """Tests the GEO cases with continuous visibility and no visibility."""
    # test init time
    init_time = Time("2015-10-04T08:00:00.000", scale="utc")

    duration = 1 * u.day

    steps = 4000

    # GEO with full visibility
    events = event_engine(generate_geo_42e_trajectory, init_time, duration, steps)

    print(f"full vis: {events.start_end_intervals}")

    assert (
        str(events.start_end_intervals.get_interval(0).start)
        == "2015-10-04T08:00:00.000"
    )

    assert (
        str(events.start_end_intervals.get_interval(0).end) == "2015-10-05T07:59:38.400"
    )

    # GEO with no visibility
    events = event_engine(generate_geo_140w_trajectory, init_time, duration, steps)

    print(f"no vis: {events.start_end_intervals}")

    assert not events.start_end_intervals.intervals


def test_init_array_mismatch():
    """Tests the input error due to array size mismatch."""
    with pytest.raises(ValueError):
        time_list = prepare_timelist(
            Time("2015-10-04T00:00:00.000", scale="utc"), 1.0 * u.day, 1000
        )["time_list"]
        DiscreteTimeEvents(time_list, [0.1, 1.2, 2.3])
