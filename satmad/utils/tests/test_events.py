# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
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
from satmad.core.celestial_bodies_lib import EARTH
from satmad.propagation.classical_orb_elems import OsculatingKeplerianOrbElems
from satmad.propagation.numerical_propagators import NumericalPropagator
from satmad.propagation.tle import TLE
from satmad.utils.discrete_time_events import DiscreteTimeEvents
from satmad.utils.timeinterval import TimeInterval


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


def check_instantaneous(trajectory_generator, init_time, duration, steps, time):
    """Generates az-el for a single time."""
    time_list = prepare_timelist(init_time, duration, steps)["time_list"]
    trajectory = trajectory_generator(time_list)

    # Init location - default Ellipsoid is WGS84
    ground_loc = EarthLocation(lat=41.015137, lon=28.979530, height=100 * u.m)

    # generate the alt-az
    sat_alt_az = trajectory(time).transform_to(AltAz(location=ground_loc, obstime=time))

    return sat_alt_az


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

    print(events.start_end_intervals)
    print(events.max_min_table)

    # Check max events
    # ----------------
    # Case 1: "2015-10-04T11:51:44.965499023", alt: 6.921085405952291 deg
    diff = (events.max_min_table[2]["value"] - 6.921085405952291) * u.deg

    print(f"diff : {diff} --- (stepsize: {(duration/steps).to(u.s)})")
    print(f'at time : {events.max_min_table[2]["time"].datetime64}')

    assert pytest.approx(diff.to_value(u.deg), abs=ang_error_1.to_value(u.deg)) == 0.0

    # Case 2: "2015-10-04T10:15:28.948189027", alt: 63.427834618376515 deg
    diff = (events.max_min_table[1]["value"] - 63.427834618376515) * u.deg

    # print(f"diff : {diff} --- (stepsize: {(duration/steps).to(u.s)})")
    # print(f'at time : {events.max_min_table[1]["time"].datetime64}')
    # time error:
    # - exact match for 21 sec timestep
    # - 2.5 sec off for 43 sec timestep
    # - 6.2 sec off for 86 sec timestep

    # altaz = check_instantaneous(
    #     generate_leo_trajectory,
    #     init_time,
    #     duration,
    #     steps,
    #     events.max_min_table[2]["time"],
    # )
    assert pytest.approx(diff.to_value(u.deg), abs=ang_error_2.to_value(u.deg)) == 0.0

    # Check start / end intervals
    # ---------------------------
    diff = (
        events.start_end_intervals.get_interval(2).start
        - Time("2015-10-04T11:49:36.549018004")
    ).to(u.ms)

    # print(f"diff : {diff} --- (stepsize: {(duration/steps).to(u.s)})")

    assert pytest.approx(diff.to_value(u.ms), abs=time_error.to_value(u.ms)) == 0.0

    diff = (
        events.start_end_intervals.get_interval(4).end
        - Time("2015-10-04T21:23:19.43270995")
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
        DiscreteTimeEvents(time_list, [0.1, 1.2, 2.3], crossing_value=1)


def test_min_max_event():
    """Tests the min max event finding with apoapsis and periapsis."""

    time = Time("2020-01-11T11:00:00.000", scale="utc")
    central_body = EARTH

    sm_axis = 7056.0 * u.km
    ecc = 0.02 * u.dimensionless_unscaled
    incl = 0 * u.deg
    raan = 0 * u.deg
    arg_perigee = 90 * u.deg
    true_an = 20 * u.deg

    init_orb_elems = OsculatingKeplerianOrbElems(
        time, sm_axis, ecc, incl, raan, arg_perigee, true_an, central_body
    )

    # generate cartesian initial conditions
    init_pvt = init_orb_elems.to_cartesian()

    # Set up propagation config
    stepsize = 10 * u.s

    prop_start = init_pvt.obstime
    prop_duration = init_orb_elems.period

    # init propagator with defaults - run propagation and get trajectory
    trajectory = NumericalPropagator(stepsize).gen_trajectory(
        init_pvt, TimeInterval(prop_start, prop_duration)
    )

    # Extract search range
    time_list = trajectory.coord_list.obstime
    r_list = trajectory.coord_list.cartesian.without_differentials().norm()

    # Find time events
    events = DiscreteTimeEvents(time_list, r_list)

    # Min / Max Event times
    # print(events.max_min_table)

    # Cross-check with orbital elems
    # print("\nCheck with initial conditions:")
    # print(f"Apoapsis : {init_orb_elems.apoapsis}")
    # print(f"Periapsis: {init_orb_elems.periapsis}")

    assert (
        pytest.approx(
            (events.max_min_table[0]["value"] - init_orb_elems.apoapsis).to_value(u.mm),
            abs=(1e-3 * u.mm).to_value(u.mm),
        )
        == 0.0
    )

    assert (
        pytest.approx(
            (events.max_min_table[1]["value"] - init_orb_elems.periapsis).to_value(
                u.mm
            ),
            abs=(5e-3 * u.mm).to_value(u.mm),
        )
        == 0.0
    )

    assert events.max_min_table[0]["type"] == "max"
    assert events.max_min_table[1]["type"] == "min"
