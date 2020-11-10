# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Trajectory Class tests.
"""


import numpy as np
from astropy import units as u
from astropy.coordinates import (
    TEME,
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord,
)
from astropy.time import Time, TimeDelta
from matplotlib import pyplot as plt
from sgp4.api import Satrec

from satmad.coordinates.trajectory import Trajectory
from satmad.propagation.tle import TLE


def generate_timelist(init_time, duration, steps):
    """generates the time lists for propagation."""
    # ****** Generate the discrete time instants through the propagation duration ******
    dt = duration / steps

    # Generate time list
    dt_list = dt * np.arange(0, steps, 1)
    time_list = init_time + dt_list

    return dt, time_list


def propagation_engine(line1, line2, time_list):
    """Tests the interpolation accuracy for a given TLE."""

    # Init satellite object from the TLE
    sat = Satrec.twoline2rv(line1, line2)

    # ****** Generate the pos, vel vectors for each time instant ******

    # Run the propagation and init pos and vel vectors in TEME
    e, r_list, v_list = sat.sgp4_array(time_list.jd1, time_list.jd2)

    # Load the time, pos, vel info into astropy objects (shallow copied)
    vel_list = CartesianDifferential(v_list, unit=u.km / u.s, xyz_axis=1)
    pos_list = CartesianRepresentation(
        r_list, unit=u.km, xyz_axis=1
    ).with_differentials(vel_list)

    # trajectory in astropy
    traj_astropy = SkyCoord(
        pos_list,
        obstime=time_list,
        frame=TEME.name,
        representation_type="cartesian",
        differential_type="cartesian",
    )

    # Init trajectory in Trajectory object
    trajectory = Trajectory(traj_astropy)

    return {"traj": trajectory, "sat": sat}


def check_errors(init_time, dt, sat, trajectory, generate_plots=False):
    """Checks the errors for the given trajectory."""
    # **** Check the results ****
    step_offset = 0
    test_stepsize = 1.0 * u.s
    test_end = 0.01 * u.day

    # print(f"Test duration: {(test_end - step_offset * dt).to(u.day)}")
    # print(f"Step size: {test_stepsize}")

    dt_list_hires = TimeDelta(
        np.arange(
            step_offset * dt.to_value(u.s),
            test_end.to_value(u.s),
            test_stepsize.to_value(u.s),
        ),
        format="sec",
    )
    time_list_hires = init_time + dt_list_hires

    # Generate the high-res test trajectory
    e_test_list, r_test_list, v_test_list = sat.sgp4_array(
        time_list_hires.jd1, time_list_hires.jd2
    )

    r_test_list = CartesianRepresentation(r_test_list.transpose(), unit=u.km)

    # Generate the error list
    r_err_list = (
        trajectory(time_list_hires).cartesian.without_differentials() - r_test_list
    )

    # Error stats
    max_err = r_err_list.norm().max()  # max expected error
    max_nominal_err = r_err_list[-100:].norm().max()  # max expected error nominal

    if generate_plots:
        print(f"max error at the beginning  : {max_err.to(u.mm)}")
        print(f"max error after settle down : {max_nominal_err.to(u.mm)}")

        dt_list_hires_plot = dt_list_hires.to_value(u.day)

        # Plot the error per axis
        fr_ticker = np.arange(
            dt_list_hires_plot[0], dt_list_hires_plot[-1], dt.to_value(u.day)
        )
        y_ticker = np.zeros(len(fr_ticker))
        plt.figure()
        plt.plot(
            dt_list_hires_plot,
            np.asarray(r_err_list.get_xyz(xyz_axis=1)),
            fr_ticker,
            y_ticker,
            "x",
        )
        plt.legend(["x", "y", "z", "stepsize"])
        plt.title(f"Interpolation Error\n({trajectory.interpolator_name})")
        plt.xlabel("Time [days]")
        plt.ylabel("Error [km]")
        # plt.yscale("log")  # Uncomment to get y axis in log scale
        plt.show()

    return max_err, max_nominal_err


def test_trajectory_export():
    """Tests the interpolation accuracy for ISS."""

    # Init TLE
    line1 = "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991"
    line2 = "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482"

    # test init time
    init_time = Time(2458826.3, format="jd", scale="utc")
    init_time.format = "isot"

    duration = TimeDelta(3.0, format="jd")
    steps = 4000  # number of steps within the propagation duration

    dt, time_list = generate_timelist(init_time, duration, steps)
    prop_output = propagation_engine(line1, line2, time_list)

    ts = prop_output["traj"].to_timeseries()

    pvt = prop_output["traj"].coord_list[10]

    # print(ts)
    # print(ts["time", "v_X"])

    assert (ts.time[10] - pvt.obstime) < 17 * u.us
    assert (ts[10]["r_X"] - pvt.cartesian.x) < 1 * u.um


def test_interpolation_err_iss():
    """Tests the interpolation accuracy for ISS."""

    # Init TLE
    line1 = "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991"
    line2 = "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482"

    # test init time
    init_time = Time(2458826.3, format="jd", scale="utc")

    duration = TimeDelta(3.0, format="jd")
    steps = 4000  # number of steps within the propagation duration

    dt, time_list = generate_timelist(init_time, duration, steps)
    prop_output = propagation_engine(line1, line2, time_list)
    max_begin_err, max_nominal_err = check_errors(
        init_time, dt, prop_output["sat"], prop_output["traj"], False
    )

    assert max_begin_err < 16 * u.mm
    assert max_nominal_err < 0.08 * u.mm


def test_interpolation_err_geo():
    """Tests the interpolation accuracy for GEO."""

    # Init TLE
    line1 = "1 99999U 12345A   20162.50918981  .00000000  00000-0  00000-0 0 00005"
    line2 = "2 99999 000.0000 124.6202 0000000 000.0000 000.0000 01.00273791000004"

    # test init time
    init_time = Time("2020:162:09:00:00", format="yday", scale="utc")

    duration = TimeDelta(3.0, format="jd")
    steps = 4000  # number of steps within the propagation duration

    dt, time_list = generate_timelist(init_time, duration, steps)
    prop_output = propagation_engine(line1, line2, time_list)
    max_err, max_nominal_err = check_errors(
        init_time, dt, prop_output["sat"], prop_output["traj"], False
    )

    assert max_err < 0.0005 * u.mm
    assert max_nominal_err < 0.0003 * u.mm


def test_interpolation_err_heo():
    """Tests the interpolation accuracy for Highly Elliptic Orbits."""

    # Init TLE
    line1 = "1 99999U 12345A   20162.50918981  .00000000  00000-0  00000-0 0 00005"
    line2 = "2 99999 000.0000 124.6202 0000000 000.0000 000.0000 01.00273791000004"

    # test init time
    init_time = Time("2020:162:00:00:00", format="yday", scale="utc")

    tle = TLE.from_tle(line1, line2)

    tle.epoch = init_time
    tle.mean_anomaly = 0
    tle.eccentricity = 0.4

    heo_line1, heo_line2 = tle.export_tle()

    # print(tle)

    duration = TimeDelta(3.0, format="jd")
    steps = 4000  # number of steps within the propagation duration

    dt, time_list = generate_timelist(init_time, duration, steps)
    prop_output = propagation_engine(heo_line1, heo_line2, time_list)
    max_err, max_nominal_err = check_errors(
        init_time, dt, prop_output["sat"], prop_output["traj"], False
    )

    assert max_err < 0.040 * u.mm
    assert max_nominal_err < 0.00030 * u.mm
