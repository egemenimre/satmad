# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Trajectory Class tests.
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, SkyCoord
from astropy.time import Time, TimeDelta
from sgp4.api import Satrec

from satmad.coordinates.frames import TEME
from satmad.coordinates.trajectory import Trajectory

_generate_plots = False


def test_interpolation_err():
    """Tests the interpolation accuracy."""

    # Init TLE
    line1 = "1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991"
    line2 = "2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482"

    # Init satellite object from the TLE
    sat = Satrec.twoline2rv(line1, line2)

    # ****** Generate the discrete time instants through the propagation duration ******
    init_time = Time(2458826, format="jd", scale="utc")
    duration = TimeDelta(3.0, format="jd")
    steps = 4000  # number of steps within the propagation duration
    dt = duration / steps

    # Generate time list
    dt_list = dt * np.arange(0, steps, 1)
    time_list = init_time + dt_list

    # SGP4 module requires time instances as jd and fraction arrays
    jd_list = time_list.jd1
    fr_list = time_list.jd2

    # ****** Generate the pos, vel vectors for each time instant ******

    # Run the propagation and init pos and vel vectors in TEME
    e, r_list, v_list = sat.sgp4_array(jd_list, fr_list)

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

    # **** Check the results ****
    step_offset = 0
    test_stepsize = 1.0 * u.s
    test_end = 0.01 * u.day

    print(f"Test duration: {(test_end - step_offset * dt).to(u.day)}")
    print(f"Step size: {test_stepsize}")

    time_list_hires = init_time + TimeDelta(
        np.arange(
            step_offset * dt.to_value(u.s),
            test_end.to_value(u.s),
            test_stepsize.to_value(u.s),
        ),
        format="sec",
    )
    jd_test_list = time_list_hires.jd1
    fr_test_list = time_list_hires.jd2

    # Generate the high-res test trajectory
    e_test_list, r_test_list, v_test_list = sat.sgp4_array(jd_test_list, fr_test_list)

    r_test_list = CartesianRepresentation(r_test_list.transpose(), unit=u.km)

    # Generate the error list
    r_err_list = (
        trajectory.get_coords(time_list_hires).cartesian.without_differentials()
        - r_test_list
    )

    # Error stats
    max_begin_err = r_err_list.norm().max()  # max expected error at begin
    max_nominal_err = r_err_list[-100:].norm().max()  # max expected error nominal

    assert max_begin_err < 1.4e-5 * u.km
    assert max_nominal_err < 7.4e-8 * u.km

    if _generate_plots:
        print(f"max error at the beginning  : {max_begin_err.to(u.mm)}")
        print(f"max error after settle down : {max_nominal_err.to(u.mm)}")

        # Plot the error per axis
        fr_ticker = np.arange(fr_test_list[0], fr_test_list[-1], dt.to_value(u.day))
        y_ticker = np.zeros(len(fr_ticker))
        plt.figure()
        plt.plot(
            fr_test_list,
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
