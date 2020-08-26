# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Methods to analyse and measure numerical propagation performance.

"""
from time import perf_counter

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import quantity_support, time_support
from matplotlib import pyplot as plt

from satmad.coordinates.trajectory import Trajectory
from satmad.core.celestial_bodies import EARTH, GM_earth
from satmad.propagation.force_models import two_body_energy
from satmad.propagation.numerical_propagators import NumericalPropagator
from satmad.propagation.sgp4_propagator import SGP4Propagator
from satmad.propagation.tle import TLE
from satmad.utils.timeinterval import TimeInterval


def tle_heo():
    """Generates the TLE with HEO test setup (e=0.70)."""

    name = "ATLAS 2A CENTAUR R/B"
    line1 = "1 23840U 96020B   18198.38669861 -.00000081  00000-0  00000+0 0  9997"
    line2 = "2 23840  21.5075 284.9295 6975939  94.1963 356.6140  2.24211004182606"

    return TLE.from_tle(line1, line2, name)


def tle_leo():
    """Generates the TLE with near-circular LEO test setup."""

    name = "RASAT"
    line1 = "1 37791U 11044D   18198.20691930 -.00000011  00000-0  70120-5 0  9992"
    line2 = "2 37791  98.1275 290.4108 0021116 321.0704  38.8990 14.64672859369594"

    return TLE.from_tle(line1, line2, name)


def tle_geo():
    """Generates the TLE with near-circular LEO test setup."""

    name = "TURKSAT 4B"
    line1 = "1 40984U 15060A   18198.04228921  .00000095  00000-0  00000+0 0  9995"
    line2 = "2 40984   0.0142 325.6420 0001708 151.3463 243.0113  1.00269765 10118"

    return TLE.from_tle(line1, line2, name)


def init_rv(tle, epoch):
    """Initialises the position and velocity vectors at the requested epoch."""
    # init SGP4
    return SGP4Propagator().propagate(tle, epoch)


def propagation_engine(
    rv_init, stepsize, solver_type, init_time_offset, duration, rtol, atol
):
    """Initialises and runs the propagation engine to yield the resulting trajectory."""

    # init propagator
    prop = NumericalPropagator(
        stepsize,
        solver_type=solver_type,
        rtol=rtol,
        atol=atol,
        name="",
        central_body=EARTH,
    )

    prop_start = rv_init.obstime + init_time_offset
    prop_duration = duration

    trajectory = prop.gen_trajectory(rv_init, TimeInterval(prop_start, prop_duration))

    return trajectory


def energy_along_trajectory(
    coord_list: SkyCoord, mu=GM_earth.to_value(u.km ** 3 / u.s ** 2)
):
    """Computes the absolute and relative (to initial) energy along the trajectory."""

    time_list = coord_list.obstime - coord_list[0].obstime

    init_coords = coord_list[0]
    init_energy = two_body_energy(
        init_coords.cartesian.without_differentials().xyz.to_value(u.km),
        init_coords.velocity.d_xyz.to_value(u.km / u.s),
        mu,
    )

    energy_list = []
    energy_diff_list = []
    for coords in coord_list:
        r = coords.cartesian.without_differentials().xyz.to_value(u.km)
        v = coords.velocity.d_xyz.to_value(u.km / u.s)
        energy = two_body_energy(r, v, mu)
        energy_list.append(energy)
        energy_diff_list.append(energy - init_energy)

    return time_list, energy_list, energy_diff_list


def plot_energy_diff(cases, title, logscale_plot=True):
    """Plots the energy difference."""
    quantity_support()
    time_support()

    fig, ax = plt.subplots()

    ax.grid()
    # ax.xaxis.set_major_formatter(FormatStrFormatter("% 3.3f"))
    # ax.yaxis.set_major_formatter(FormatStrFormatter("% 3.3f"))
    plt.title(f"Specific Energy Evolution wrt. Initial State\n({title})")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Specific Energy Diff (m^2/s^2)")

    if logscale_plot:
        ax.set_yscale("log")

    fig.set_size_inches(8, 6)

    for case in cases:
        if logscale_plot:
            # Make sure everything is positive for plotting
            ax.plot(
                cases[case]["time_list"].to(u.day),
                np.abs(cases[case]["energy_diff_list"]),
                label=case,
            )
        else:
            ax.plot(
                cases[case]["time_list"].to(u.day),
                cases[case]["energy_diff_list"],
                label=case,
            )

    plt.legend()
    plt.show()


def plot_runtimes(cases, title, logscale_plot=True):
    """Plots the runtimes."""
    quantity_support()
    time_support()

    fig, ax = plt.subplots()

    # ax.grid()
    # ax.xaxis.set_major_formatter(FormatStrFormatter("% 3.3f"))
    # ax.yaxis.set_major_formatter(FormatStrFormatter("% 3.3f"))
    plt.title(f"Runtime Analysis\n({title})")
    ax.set_ylabel("Runtime (s)")

    if logscale_plot:
        ax.set_yscale("log")

    fig.set_size_inches(8, 6)

    # generate data lists
    names = [case for case in cases]
    runtimes = [cases[case]["runtime"] for case in cases]

    # generate the plots
    if logscale_plot:
        # Make sure everything is positive for plotting
        ax.bar(names, [np.abs(runtime) for runtime in runtimes])
    else:
        ax.bar(names, runtimes)

    # Annotate the plot
    for x, y in zip(names, runtimes):
        label = "{:.2f}".format(y)
        plt.annotate(
            label,  # this is the text
            (x, y),  # this is the point to label
            textcoords="offset points",  # how to position the text
            xytext=(0, 5),  # distance from text to points (x,y)
            ha="center",
        )  # horizontal alignment can be left, right or center

    # ax.legend([case for case in cases])
    plt.show()


def analyse_energy_along_trajectory(
    cases, name, rv_init, stepsize, solver_type, init_time_offset, duration, rtol, atol
):
    """Analyses the energy evolution along the trajectory."""

    print(f"**** Case: {name} ****")

    tic = perf_counter()

    # run propagation and get trajectory
    trajectory = propagation_engine(
        rv_init, stepsize, solver_type, init_time_offset, duration, rtol, atol
    )

    toc = perf_counter()

    runtime = (toc - tic) * u.s

    # Compute the energy along the trajectory
    time_list, energy_list, energy_diff_list = energy_along_trajectory(
        trajectory.coord_list
    )
    cases[name] = {
        "traj": trajectory,
        "time_list": time_list,
        "energy_list": energy_list,
        "energy_diff_list": energy_diff_list,
        "runtime": runtime,
    }

    print(f"Runtime: {runtime:0.4f} seconds")

    # Print summary results

    abs_energy_diff_list = abs(np.array(energy_diff_list))

    # Max energy diff (abs value) - should be low
    energy_diff_max = abs_energy_diff_list.max()
    # Mean energy diff - should be zero to have zero mean energy
    energy_diff_mean = np.array(energy_diff_list).mean()

    # initial energy diff and final energy diff to assess energy leak
    # First and last 100 elements taken and they should also ideally be zero
    init_energy_diff_mean = np.array(energy_diff_list)[:100].mean()
    final_energy_diff_mean = np.array(energy_diff_list)[-100:].mean()

    print(f"max energy diff: {energy_diff_max}")
    print(f"mean energy    : {energy_diff_mean}")

    print(f"init mean energy diff : {init_energy_diff_mean}")
    print(f"final mean energy diff: {final_energy_diff_mean}")


def pos_along_trajectory(coord_list: SkyCoord, truth_traj: Trajectory):
    """Computes the absolute and relative (to "truth") position along the trajectory."""

    time_list = coord_list.obstime - coord_list[0].obstime

    pos_diff_list = []
    for coords in coord_list:
        r = coords.cartesian.without_differentials()
        # v = coords.velocity.d_xyz.to_value(u.km / u.s)

        r_truth = truth_traj(coords.obstime).cartesian.without_differentials()

        pos_diff_list.append((r - r_truth).norm().to_value(u.m))

    return time_list, pos_diff_list


def analyse_pos_along_trajectory(
    truth_traj,
    cases,
    name,
    rv_init,
    stepsize,
    solver_type,
    init_time_offset,
    duration,
    rtol,
    atol,
):
    """Analyses the position variation along the trajectory."""

    print(f"**** Case: {name} ****")

    tic = perf_counter()

    # run propagation and get trajectory
    trajectory = propagation_engine(
        rv_init, stepsize, solver_type, init_time_offset, duration, rtol, atol
    )

    toc = perf_counter()

    runtime = (toc - tic) * u.s

    # Compute the energy along the trajectory
    time_list, pos_diff_list = pos_along_trajectory(trajectory.coord_list, truth_traj)
    cases[name] = {
        "traj": trajectory,
        "time_list": time_list,
        "pos_diff_list": pos_diff_list,
        "runtime": runtime,
    }

    print(f"Runtime: {runtime:0.4f} seconds")

    # Print summary results

    abs_pos_diff_list = abs(np.array(pos_diff_list))

    # Max pos diff (abs value) - should be low
    pos_diff_max = abs_pos_diff_list.max()
    # # Mean pos diff - should be zero to have zero mean energy
    # pos_diff_mean = np.array(pos_diff_list).mean()
    #
    # # initial pos diff and final pos diff to assess energy leak
    # # First and last 100 elements taken and they should also ideally be zero
    # init_pos_diff_mean = np.array(pos_diff_list)[:100].mean()
    # final_pos_diff_mean = np.array(pos_diff_list)[-100:].mean()

    print(f"max pos diff: {pos_diff_max} m")
    # print(f"mean pos diff   : {pos_diff_mean}")
    #
    # print(f"init mean energy diff : {init_pos_diff_mean}")
    # print(f"final mean energy diff: {final_pos_diff_mean}")


def plot_pos_diff(cases, title, logscale_plot=True):
    """Plots the pos difference."""
    quantity_support()
    time_support()

    fig, ax = plt.subplots()

    ax.grid()
    # ax.xaxis.set_major_formatter(FormatStrFormatter("% 3.3f"))
    # ax.yaxis.set_major_formatter(FormatStrFormatter("% 3.3f"))
    plt.title(f"Position Difference Variation\n({title})")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Pos Diff (m)")

    if logscale_plot:
        ax.set_yscale("log")

    fig.set_size_inches(8, 6)

    for case in cases:
        if logscale_plot:
            # Make sure everything is positive for plotting
            ax.plot(
                cases[case]["time_list"].to(u.day),
                np.abs(cases[case]["pos_diff_list"]),
                label=case,
            )
        else:
            ax.plot(
                cases[case]["time_list"].to(u.day),
                cases[case]["pos_diff_list"],
                label=case,
            )

    ax.legend()
    plt.show()
