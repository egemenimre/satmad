# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Analysis of numerical propagation errors with DOP853 numerical propagation method.
Case: Eccentricity and pos accuracy.
"""
from typing import Dict

from astropy import units as u

from satmad.propagation.numerical_propagators import ODESolverType
from satmad.propagation.tests.num_prop_analysis.num_prop_analysis_engine import (
    analyse_pos_along_trajectory,
    init_rv,
    plot_pos_diff,
    plot_runtimes,
    propagation_engine,
    tle_heo,
)

if __name__ == "__main__":

    # Set up propagation config
    stepsize = 120 * u.s
    solver_type = ODESolverType.DOP853
    rtol = 1e-11
    atol = 1e-11
    init_time_offset = 0 * u.day
    duration = 10.0 * u.day

    ecc_list = [0.001, 0.1, 0.2, 0.4, 0.5, 0.7, 0.75]
    tle_list = []
    for ecc in ecc_list:
        tle = tle_heo()
        tle.eccentricity = ecc
        tle_list.append(tle)

    cases: Dict[str, Dict] = {}
    type(cases)

    # run all cases
    for tle in tle_list:
        name = str(tle.eccentricity)

        rv_init = init_rv(tle, tle.epoch)

        # run propagation and get truth trajectory
        truth_traj = propagation_engine(
            rv_init,
            stepsize,
            ODESolverType.DOP853,
            init_time_offset,
            duration,
            rtol=3e-14,
            atol=1e-15,
        )

        analyse_pos_along_trajectory(
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
        )

    session_name = f"[DOP853, atol = {atol} , rtol = {rtol}]"

    # Plot pos differences
    plot_pos_diff(cases, title=session_name, logscale_plot=True)

    # Plot runtimes
    plot_runtimes(
        cases,
        title=f"{session_name} [duration: {duration.to_value(u.day)} days]",
        logscale_plot=False,
    )
