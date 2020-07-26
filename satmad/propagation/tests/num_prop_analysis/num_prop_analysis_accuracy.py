# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Analysis of numerical propagation errors with DOP853 numerical propagation method.
Case: Position accuracy
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
    tle_leo,
)

if __name__ == "__main__":
    # Set up propagation config
    tle = tle_leo()
    rv_init = init_rv(tle, tle.epoch)
    stepsize = 120 * u.s
    solver_type = ODESolverType.DOP853
    rtol = 3e-14
    # atols = [1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10]
    rtols = [1e-14, 1e-13, 1e-12, 1e-11, 1e-10]
    atol = 1e-14
    init_time_offset = 0 * u.day
    duration = 10.0 * u.day

    # atol effect for fixed low rtol: no impact on pos error
    # atol effect for fixed high rtol: no impact on pos error, lower accuracy
    # rtol effect for fixed low atol: pos error varies with rtol
    # rtol effect for fixed high atol: pos error varies with rtol, same accuracy

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

    cases: Dict[str, Dict] = {}
    type(cases)

    # run all cases
    for rtol in rtols:
        name = str(rtol)
        # rtol = atol

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

    session_name = f"[DOP853, atol = {atol} , varying rtol]"

    # Plot pos diff differences
    plot_pos_diff(cases, title=session_name, logscale_plot=True)

    # Plot runtime for one case to see its behaviour
    plot_pos_diff(
        {"1e-10": cases["1e-10"]}, title=f"{session_name}", logscale_plot=False
    )

    # Plot runtimes
    plot_runtimes(
        cases,
        title=f"{session_name} [duration: {duration.to_value(u.day)} days]",
        logscale_plot=True,
    )
