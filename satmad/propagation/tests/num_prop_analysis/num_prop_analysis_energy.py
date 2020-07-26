# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Analysis of numerical propagation errors with various numerical propagation methods.
Case: Energy error and runtimes.
"""
from typing import Dict

from astropy import units as u

from satmad.propagation.numerical_propagators import ODESolverType
from satmad.propagation.tests.num_prop_analysis.num_prop_analysis_engine import (
    analyse_energy_along_trajectory,
    init_rv,
    plot_energy_diff,
    plot_runtimes,
    tle_heo,
)

if __name__ == "__main__":
    # Set up propagation config
    # tle = tle_leo()
    tle = tle_heo()
    rv_init = init_rv(tle, tle.epoch)
    stepsize = 120 * u.s
    solver_types = [solver_type for solver_type in ODESolverType]
    # solver_types.remove(ODESolverType.BDF)
    # solver_types.remove(ODESolverType.RK23)
    # solver_types.clear()
    # solver_types.append(ODESolverType.DOP853)
    rtol = 3e-14
    # atol = [1e-12, 1e-12, 1e-12, 1e-15, 1e-15, 1e-15]
    atol = 1e-12
    init_time_offset = 0 * u.day
    duration = 1.0 * u.day

    cases: Dict[str, Dict] = {}
    type(cases)

    # run all cases
    for solver_type in solver_types:
        name = solver_type.value

        analyse_energy_along_trajectory(
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

    session_name = f"[rtol: {rtol}, atol: {atol}]"

    # Plot energy differences
    plot_energy_diff(cases, title=session_name, logscale_plot=True)

    # Plot DOP853 on linear scale to see its behaviour
    plot_energy_diff(
        {"DOP853": cases["DOP853"]}, title=session_name, logscale_plot=False
    )

    # Plot runtimes
    plot_runtimes(
        cases, title=f"{session_name} [duration: {duration.to_value(u.day)} days]"
    )
