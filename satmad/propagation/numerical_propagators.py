# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Base module for numerical propagators.

"""
from enum import Enum

import numpy as np
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord,
)
from scipy.integrate import solve_ivp

from satmad.coordinates.trajectory import Trajectory
from satmad.core.celestial_bodies import EARTH
from satmad.propagation.base_propagator import AbstractPropagator
from satmad.propagation.force_models import two_body_accel


class ODESolverType(Enum):
    """
    Built-in solver types for :meth:`~scipy.integrate.solve_ivp`.
    """

    RK45 = "RK45"
    RK23 = "RK23"
    DOP853 = "DOP853"
    RADAU = "Radau"
    BDF = "BDF"
    LSODA = "LSODA"


class NumericalPropagator(AbstractPropagator):
    """
    Numerical orbit propagator based on Scipy ODE Solvers.

    Analysis shows that `DOP853` is currently the best (and default) choice for
    accuracy and runtime performance.

    Parameters
    ----------
    name : str
        Name or identifier for the propagator
    stepsize : Quantity or `~astropy.time.TimeDelta`
        Output stepsize for the propagator
    solver_type : ODESolverType
        Solver type to be used
    rtol : float
        relative tolerance value
    atol : float
        absolute tolerance value
    central_body : CelestialBody
        Centre Celestial Body for the propagator
    """

    def __init__(
        self,
        stepsize=60 * u.s,
        solver_type=ODESolverType.DOP853,
        rtol=1e-12,
        atol=1e-14,
        name="",
        central_body=EARTH,
    ):
        super().__init__(stepsize, name, central_body)

        self._solver_type = solver_type

        self.rtol = rtol
        self.atol = atol

        if not name:
            self._name = f"{solver_type.value} (rtol: {rtol}, atol: {atol} )"

    def gen_trajectory(self, init_coords, interval):
        """
        Generates the trajectory (time and coordinate information) for the given
        interval with the internal output stepsize.

        If the initial coordinate is in a different frame than the propagation frame,
        it is automatically converted to the proper frame for propagation.

        Parameters
        ----------
        init_coords : SkyCoord
            Initial coordinates (the first value is used)
        interval : TimeInterval
            Time interval for which the ephemeris will be generated
        Returns
        -------
        Trajectory
            The output trajectory
        """
        # generate the output timelist
        time_list = self._generate_time_list(interval)

        # Propagate the orbit
        if init_coords.isscalar:
            # Single SkyCoord element available
            init_coords_converted = init_coords
        else:
            # Multiple SkyCoord elements available, use first
            init_coords_converted = init_coords[0]

        # Check for the initial coordinate system, convert if necessary
        if self.central_body.inert_coord != init_coords_converted.frame.name:
            init_coords_converted = init_coords_converted.transform_to(
                self.central_body.inert_coord
            )

        # Run propagation
        coords_list = self._run_propagation(init_coords_converted, time_list)

        # generate the trajectory
        return Trajectory(coords_list)

    def _run_propagation(self, init_coords, time_list):
        """Solves the ODE or runs the actual propagation.

        Parameters
        ----------
        init_coords : SkyCoord
            Initial coordinates (the first value is used)
        time_list : Time
            Time list for outputs
        Returns
        -------
        SkyCoord
            The output coordinates corresponding to `time_list`
        """

        # generate the time list since epoch (in seconds)
        t_list_epoch = (time_list - time_list[0]).to_value(u.s)

        # concat the init vector
        rv_init = [
            *init_coords.cartesian.without_differentials().xyz.to_value(u.km),
            *init_coords.velocity.d_xyz.to_value(u.km / u.s),
        ]

        solver = solve_ivp(
            _ode_diff_eqns,
            (t_list_epoch[0], t_list_epoch[-1]),
            rv_init,
            method=self._solver_type.value,
            dense_output=True,
            rtol=self.rtol,
            atol=self.atol,
        )

        r_list = np.empty((3, len(t_list_epoch)))
        v_list = np.empty((3, len(t_list_epoch)))
        # Fill the raw output data
        i = 0
        for t in t_list_epoch:
            rv = solver.sol(t)
            r_list[:, i] = rv[:3]
            v_list[:, i] = rv[3:]
            i = i + 1

        # Load the time, pos, vel info into astropy objects (shallow copied)
        rv_list = CartesianRepresentation(r_list, unit=u.km).with_differentials(
            CartesianDifferential(v_list, unit=u.km / u.s)
        )

        coords_list = SkyCoord(
            rv_list,
            obstime=time_list,
            frame=GCRS,
            representation_type="cartesian",
            differential_type="cartesian",
            copy=False,
        )
        return coords_list


__km3s2 = u.km ** 3 / u.s ** 2


def _ode_diff_eqns(t, rv, mu=EARTH.mu):
    """
    Defines the Ordinary Differential Equations that govern the motion of the
    satellite. In practice, this method computes the instantaneous accelerations
    acting on the satellite.

    Note that the input `rv` store the position and velocity vectors in the
    inertial frame of the central body. For the Earth this is GCRS.

    This method is used by the scipy `solve_ivp()` solver to evaluate the differential
    equations at a given time.

    Parameters
    ----------
    t : Time
        Time where the accelerations are defined
    rv : ndarray
        array of 6 values, the first three is the position vector (km) and
        the second three is the velocity vector (km/s)
    mu : Quantity
        GM value in :math:`km^3 / s^2`

    Returns
    -------
    ndarray
        The array with velocity and acceleration information.
    """

    r = rv[:3]
    v = rv[3:]

    # two-body acceleration
    a = two_body_accel(r, mu.to_value(__km3s2))

    # concat the two lists using PEP-448 unpack
    return [*v, *a]
