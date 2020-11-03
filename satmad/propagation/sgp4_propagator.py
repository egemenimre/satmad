# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
SGP4 Propagator to propagate Earth bound satellites.

"""
import numpy as np
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    TEME,
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord,
)

from satmad.coordinates.trajectory import Trajectory
from satmad.core.celestial_bodies import EARTH
from satmad.propagation.base_propagator import AbstractPropagator


class SGP4Propagator(AbstractPropagator):
    """
    SGP4 Analytical propagator using `TEME` mean orbital elements.

    While `central_body` of the propagator is `EARTH`, its constants are not used.

    Parameters
    ----------
    name : str
        Name of the propagator
    stepsize : Quantity or TimeDelta
        output stepsize for the propagator
    """

    def __init__(self, stepsize=60 * u.s, name="SGP4"):
        super().__init__(stepsize, name, central_body=EARTH)

    @staticmethod
    def propagate(tle, time):
        """
        Computes the coordinates at the target `time` using the source `TLE` mean
        orbital elements.

        Parameters
        ----------
        tle : TLE
            Mean orbital elements (Two-Line-Element)
        time : Time
            Target time for propagation

        Returns
        -------
        coords : SkyCoord
            Coordinates at target time (in `GCRS`)

        Raises
        ------
        SGP4GeneralError
            Errors in SGP4 propagation
        SGP4SatDecayedError
            Satellite coordinates at required coordinates are below decay altitude.
        """
        # compute coords at time
        e, r, v = tle.satrec.sgp4(time.utc.jd1, time.utc.jd2)

        # check error code
        SGP4Propagator.__handle_err_code(e)

        v_teme = CartesianDifferential(np.asarray(v), unit=u.km / u.s)
        r_teme = CartesianRepresentation(np.asarray(r), unit=u.km)
        rv_gcrs = TEME(r_teme.with_differentials(v_teme), obstime=time).transform_to(
            GCRS(obstime=time)
        )

        return SkyCoord(
            rv_gcrs, representation_type="cartesian", differential_type="cartesian"
        )

    def gen_trajectory(self, tle, interval):
        """
        Generates the trajectory (time and coordinate information) for the given
        interval with the internal stepsize.
        Parameters
        ----------
        tle : TLE
            Two-Line-Element initial orbit information (TEME mean orbital elements)
        interval : TimeInterval
            Time interval for which the ephemeris will be generated
        Returns
        -------
        Trajectory
            The output trajectory (in `GCRS`)
        """

        # ****** Generate the pos, vel vectors for each time instant ******

        # generate the output timelist
        time_list = self._generate_time_list(interval)

        # Run the propagation and init pos and vel vectors in TEME
        e, r_list, v_list = tle.satrec.sgp4_array(time_list.jd1, time_list.jd2)

        # Load the time, pos, vel info into astropy objects (shallow copied)
        rv_list_teme = CartesianRepresentation(
            r_list, unit=u.km, xyz_axis=1
        ).with_differentials(CartesianDifferential(v_list, unit=u.km / u.s, xyz_axis=1))

        rv_list_gcrs = TEME(
            rv_list_teme,
            obstime=time_list,
            representation_type="cartesian",
            differential_type="cartesian",
        ).transform_to(GCRS(obstime=time_list))

        # trajectory in astropy
        traj_astropy = SkyCoord(
            rv_list_gcrs,
            obstime=time_list,
            frame="gcrs",
            representation_type="cartesian",
            differential_type="cartesian",
        )

        # Init trajectory in Trajectory object
        trajectory = Trajectory(traj_astropy)

        return trajectory

    @staticmethod
    def __handle_err_code(e):
        """
        Handle the error code out of the propagation.

        Error codes are:
            0. no error
            1. mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 earth radius
            2. mean motion less than 0.0
            3. perturbed eccentricity out of limits, ecc < 0.0  or  ecc > 1.0
            4. semi-latus rectum < 0.0
            5. epoch elements are sub-orbital
            6. satellite has decayed

        Parameters
        ----------
        e : int
            error code

        Raises
        ------
        SGP4GeneralError
            Errors in SGP4 propagation - satellite orbit evolved to an invalid state.
        SGP4SatDecayedError
            Satellite coordinates at required coordinates are below decay altitude.
        """
        if e == 0:
            pass
        if e == 1:
            raise SGP4GeneralError(
                "SGP4 Propagation Error: "
                "ecc >= 1.0 or ecc < -0.001 or a < 0.95 earth radius"
            )
        if e == 2:
            raise SGP4GeneralError("SGP4 Propagation Error: mean motion less than 0.0")
        if e == 3:
            raise SGP4GeneralError(
                "SGP4 Propagation Error: "
                "perturbed eccentricity out of limits, ecc < 0.0  or  ecc > 1.0"
            )
        if e == 4:
            raise SGP4GeneralError("SGP4 Propagation Error: semi-latus rectum < 0.0")
        if e == 5:
            raise SGP4GeneralError(
                "SGP4 Propagation Error: epoch elements are sub-orbital"
            )
        if e == 6:
            raise SGP4SatDecayedError(
                "SGP4 Propagation Error:"
                "Satellite has decayed (Sat number: {tle.sat_number})"
            )


class SGP4GeneralError(Exception):
    """General SGP4 Propagation errors, usually orbital elements out of range."""

    pass


class SGP4SatDecayedError(Exception):
    """Satellite coordinates at required coordinates are below decay altitude."""

    pass
