# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
SGP4 Propagator.

"""
import numpy as np
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord,
)

from satmad.coordinates.frames import TEME
from satmad.propagation.base_propagator import AbstractPropagator


class SGP4Propagator(AbstractPropagator):
    """
    SGP4 Analytical propagator using `TEME` mean orbital elements.

    Parameters
    ----------
    name : str
        Name of the propagator
    stepsize : Quantity or TimeDelta
        output stepsize for the propagator
    """

    def __init__(self, name="SGP4", stepsize=60 * u.s):
        super().__init__(name, stepsize)

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
        SGP4Propagator.__handle_err_code(e, tle)

        v_teme = CartesianDifferential(np.asarray(v), unit=u.km / u.s)
        r_teme = CartesianRepresentation(np.asarray(r), unit=u.km)
        rv_gcrs = TEME(r_teme.with_differentials(v_teme), obstime=time).transform_to(
            GCRS(obstime=time)
        )

        return SkyCoord(
            rv_gcrs, representation_type="cartesian", differential_type="cartesian",
        )

    @staticmethod
    def __handle_err_code(e, tle):
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

    # def generate_eph(self, tle, interval):
    #     pass


class SGP4GeneralError(Exception):
    """General SGP4 Propagation errors, usually orbital elements out of range."""

    pass


class SGP4SatDecayedError(Exception):
    """Satellite coordinates at required coordinates are below decay altitude."""

    pass
