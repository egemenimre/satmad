# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Base module for analytical and numerical propagators.

"""
from abc import ABC, abstractmethod

import numpy as np
from astropy import units as u
from astropy.units import Quantity

from satmad.coordinates.trajectory import Trajectory
from satmad.core.celestial_bodies_lib import EARTH


class AbstractPropagator(ABC):
    """
    Base class for the propagators.

    Parameters
    ----------
    name : str
        Name of the propagator
    stepsize : Quantity or TimeDelta
        output stepsize for the propagator
    central_body : CelestialBody
        Centre Celestial Body for the propagator
    """

    def __init__(self, stepsize=60 * u.s, name="Propagator", central_body=EARTH):
        self.name = name
        self._stepsize = stepsize
        self.central_body = central_body

    @abstractmethod
    def gen_trajectory(self, init_orbit, interval, **kwargs):
        """
        Generates the trajectory (time and coordinate information) for the given
        interval with the internal output stepsize.

        Parameters
        ----------
        init_orbit
            Initial coordinates or orbit
        interval : TimeInterval
            Time interval for which the ephemeris will be generated
        kwargs
            Keywords specific to the implementation
        Returns
        -------
        Trajectory
            The output trajectory
        """
        pass

    @property
    def stepsize(self) -> Quantity:
        """Output stepsize for the propagator."""
        return self._stepsize

    @stepsize.setter
    def stepsize(self, stepsize):
        """
        Parameters
        ----------
        stepsize : TimeDelta or Quantity
        """
        self._stepsize = stepsize

    def _generate_time_list(self, interval):
        # generate number of steps (forced to rounded up int and added one to ensure
        # time list covers the entire interval)
        no_of_steps = np.ceil((interval.duration / self.stepsize).decompose()) + 1

        # make sure there are enough elements for interpolation
        if no_of_steps < Trajectory.reqd_min_elements():
            no_of_steps = Trajectory.reqd_min_elements()

        # output time list
        time_list = interval.start + self.stepsize * np.arange(0, no_of_steps)
        time_list.format = "isot"

        return time_list
