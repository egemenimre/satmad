# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Base module for analytical and numerical propagators.

"""
from abc import ABC

from astropy.units import Quantity
import numpy as np
from astropy import units as u


from satmad.coordinates.trajectory import Trajectory


class AbstractPropagator(ABC):
    """
    Base class for the propagators.

    Parameters
    ----------
    name : str
        Name of the propagator
    stepsize : Quantity or TimeDelta
        output stepsize for the propagator
    """

    def __init__(self, stepsize=60 * u.s, name="Propagator"):
        self.name = name
        self._stepsize = stepsize

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
