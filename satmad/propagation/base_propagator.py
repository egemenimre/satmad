# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Base classes for analytical and numerical propagators.

"""
from abc import ABC

from astropy import units as u
from astropy.units import Quantity


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

    def __init__(self, name, stepsize=60 * u.s):
        self._name = name
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

    @property
    def name(self) -> str:
        """Name of the propagator."""
        return self._name
