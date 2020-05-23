# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Time interval module.

`TimeInterval` class stores time intervals. Module also offers additional functionality
such as operations between `TimeInterval` objects (union, intersection etc.).

"""

from astropy import units as u
from astropy.time import Time, TimeDelta


class TimeInterval:
    """
    Represent and manipulate time intervals.

    Uses Astropy `Time` classes under the hood.

    Parameters
    ----------
    init_times : Time
    end_times : Time or TimeDelta
    replicate : bool
        `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False` to
        use a shallow copy to save memory

    """

    # TODO
    #     flag end inclusive
    #     flag begin inclusive

    _EPS_TIME = 10 * u.us
    """Allowable time threshold, this much 'out of bounds' is allowed when handling
    requested interpolation times. This helps with floating point artifacts such as
    round-off errors."""

    def __init__(self, init_times, end_times, replicate=False):

        if replicate:
            # replicate internal time data
            self._init_times = init_times.replicate()
        else:
            # shallow copy the internal time data
            self._init_times = init_times

        if isinstance(end_times, TimeDelta):
            # End times are of type `TimeDelta`, make sure array sizes match
            # - or that there is a single duration
            try:
                self._end_times = self._init_times + end_times
            except Exception as err:
                raise ValueError(
                    f"Array size mismatch between "
                    f"Initial Times (size: {len(self._init_times)}) and "
                    f"Durations (size: {len(end_times)})."
                ) from err

        elif isinstance(end_times, Time):
            # End times are of type `Time`, make sure array sizes match
            if len(end_times) == len(self._init_times):
                if replicate:
                    # replicate internal time data
                    self._end_times = end_times.replicate()
                else:
                    # shallow copy the internal time data
                    self._end_times = end_times
            else:
                raise ValueError(
                    f"Array size mismatch between "
                    f"Initial Times (size: {len(self._init_times)}) and "
                    f"End Times (size: {len(end_times)})."
                )

    def __repr__(self):
        return (
            f"<{self.__class__.__name__} object "
            f"(stores {self._init_times.__class__.__name__} objects) : "
            f"scale='{self._init_times.scale}' "
            f"format='{self._init_times.format}' "
            f"init times={getattr(self._init_times, self._init_times.format)}>"
            f"end times={getattr(self._end_times, self._end_times.format)}>"
        )

    def get_interval(self, index):
        """
        Gets the time interval for the given index.

        Parameters
        ----------
        index : int
            requested index

        Returns
        -------
        init_time, end_time
            the (initial time, end time) tuple corresponding to the index

        Raises
        ------
        IndexError
            Requested index is out of bounds

        """
        return self._init_times[index], self._end_times[index]

    def __str__(self):
        txt = ""
        for i, init_time in enumerate(self._init_times):
            txt += f"[ {str(getattr(self._init_times[i], self._init_times.format))}"
            txt += f"  {str(getattr(self._end_times[i], self._end_times.format))} ]\n"

        return txt
