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

_EPS_TIME = 10 * u.us
"""Allowable time threshold, this much 'out of bounds' is allowed when assuming two
instances of time are *practically* equal. This helps with floating point artifacts such as
round-off errors."""


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

    # TODO check for out of validity errors
    #  when doing Union with an external interval
    #
    # TODO  take union etc methods to within class
    #  with thin wrappers and checks and multiple items

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


def _is_within_interval(time, ref_interval, inclusive=True):
    """
    Checks whether the requested time is within the reference interval.

    - If the `inclusive` flag is True: ref_interval_begin <= time <= ref_interval_end
    - If the `inclusive` flag is False: ref_interval_begin < time < ref_interval_end

    Parameters
    ----------
    time : Time
        time to be checked
    ref_interval : (Time, Time)
        Reference interval (Tuple of initial time and end time)
    inclusive : bool
        True if max and min boundaries are included in the check, False otherwise

    Returns
    -------
    bool
        True if time is within the reference interval, False otherwise
    """
    if inclusive:
        if ref_interval[0].value <= time.value <= ref_interval[1].value:
            return True
        else:
            return False
    else:
        if ref_interval[0].value < time.value < ref_interval[1].value:
            return True
        else:
            return False


def _is_intersecting(interval, ref_interval):
    """
    Checks whether the requested interval intersects (or is contained within)
    the reference interval.

    Parameters
    ----------
    interval : (Time, Time)
        Check interval (Tuple of initial time and end time)
    ref_interval : (Time, Time)
        Reference interval (Tuple of initial time and end time)

    Returns
    -------
    bool
        True if check interval intersects with the reference interval, False otherwise
    """
    if _contains(interval, ref_interval):
        # interval is contained or the same as ref interval
        return True
    elif _is_within_interval(interval[0], ref_interval, False) or _is_within_interval(
        interval[1], ref_interval, False
    ):
        # at least interval_begin or interval_end is within the interval
        return True
    else:
        return False


def _contains(interval, ref_interval):
    """
    Checks whether the requested interval is contained within the reference interval.

    Parameters
    ----------
    interval : (Time, Time)
        Check interval (Tuple of initial time and end time)
    ref_interval : (Time, Time)
        Reference interval (Tuple of initial time and end time)

    Returns
    -------
    bool
        True if check interval is contained within the reference interval,
        False otherwise
    """
    if _is_within_interval(interval[0], ref_interval, True) and _is_within_interval(
        interval[1], ref_interval, True
    ):
        return True
    else:
        return False


def _intersect(interval, ref_interval):
    """
    Intersection operator for two time intervals.

    Returns a new interval that is
    the Intersection of two intervals, or None if there is no intersection.

    Parameters
    ----------
    interval : (Time, Time)
        Check interval (Tuple of initial time and end time)
    ref_interval : (Time, Time)
        Reference interval (Tuple of initial time and end time)

    Returns
    -------
    (Time, Time)
        A new interval (Tuple of initial time and end time) that is
        the Intersection of two intervals, or `None` if there is no intersection

    """
    if _is_intersecting(interval, ref_interval):

        if interval[0] < ref_interval[0]:
            # interval begins before ref_interval - set ref_interval begin
            new_init = ref_interval[0]
        else:
            # interval begins after ref_interval - set interval begin
            new_init = interval[0]

        if ref_interval[1] < interval[1]:
            # interval ends after ref_interval - set ref_interval end
            new_end = ref_interval[1]
        else:
            # interval ends before ref_interval - set interval end
            new_end = interval[1]

        return new_init, new_end
    else:
        # No intersection found
        return None


def _union(interval, ref_interval):
    """
    Union operator for two time intervals.

    Returns a new interval that is
    the Union of two intervals, or None if there is no intersection.

    Parameters
    ----------
    interval : (Time, Time)
        Check interval (Tuple of initial time and end time)
    ref_interval : (Time, Time)
        Reference interval (Tuple of initial time and end time)

    Returns
    -------
    (Time, Time)
        A new interval (Tuple of initial time and end time) that is
        the Union of two intervals, or `None` if there is no intersection

    """
    if _is_intersecting(interval, ref_interval):

        if interval[0] < ref_interval[0]:
            # interval begins before ref_interval - set interval begin
            new_init = interval[0]
        else:
            # interval begins after ref_interval - set ref_interval begin
            new_init = ref_interval[0]

        if ref_interval[1] < interval[1]:
            # interval ends after ref_interval - set interval end
            new_end = interval[1]
        else:
            # interval ends before ref_interval - set ref_interval end
            new_end = ref_interval[1]

        return new_init, new_end
    else:
        # No intersection found
        return None


def _duration(interval):
    """
    Computes the duration of an interval.

    Parameters
    ----------
    interval : (Time, Time)
        Interval (Tuple of initial time and end time)

    Returns
    -------
    duration : TimeDelta
        Duration of the interval (can be negative)
    """
    return interval[1] - interval[0]


def _expand(interval, delta_begin=0, delta_end=0):
    """
    Expands (or shrinks) the given interval.

    Generates a new tuple, where:
    - new interval begin: interval[0] - delta_begin
    - new interval end: interval[1] + delta_end

    Negative begin and/or end times are possible (to shrink the interval),
    though values ending in a negative interval will throw a `ValueError`.

    Parameters
    ----------
    interval : (Time, Time)
        Interval (Tuple of initial time and end time)
    delta_begin : TimeDelta
        The delta time to expand the beginning of the interval
        (or to shrink, with negative values)
    delta_end : TimeDelta
        The delta time to expand the end of the interval
        (or to shrink, with negative values)
    Returns
    -------
    (Time, Time)
        A new interval (Tuple of initial time and end time) that is
        the result of the requested change

    Raises
    ------
    ValueError
        Raised if the requested expansion results in a negative duration interval
    """
    new_interval = (interval[0] - delta_begin, interval[1] + delta_end)

    if _duration(new_interval) > _EPS_TIME:
        return new_interval
    else:
        raise ValueError(
            "Duration of the expanded/shrunk interval is negative or zero."
        )


def _interval_equals(interval, ref_interval):
    """
    Checks whether two intervals are (almost) equal in value.

    If the two begin or end values are as close as `_EPS_TIME`, then they are assumed
    equal in value.

    Parameters
    ----------
     interval : (Time, Time)
        Check interval (Tuple of initial time and end time)
    ref_interval : (Time, Time)
        Reference interval (Tuple of initial time and end time)

    Returns
    -------
    bool
        True if interval begin and end are (almost) equal, False otherwise

    """
    begin_equal = abs(interval[0] - ref_interval[0]) < _EPS_TIME
    end_equal = abs(interval[1] - ref_interval[1]) < _EPS_TIME

    if begin_equal and end_equal:
        return True
    else:
        return False
