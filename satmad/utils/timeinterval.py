# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) ${YEAR} Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Time interval module.

`TimeInterval` class stores time intervals and `TimeIntervalList` class stores
lists of `TimeInterval` objects.
"""

from astropy import units as u
from astropy.time import Time, TimeDelta

_EPS_TIME = 10 * u.us
"""Allowable time threshold, this much 'out of bounds' is allowed when assuming two
instances of time are *practically* equal. This helps with floating point artifacts
such as round-off errors."""


class TimeInterval:
    """
    Represent and manipulate a single time interval.

    Uses Astropy `Time` classes under the hood. Note that only the first element
    of `start_time` and `end_time` objects are used, even if they are defined as arrays.
    Note that `end_time` should be greater than (or later than) `start_time`. Zero
    duration time intervals are not allowed and will raise `ValueError`.

    Parameters
    ----------
    start_time : Time
        Initial time to mark the start of the interval (only the first
        instance contained in the `Time` object is used)
    end_time : Time or TimeDelta
        End time or duration to mark the end of the interval
        (only the first instance contained in the `Time` or `TimeDelta` object is used)
    replicate : bool
        `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False` to
        use a shallow copy to save memory

    Raises
    ------
    ValueError
        Raised if `end_time` is not an instance of `Time` or `TimeDelta`. Also raised
        if `end_time` <= `start_time`.
    """

    def __init__(self, start_time, end_time, replicate=False):

        # Check whether init and end times are scalar or array
        if start_time.isscalar:
            start = start_time
        else:
            start = start_time[0]

        if end_time.isscalar:
            end = end_time
        else:
            end = end_time[0]

        # Check the type of end_time and fill values accordingly
        if isinstance(end_time, TimeDelta):
            # end time is a TimeDelta (duration)

            if replicate:
                self._start = start.replicate()
            else:
                self._start = start
            self._end = start + end

        elif isinstance(end_time, Time) and not isinstance(end_time, TimeDelta):
            # end time is a Time

            # replicate or shallow copy start and end values - use the first
            # time instances only
            if start < end:
                # end time is later than start time - no probs
                if replicate:
                    self._start = start.replicate()
                    self._end = end.replicate()
                else:
                    self._start = start
                    self._end = end
            else:
                # end time is earlier than start time - raise error
                raise ValueError(
                    f"End time ({end.isot}) is earlier than start time ({start.isot})"
                )
        else:
            raise ValueError(
                f"End time is an instance of {end.__class__()}, "
                f"only Time or TimeDelta classes are allowed."
            )

        # Raise exception if the interval duration is too short
        # Duration guaranteed to be positive thanks to the juggling above
        if self.duration() <= _EPS_TIME:
            raise ValueError("Duration of the interval is negative or zero.")

    def is_within_interval(self, time, inclusive=True):
        """
        Checks whether the requested time is within the time interval.

        - If the `inclusive` flag is True: interval_start <= time <= interval_end
        - If the `inclusive` flag is False: interval_start < time < interval_end

        Parameters
        ----------
        time : Time
            Time to be checked
        inclusive : bool
            True if max and min boundaries are included in the check, False otherwise

        Returns
        -------
        bool
            True if time is within the reference interval, False otherwise
        """
        if inclusive:
            if self._start.value <= time.value <= self._end.value:
                return True
            else:
                return False
        else:
            if self._start.value < time.value < self._end.value:
                return True
            else:
                return False

    def is_equal(self, interval):
        """
        Checks whether two intervals are (almost) equal in value.

        If the two start or end values are as close as `_EPS_TIME`, then they are
        assumed to be equal in value.

        Parameters
        ----------
         interval : TimeInterval
            Time interval to be checked

        Returns
        -------
        bool
            True if interval start and end are (almost) equal, False otherwise

        """
        start_equal = abs(interval._start - self._start) < _EPS_TIME
        end_equal = abs(interval._end - self._end) < _EPS_TIME

        if start_equal and end_equal:
            return True
        else:
            return False

    def expand(self, start_delta=0, end_delta=0):
        """
        Expands (or shrinks) the interval.

        Generates a new `TimeInterval`, where:

        - new interval start: interval_start - start_delta
        - new interval end: interval_end + end_delta

        Negative start and/or end times are possible (to shrink the interval),
        though values ending in a negative interval will throw a `ValueError`.

        Parameters
        ----------
        start_delta : TimeDelta
            The delta time to expand the start of the interval
            (or to shrink, with negative values)
        end_delta : TimeDelta
            The delta time to expand the end of the interval
            (or to shrink, with negative values)
        Returns
        -------
        TimeInterval
            A new `TimeInterval` that is the result of the requested change

        Raises
        ------
        ValueError
            Raised if the requested expansion results in a negative duration interval
        """
        start = self._start - start_delta
        end = self._end + end_delta

        duration = end - start

        if duration > _EPS_TIME:
            return TimeInterval(start, end)
        else:
            raise ValueError(
                "Duration of the expanded/shrunk interval is negative or zero."
            )

    def duration(self):
        """
        Computes the duration of the interval.

        Parameters
        ----------

        Returns
        -------
        duration : TimeDelta
            Duration of the interval (always positive)
        """
        return self._end - self._start

    def is_intersecting(self, interval):
        """
        Checks whether the requested interval intersects (or is contained within)
        the reference interval.

        Parameters
        ----------
        interval : TimeInterval
            Time interval to be checked

        Returns
        -------
        bool
            True if check interval intersects with the reference interval,
            False otherwise
        """
        if self.contains(interval):
            # interval is contained or the same as this (self) interval
            return True
        elif self.is_within_interval(interval._start, False) or self.is_within_interval(
            interval._end, False
        ):
            # at least interval_start or interval_end is within this interval
            return True
        else:
            return False

    def contains(self, interval):
        """
        Checks whether the requested interval is contained within this (self) interval.

        Parameters
        ----------
        interval : TimeInterval
            Time interval to be checked

        Returns
        -------
        bool
            True if check interval is contained within this interval,
            False otherwise
        """
        if self.is_within_interval(interval._start, True) and self.is_within_interval(
            interval._end, True
        ):
            return True
        else:
            return False

    def intersect(self, interval):
        """
        Intersection operator for a time interval and this time interval.

        Returns a new interval that is
        the Intersection of two intervals, or None if there is no intersection.

        Parameters
        ----------
        interval : TimeInterval
            Time interval to be checked

        Returns
        -------
        TimeInterval
            A new interval that is the Intersection of two intervals,
            or `None` if there is no intersection

        """
        if self.is_intersecting(interval):

            if interval._start < self._start:
                # other interval starts before this interval - set as
                # this interval start
                new_init = self._start
            else:
                # other interval starts after this interval - set as
                # other interval start
                new_init = interval._start

            if self._end < interval._end:
                # other interval ends after this interval - set as this interval end
                new_end = self._end
            else:
                # other interval ends before this interval - set as other interval end
                new_end = interval._end

            return TimeInterval(new_init, new_end)
        else:
            # No intersection found
            return None

    def union(self, interval):
        """
        Union operator for a time interval and this time interval.

        Returns a new interval that is
        the Union of two intervals, or None if there is no intersection.

        Parameters
        ----------
        interval : TimeInterval
            Time interval to be checked

        Returns
        -------
        TimeInterval
            A new interval that is the Union of two intervals,
            or `None` if there is no intersection

        """
        if self.is_intersecting(interval):

            if interval._start < self._start:
                # other interval starts before this interval - set as
                # other interval start
                new_init = interval._start
            else:
                # other interval starts after this interval - set as this interval start
                new_init = self._start

            if self._end < interval._end:
                # other interval ends after this interval - set as other interval end
                new_end = interval._end
            else:
                # other interval ends before this interval - set as this interval end
                new_end = self._end

            return TimeInterval(new_init, new_end)
        else:
            # No intersection found
            return None

    @property
    def start(self) -> Time:
        """Returns the start time of the interval."""
        return self._start

    @property
    def end(self) -> Time:
        """Returns the end time of the interval."""
        return self._end

    def __str__(self):
        txt = f"[ {str(getattr(self.start, self.start.format))}"
        txt += f"  {str(getattr(self.end, self.end.format))} ]"

        return txt


class TimeIntervalList:
    """
    Represent and manipulate time intervals.

    Uses Astropy `Time` classes under the hood. Note that the number of `init_times`
    and `end_times` objects should match.

    `start_valid` and `end_valid` values are used to mark the start and end of this
    list of time intervals. If they are not specified, the beginning and end points
    of the `TimeInterval` are used.

    Parameters
    ----------
    start_times : Time
        List of initial times to mark the beginning of each interval
    end_times : Time or TimeDelta
        List of `TimeDelta` objects or end times, to mark the end of each interval
    start_valid : Time
        Time at which the validity of this `TimeInterval` starts (only the first
        instance contained in the `Time` object is used)
    end_valid : Time or TimeDelta
        Time or duration at which the validity of this `TimeInterval` ends (only the
        first instance contained in the `Time` object is used)
    replicate : bool
        `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False` to
        use a shallow copy to save memory

    Raises
    ------
    ValueError
        Raised if `init_times` and `end_times` objects sizes mismatch.
    """

    # TODO check for out of validity errors
    #  when doing Union with an external interval
    #
    # TODO  take union etc methods to within class
    #  with thin wrappers and checks and multiple items

    def __init__(
        self, start_times, end_times, start_valid=None, end_valid=None, replicate=False
    ):

        end_times_tmp = None
        # Fill end times depending on the type
        if isinstance(end_times, TimeDelta):
            # End times are of type `TimeDelta`, make sure array sizes match
            # - or that there is a single duration
            try:
                end_times_tmp = start_times + end_times
            except Exception as err:
                raise ValueError(
                    f"Array size mismatch between "
                    f"Initial Times (size: {len(start_times)}) and "
                    f"Durations (size: {len(end_times)})."
                ) from err

        elif isinstance(end_times, Time) and not isinstance(end_times, TimeDelta):
            # End times are of type `Time`, make sure array sizes match
            if len(end_times) == len(start_times):
                end_times_tmp = end_times
            else:
                raise ValueError(
                    f"Array size mismatch between "
                    f"Initial Times (size: {len(start_times)}) and "
                    f"End Times (size: {len(end_times)})."
                )

        # Fill the `TimeInterval` list
        self._interval_list: list = []
        for i, start_time in enumerate(start_times):
            self._interval_list.append(TimeInterval(start_times[i], end_times_tmp[i]))

        # Init range of validity
        self._valid_interval = self.__init_validity_rng(
            start_valid, end_valid, replicate
        )

    def __init_validity_rng(self, start_valid=None, end_valid=None, replicate=False):
        """
        Initialises the beginning and end of the range of validity.

        `start_valid` and `end_valid` values are used to mark the start and end of this
        list of time intervals. If they are not specified, the beginning and end points
        of the `TimeInterval` are used.

        Parameters
        ----------
        start_valid : Time
            Time at which the validity of this `TimeInterval` starts (only the first
            instance contained in the `Time` object is used)
        end_valid : Time or TimeDelta
            Time or duration at which the validity of this `TimeInterval` ends (only the
            first instance contained in the `Time` object is used)
        replicate : bool
            `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False`
            to use a shallow copy to save memory

        Returns
        -------
        TimeInterval
            A new `TimeInterval` instance containing the start and end of validity
        """
        if start_valid is None:
            # No init time for validity defined, use the first interval init
            start_valid_tmp = self.get_interval(0).start
        else:
            start_valid_tmp = start_valid

        if end_valid is None:
            # No end time for validity defined, use the last interval end
            end_valid_tmp = self.get_interval(-1).end
        else:
            end_valid_tmp = end_valid

        return TimeInterval(start_valid_tmp, end_valid_tmp, replicate)

    def get_interval(self, index):
        """
        Gets the time interval for the given index.

        Parameters
        ----------
        index : int
            requested index

        Returns
        -------
        TimeInterval
            `TimeInterval` corresponding to the index

        Raises
        ------
        IndexError
            Requested index is out of bounds

        """
        return self._interval_list[index]

    def validity_interval(self):
        """
        Gets the time interval of validity for the `TimeIntervalList`.

        Returns
        -------
        TimeInterval
            `TimeInterval` of interval of validity

        """
        return self._valid_interval

    def __str__(self):
        txt = ""
        for interval in self._interval_list:
            txt += str(interval) + "\n"

        return txt
