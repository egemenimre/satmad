# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Time interval module.

`TimeInterval` class stores time intervals and `TimeIntervalList` class stores
lists of `TimeInterval` objects.
"""

import portion as p
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.units import Quantity

_EPS_TIME = 10 * u.ns
"""Allowable time threshold, this much 'out of bounds' is allowed when assuming two
instances of time are *practically* equal. This helps with floating point artifacts
such as round-off errors."""


class TimeInterval:
    """
    Represent and manipulate a single time interval.

    This is a thin wrapper around the  :class:`portion.interval.Interval` class
    from `portion` package (for the Atomic intervals), using Astropy `Time` classes
    under the hood.


    Note that:

    * only the first element of `start_time` and `end_time` objects are used,
      even if they are defined as arrays.

    * `end_time` should be greater than (or later than) `start_time`.
      Zero duration time intervals are not allowed and will raise `ValueError`.

    Parameters
    ----------
    start_time : Time or TimeInterval
        Initial time to mark the start of the interval (only the first
        instance contained in the `Time` object is used). If specified as TimeInterval,
        it is copied into this and `end_time` is ignored.
    end_time : Time, TimeDelta, Quantity or None
        End time or duration to mark the end of the interval
        (only the first instance contained in the `Time` or `TimeDelta` object is used)
        None is acceptable only when `start_time` is defined as a `TimeInterval`
    replicate : bool
        `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False` to
        use a shallow copy to save memory
    start_inclusive : bool
        True if the start of the interval is inclusive (closed), False if exclusive
        (open)
    end_inclusive : bool
        True if the start of the interval is inclusive (closed), False if exclusive
        (open)

    Raises
    ------
    ValueError
        Raised if `end_time` is not an instance of `Time` or `TimeDelta`. Also raised
        if `end_time` <= `start_time`.
    """

    def __init__(
        self,
        start_time,
        end_time,
        replicate=False,
        start_inclusive=True,
        end_inclusive=True,
    ):

        start_interval = None
        end_interval = None
        if isinstance(start_time, TimeInterval):
            # start_time is a valid TimeInterval, just copy its elements
            start_inclusive = start_time.p_interval.left
            end_inclusive = start_time.p_interval.right

            if replicate:
                start_interval = start_time.start.replicate()
                end_interval = start_time.end.replicate()
            else:
                start_interval = start_time.start
                end_interval = start_time.end
        else:
            # Check whether init and end times are scalar or array
            if start_time.isscalar:
                start_tmp = start_time
            else:
                start_tmp = start_time[0]

            if end_time.isscalar:
                end_tmp = end_time
            else:
                end_tmp = end_time[0]

            # Check the type of end_time and fill values accordingly
            if isinstance(end_tmp, TimeDelta):
                # end time is a TimeDelta (duration)

                if replicate:
                    start_interval = start_tmp.replicate()
                    end_interval = start_tmp + end_tmp
                else:
                    start_interval = start_tmp
                    end_interval = start_tmp + end_tmp

            elif isinstance(end_tmp, Quantity):
                # end time is a Quantity, but check whether it is a valid Time Quantity
                if end_tmp.unit.decompose().physical_type == "time":
                    start_interval = start_tmp
                    end_interval = start_tmp + end_tmp

            elif isinstance(end_tmp, Time) and not isinstance(end_tmp, TimeDelta):
                # end time is a Time

                # replicate or shallow copy start and end values - use the first
                # time instances only
                if start_tmp < end_tmp:
                    # end time is later than start time - no probs
                    if replicate:
                        start_interval = start_tmp.replicate()
                        end_interval = end_tmp.replicate()
                    else:
                        start_interval = start_tmp
                        end_interval = end_tmp
                else:
                    # end time is earlier than start time - raise error
                    raise ValueError(
                        f"End time ({end_tmp.isot}) is earlier than start time"
                        f"({start_tmp.isot})"
                    )
            else:
                raise ValueError(
                    f"End time is an instance of {end_tmp.__class__()}, "
                    f"only Time, Quantity (Temporal) or TimeDelta classes are allowed."
                )

        # Initialise the interval
        self._interval = p.closed(start_interval, end_interval).replace(
            left=start_inclusive, right=end_inclusive
        )

        # Raise exception if the interval duration is too short
        # Duration guaranteed to be positive thanks to the juggling above
        if self.duration <= _EPS_TIME:
            raise ValueError("Duration of the interval is negative or zero.")

    def is_in_interval(self, time):
        """
        Checks whether the requested time is within the time interval.

        Parameters
        ----------
        time : Time
            Time to be checked

        Returns
        -------
        bool
            True if time is within the reference interval, False otherwise
        """
        # check upper and lower boundaries for tolerance
        if _are_times_almost_equal(self.start, time):
            # time at starting edge - is edge closed?
            if self._interval.left:
                return True
            else:
                return False

        if _are_times_almost_equal(self.end, time):
            # time at end edge - is edge closed?
            if self._interval.right:
                return True
            else:
                return False

        # Time not edges, do a regular check
        return self._interval.contains(time)

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
        start_equal = _are_times_almost_equal(interval.start, self.start)
        end_equal = _are_times_almost_equal(interval.end, self.end)

        if start_equal and end_equal:
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
        if self.is_in_interval(interval.start) and self.is_in_interval(interval.end):
            return True
        else:
            return False

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
        intersection = self._interval.intersection(interval._interval)

        if intersection.empty:
            # There is absolutely no intersection
            return False

        if _are_times_almost_equal(intersection.upper, intersection.lower):
            # intersection below tolerance - practically empty intersection
            return False
        else:
            # see what the underlying function says
            return self._interval.overlaps(interval._interval)

    def intersect(self, interval):
        """
        Intersection operator for a time interval and this time interval.

        Returns a new interval that is the Intersection of two intervals,
        or None if there is no intersection.

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
            intersection = self._interval.intersection(interval._interval)

            return _create_interval_from_portion(intersection)
        else:
            return None

    def union(self, interval):
        """
        Union operator for a time interval and this time interval.

        Returns a new interval that is the Union of two intervals,
        or None if there is no intersection.

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
            union = self._interval.union(interval._interval)
            return _create_interval_from_portion(union)
        else:
            return None

    def expand(
        self,
        start_delta=0,
        end_delta=0,
        replicate=False,
        start_inclusive=True,
        end_inclusive=True,
    ):
        """
        Expands (or shrinks) the interval.

        Generates a new, expanded (or shrunk) `TimeInterval`, where:

        - new interval start: interval_start - start_delta
        - new interval end: interval_end + end_delta

        Negative start and/or end times are possible (to shrink the interval),
        though values ending in a negative interval will throw a `ValueError`.
        This method can be used to modify the ends of the interval (open or closed)
        as well.

        Parameters
        ----------
        start_delta : TimeDelta
            The delta time to expand the start of the interval
            (or to shrink, with negative values)
        end_delta : TimeDelta
            The delta time to expand the end of the interval
            (or to shrink, with negative values)
        replicate : bool
            `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False`
            to use a shallow copy to save memory
        start_inclusive : bool
            True if the start of the new interval is inclusive (closed), False if
            exclusive (open)
        end_inclusive : bool
            True if the start of the new interval is inclusive (closed), False if
            exclusive (open)

        Returns
        -------
        TimeInterval
            A new `TimeInterval` that is the result of the requested change

        Raises
        ------
        ValueError
            Raised if the requested expansion results in a negative duration interval
        """
        start = self.start - start_delta
        end = self.end + end_delta

        duration = end - start

        if duration > _EPS_TIME:
            return TimeInterval(
                start,
                end,
                replicate=replicate,
                start_inclusive=start_inclusive,
                end_inclusive=end_inclusive,
            )
        else:
            raise ValueError(
                "Duration of the expanded/shrunk interval is negative or zero."
            )

    @property
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
        return self.end - self.start

    @property
    def start(self) -> Time:
        """Returns the start time of the interval."""
        return self._interval.lower

    @property
    def end(self) -> Time:
        """Returns the end time of the interval."""
        return self._interval.upper

    def __str__(self):
        txt = ""
        if self._interval.left:
            txt += "["
        else:
            txt += "("

        txt += f" {str(getattr(self.start, self.start.format))}"
        txt += f"  {str(getattr(self.end, self.end.format))} "

        if self._interval.right:
            txt += "]"
        else:
            txt += ")"

        return txt

    @property
    def p_interval(self):
        """
        Returns the underlying `Interval` object.

        .. warning:: Most users will not need to access this object. Intended for
            developer use only.
        """
        return self._interval


def _are_times_almost_equal(t1, t2):
    """
    Checks whether two time instances are are equal to a tolerance given by `_EPS_TIME`.

    Parameters
    ----------
    t1, t2 : Time
        Times to be checked

    Returns
    -------
    True if times are almost equal, False otherwise

    """
    return abs(t1 - t2) < _EPS_TIME


def _create_interval_from_portion(interval, replicate=False):
    """
    Create a new `TimeInterval` from a given :class:`portion.interval.Interval`
    instance.

    Parameters
    ----------
    interval : Interval
        input `Interval` instance
    replicate : bool
        `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False`
        to use a shallow copy to save memory

    Returns
    -------
    TimeInterval
        New `TimeInterval` instance
    """
    return TimeInterval(
        interval.lower,
        interval.upper,
        start_inclusive=interval.left,
        end_inclusive=interval.right,
        replicate=replicate,
    )


class TimeIntervalList:
    """
    Represent and manipulate time intervals.

    This is a thin wrapper around the :class:`portion.interval.Interval` class,
    using Astropy `Time` classes under the hood.

    `start_valid` and `end_valid` values are used to mark the start and end of this
    list of time intervals. If they are not specified, the beginning and end points
    of the list of `TimeInterval` instances are used.

    Parameters
    ----------
    intervals : list[TimeInterval]
        List of intervals
    start_valid : Time or TimeInterval
        Time at which the validity of this `TimeInterval` starts (only the first
        instance contained in the `Time` object is used). If specified as TimeInterval,
        it is copied into this and `end_valid` is ignored (the only case where `None`
        is acceptable).
    end_valid : Time, TimeDelta or None
        Time or duration at which the validity of this `TimeInterval` ends (only the
        first instance contained in the `Time` object is used). None is acceptable
        only when `start_valid` is defined as a `TimeInterval`
    replicate : bool
        `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False` to
        use a shallow copy to save memory

    Raises
    ------
    ValueError
        Raised if `init_times` and `end_times` objects sizes mismatch.
    """

    def __init__(self, intervals, start_valid=None, end_valid=None, replicate=False):

        self._intervals: list = []

        if intervals:
            # if start_times is None, then there is no time interval defined

            # Fill the `Interval` list and merge as necessary
            p_intervals = self._to_p_intervals(intervals)

            # Fill the atomic `TimeInterval` objects using the merged list
            self._intervals = self._to_time_intervals(p_intervals, replicate)

        # Init range of validity
        self._valid_interval = self.__init_validity_rng(
            start_valid, end_valid, replicate
        )

    def __init_validity_rng(self, start_valid=None, end_valid=None, replicate=False):
        """
        Initialises the beginning and end of the range of validity.

        `start_valid` and `end_valid` values are used to mark the start and end of this
        list of time intervals (and are Inclusive or Closed by default).
        If they are not specified, the beginning and end points
        of the `TimeInterval` are used.

        Parameters
        ----------
        start_valid : Time or TimeInterval
            Time at which the validity of this `TimeInterval` starts (only the first
            instance contained in the `Time` object is used). If specified as
            `TimeInterval`, it is copied into this and `end_valid` is ignored
            (the only case where `None` is acceptable for it).
        end_valid : Time or TimeDelta or None
            Time or duration at which the validity of this `TimeInterval` ends (only the
            first instance contained in the `Time` object is used). None is acceptable
            only when `start_time` is defined as a `TimeInterval`
        replicate : bool
            `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False`
            to use a shallow copy to save memory

        Returns
        -------
        TimeInterval
            A new `TimeInterval` instance containing the start and end of validity
        """

        start_inclusive = True
        end_inclusive = True

        if isinstance(start_valid, TimeInterval):
            # `start_valid` is already a TimeInterval, just copy it
            return TimeInterval(start_valid, None, replicate)

        if start_valid is None:
            # No init time for validity defined, use the first interval init
            start_valid_tmp = self.get_interval(0).start
            start_inclusive = self.get_interval(0).p_interval.left
        else:
            start_valid_tmp = start_valid

        if end_valid is None:
            # No end time for validity defined, use the last interval end
            end_valid_tmp = self.get_interval(-1).end
            end_inclusive = self.get_interval(0).p_interval.right
        else:
            end_valid_tmp = end_valid

        # Generate the `TimeInterval`
        return TimeInterval(
            start_valid_tmp,
            end_valid_tmp,
            start_inclusive=start_inclusive,
            end_inclusive=end_inclusive,
            replicate=replicate,
        )

    def is_in_interval(self, time):
        """
        Checks whether the requested time is within the time interval list.

        Parameters
        ----------
        time : Time
            Time to be checked

        Returns
        -------
        bool
            `True` if time is within the interval list, `False` otherwise. Also returns
            `False` if requested time is outside validity interval.
        """
        # Is time within validity interval?
        if not self.valid_interval.is_in_interval(time):
            return False

        # Are there any events that contain this time instant?
        for interval in self._intervals:
            if interval.is_in_interval(time):
                return True

        # If we are here, then no interval contains the time
        return False

    def is_intersecting(self, interval):
        """
        Checks whether the requested interval intersects (or is contained within)
        the interval list.

        Parameters
        ----------
        interval : TimeInterval
            Time interval to be checked

        Returns
        -------
        bool
            `True` if check interval intersects with the interval list,
            `False` otherwise
        """
        if len(self._intervals) == 0:
            # No interval present in the list, hence no intersection
            return False

        # While not very elegant, loop through the interval list to check
        # for intersections
        intersect_intervals = [
            interval_member
            for interval_member in self.intervals
            if interval_member.is_intersecting(interval)
        ]

        if len(intersect_intervals) > 0:
            return True
        else:
            return False

    def intersect(self, interval):
        """
        Intersection operator for a time interval and this time interval list.

        Returns a new interval that is the Intersection of the interval and
        the time interval list, or None if there is no intersection.

        Parameters
        ----------
        interval : TimeInterval
            Time interval to be checked

        Returns
        -------
        TimeInterval
           A new interval that is the Intersection of the interval and
           the time interval list, or None if there is no intersection.

        """
        if len(self.intervals) == 0:
            # No interval present in the list, hence no intersection possible
            return None

        # While not very elegant, loop through the interval list to check
        # for intersections
        intersect_intervals = [
            interval_member.intersect(interval)
            for interval_member in self.intervals
            if interval_member.is_intersecting(interval)
        ]

        if len(intersect_intervals) > 0:
            # There can be only a single intersection
            return intersect_intervals[0]
        else:
            # no intersection
            return None

    def intersect_list(self, interval_list):
        """
        Intersection operator for a time interval list and this time interval list.

        Returns a new interval list that is the Intersection of `interval_list` and
        this time interval list, or None if there is no intersection even in the
        validity intervals.

        Parameters
        ----------
        interval_list : TimeIntervalList
            Time interval list to be checked

        Returns
        -------
        TimeIntervalList
           A new interval list that is the Intersection of `interval_list` and
           this time interval list, or None if there is no intersection.

        """
        if not self.valid_interval.is_intersecting(interval_list.valid_interval):
            # Validity intervals do not intersect, hence no intersection possible
            return None

        common_valid_interval = self.valid_interval.intersect(
            interval_list.valid_interval
        )

        # Compute the portion intervals
        p_self_intervals = TimeIntervalList._to_p_intervals(self.intervals)
        p_other_intervals = TimeIntervalList._to_p_intervals(interval_list.intervals)

        # Do the Intersection
        p_final = p_self_intervals.intersection(p_other_intervals)

        return TimeIntervalList(
            self._to_time_intervals(p_final), start_valid=common_valid_interval
        )

    def union(self, interval):
        """
        Union operator for a time interval and this time interval list.

        Returns a new interval that is the Union of the interval and
        the time interval list, or None if there is no intersection.

        Parameters
        ----------
        interval : TimeInterval
            Time interval to be checked

        Returns
        -------
        TimeInterval
           A new interval that is the Union of the interval and
           the time interval list, or None if there is no intersection.

        """
        if len(self.intervals) == 0:
            # No interval present in the list, hence no intersection possible
            return None

        # While not very elegant, loop through the interval list to check
        # for intersections
        union_intervals = [
            interval_member.union(interval)
            for interval_member in self.intervals
            if interval_member.is_intersecting(interval)
        ]

        if len(union_intervals) > 0:
            # There can be only a single intersection
            return union_intervals[0]
        else:
            # no intersection
            return None

    def union_list(self, interval_list):
        """
        Union operator for a time interval list and this time interval list.

        Returns a new interval list that is the Union of `interval_list` and
        this time interval list, or None if there is no intersection even in the
        validity intervals.

        Parameters
        ----------
        interval_list : TimeIntervalList
            Time interval list to be checked

        Returns
        -------
        TimeIntervalList
           A new interval list that is the Union of `interval_list` and
           this time interval list, or None if there is no intersection.

        """
        if not self.valid_interval.is_intersecting(interval_list.valid_interval):
            # Validity intervals do not intersect, hence no intersection possible
            return None

        common_valid_interval = self.valid_interval.intersect(
            interval_list.valid_interval
        )

        # Compute the portion intervals
        p_self_intervals = TimeIntervalList._to_p_intervals(self.intervals)
        p_other_intervals = TimeIntervalList._to_p_intervals(interval_list.intervals)

        # Do the Union
        p_union = p_self_intervals.union(p_other_intervals)

        # Reduce union to the common interval
        p_common = common_valid_interval.p_interval
        p_final = p_union.intersection(p_common)

        return TimeIntervalList(
            self._to_time_intervals(p_final), start_valid=common_valid_interval
        )

    def invert(self, replicate=False):
        """
        Creates an *inverted* (or complement) copy of this time interval list, while
        keeping the same validity range.

        For example, for a single interval of `[t0, t1]` in a validity interval
        `[T0,T1]`, the inverted interval list would be `[T0,t0]` and `[t1,T1]`. If
        there are no intervals, the inverse becomes the entire validity interval.

        Parameters
        ----------
        replicate : bool
            `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False`
            to use a shallow copy to save memory

        Returns
        -------
        A new `TimeIntervalList` that has the same validity range but the individual
        intervals are inverted.
        """
        # Convert `TimeInterval` list to `Interval`
        p_interval = self._to_p_intervals(self._intervals)

        # Do the inversion
        p_int_inverted = ~p_interval

        # Fix the ends as necessary - no inf allowed
        p_validity = self._to_p_intervals(self.valid_interval)
        p_int_inverted = p_int_inverted.intersection(p_validity)

        # Generate the `TimeInterval` list
        intervals = self._to_time_intervals(p_int_inverted, replicate=False)

        # Create the `TimeIntervalList` object
        return TimeIntervalList(
            intervals, start_valid=self.valid_interval, replicate=replicate
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
        TimeInterval
            `TimeInterval` corresponding to the index

        Raises
        ------
        IndexError
            Requested index is out of bounds

        """
        return self._intervals[index]

    @property
    def valid_interval(self) -> TimeInterval:
        """
        Gets the time interval of validity for the `TimeIntervalList`.
        """
        return self._valid_interval

    @property
    def intervals(self):
        """
        Gets the time intervals within this `TimeIntervalList`.
        """
        return self._intervals

    @staticmethod
    def _to_time_intervals(p_intervals, replicate=False):
        """
        Converts a `Interval` list to a `TimeInterval` list.

        Parameters
        ----------
        p_intervals : Interval
            List of intervals
        replicate : bool
            `True` to replicate (deep copy) the `Time` or `TimeDelta` objects, `False` to
            use a shallow copy to save memory

        Returns
        -------
        list[TimeInterval]
            `TimeInterval` object with the list of time intervals
        """
        intervals: list = []

        # Fill the atomic `TimeInterval` objects using the merged list
        for p_interval in p_intervals:
            intervals.append(
                _create_interval_from_portion(p_interval, replicate=replicate)
            )

        return intervals

    @staticmethod
    def _to_p_intervals(intervals):
        """
        Converts a `TimeInterval` instance or a list to an `Interval` list
        (portion library objects).

        This is usually done to merge and simplify the elements of the list.

        Parameters
        ----------
        intervals : TimeInterval or list[TimeInterval]
            List of intervals

        Returns
        -------
        `Interval` object with the list of time intervals
        """
        # Fill the `Interval` list and merge as necessary
        p_intervals = p.empty()
        if isinstance(intervals, list):
            for interval in intervals:
                p_intervals = p_intervals.union(interval.p_interval)
        else:
            # intervals object is a single TimeInterval
            p_intervals = intervals.p_interval

        return p_intervals

    def __str__(self):
        txt = ""
        if self._intervals:
            # List not empty
            for interval in self._intervals:
                txt += str(interval) + "\n"
        else:
            txt = "Time interval list is empty."

        return txt
