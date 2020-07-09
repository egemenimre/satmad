Time Intervals and Time Interval Lists
======================================

Introduction
------------
Time intervals are critical to define the start and end of certain events such as
start and end of communications with a groundstation, entering and exiting
the eclipse or start and end of a thruster firing.
This is particularly useful when two intervals (or sets of intervals)
can be evaluated through operations such as *union* and *intersection*. This enables
us to answer questions such as "What are the time intervals where thruster firings
occur during communications?" (an intersection operation between
'intervals of thruster firings' and
'communications interval lists') or "When can I see a satellite at night?" (an
intersection operation between intervals of 'satellite above horizon', 'sun below
horizon' and 'satellite not in eclipse').


The :mod:`.timeinterval` module provides the basic time interval functionality with
the :class:`.TimeInterval` class i.e.,
a time interval with a start and end time/date, using the high precision
:class:`astropy.time.Time` classes under the hood to represent time and
:class:`portion.interval.Interval` class to manage and manipulate the time intervals.

A :class:`.TimeInterval` can interact
with other intervals through :meth:`.TimeInterval.union` and :meth:`.TimeInterval.intersect`
methods. They can change their size through :meth:`.TimeInterval.expand` and they can
check whether they contain (:meth:`.TimeInterval.contains`) or intersect with
(:meth:`.TimeInterval.is_intersecting`) another time interval.

A list of such time intervals constitute :class:`.TimeIntervalList` class. A list also
has a start and end of validity. This usually marks the start and end of an analysis.
For example, a communications list that is valid for one day and containing no time intervals
would mean that there are no communication opportunities for that day. The list can simply be
inverted (:meth:`.TimeIntervalList.invert`) to get a list of 'no communication duration',
which would then show a list with a single :class:`.TimeInterval` that spans the entire
duration of validity.



Using the Basic :class:`.TimeInterval` Class
----------------------------------------------

A :class:`.TimeInterval` class can be simply initialised with a start time
and either with an end time (:class:`astropy.time.Time`) or with a duration (
:class:`astropy.time.TimeDelta`). These start and end times can be retrieved
by the properties :meth:`.TimeInterval.start` and :meth:`.TimeInterval.end`.

.. code-block:: python
    :linenos:

    from astropy.time import Time, TimeDelta
    from satmad.utils.timeinterval import TimeInterval

    interval_with_end_time = TimeInterval(
        Time("2020-04-11T00:00:00.000", scale="utc"),
        Time("2020-04-11T00:10:00.000", scale="utc"),
    )
    interval_with_duration = TimeInterval(
        Time("2020-04-11T00:00:00", scale="utc"),
        TimeDelta(60.0, format='sec'),
    )

The resulting time intervals can be quickly shown as:

    >>> str(interval_with_end_time)
    '[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n'
    >>> str(interval_with_duration)
    '[ 2020-04-11T00:00:00.000  2020-04-11T00:01:00.000 ]\n'

.. note:: The end time of the interval should be later than the start time.
    Otherwise a `ValueError` will be raised.

The :class:`.TimeInterval` class can answer some questions:

- :meth:`.TimeInterval.is_in_interval`: Is a given time within this interval?
- :meth:`.TimeInterval.is_equal`: Is a given interval equal to this interval?
- :meth:`.TimeInterval.is_intersecting`: Does a given interval have an intersection with this interval?
- :meth:`.TimeInterval.contains`: Does a given interval contain this interval?
- :meth:`.TimeInterval.duration`: What is the duration of this interval?

.. code-block:: python

    >>> interval_with_end_time.is_in_interval(Time("2020-04-11T00:05:00.000", scale="utc"))
    True
    >>> interval_with_end_time.is_equal(TimeInterval(Time("2020-04-09T00:00:00", scale="utc"), TimeDelta(600.0, format="sec")))
    True
    >>> interval_with_end_time.is_intersecting(TimeInterval(Time("2020-04-11T00:05:00", scale="utc"), TimeDelta(600.0, format="sec")))
    True
    >>> interval_with_end_time.contains(TimeInterval(Time("2020-04-11T00:05:00", scale="utc"), TimeDelta(60.0, format="sec")))
    True
    >>> interval_with_end_time.duration().sec
    599.9999999999931


The intervals can be expanded or shrunk through the :meth:`.TimeInterval.expand` method and
defining positive or negative :class:`astropy.time.TimeDelta` values to modify the start and
end times of the interval.

.. code-block:: python

    >>> str(interval_with_end_time)
    '[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n'
    >>> expanded = interval_with_end_time.expand(start_delta=TimeDelta(60.0, format="sec"), end_delta=TimeDelta(-120.0, format="sec"))
    >>> str(expanded)
    '[ 2020-04-10T23:59:00.000  2020-04-11T00:08:00.000 ]\n'

Time intervals can be subjected to an intersection (a new `TimeInterval` that is the intersection of
two intervals) or union (a new `TimeInterval` that is the union of
two intervals) operator. These operators are possible only when these two intervals have
some intersection - otherwise the result will be a `None`.

.. code-block:: python

    >>> str(interval_with_end_time.union(expanded))
    '[ 2020-04-10T23:59:00.000  2020-04-11T00:10:00.000 ]\n'
    >>> str(interval_with_end_time.intersect(expanded))
    '[ 2020-04-11T00:00:00.000  2020-04-11T00:08:00.000 ]\n'


List of time intervals: The :class:`.TimeIntervalList` Class
------------------------------------------------------------------

The :class:`.TimeIntervalList` class stores the :class:`.TimeInterval` objects as well as
another `TimeInterval` to represent the bounds of the validity of this list. If this validity
interval is not defined explicitly, then it is assumed to start with the beginning of the first
`TimeInterval` and end with the end of the final `TimeInterval`.

Any interval within the list can be queried through :meth:`.TimeIntervalList.get_interval`
method. Similarly, the `TimeInterval` that keeps the interval of validity can be queried
through :meth:`.TimeIntervalList.validity_interval` method.

The :class:`.TimeIntervalList` usually will not be generated explicitly by a user, except,
for example, as an external constraint such as the durations when a groundstation is not
available. Usually such lists are results of certain analyses such as eclipse intervals
for a location on ground or different attitude profiles for a satellite.

The `TimeIntervalList` will yield a new, inverted (or complementing) version of itself
through the :meth:`.TimeIntervalList.invert` method. For example, for a single interval of
`[t0, t1]` in a validity interval `[T0,T1]`, the inverted interval list would be `[T0,t0]`
and `[t1,T1]`. If there are no intervals, the inverse becomes the entire validity interval.

Reference/API
-------------
.. automodule:: satmad.utils.timeinterval
    :members:
    :undoc-members:
