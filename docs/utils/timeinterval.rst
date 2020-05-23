Time Interval
==============

Introduction
------------
The :class:`.TimeInterval` package provides a time interval functionality i.e.,
an time interval with a start and end time/date, using the high precision
:class:`astropy.time.Time` classes under the hood.

If the `Time` object passed as the initial time (`init_times`) has multiple
time instances, then the end time definitions (in :class:`astropy.time.Time` or
:class:`astropy.time.TimeDelta`) should have the same number of time instances
to match.

While the current functionality is rather limited, soon it is going to feature
operations such as Union and Intersection.

Reference/API
-------------
.. automodule:: satmad.utils.timeinterval
    :members:
    :undoc-members:
