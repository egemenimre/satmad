Finding Events and Intervals in Discrete Time Data
===================================================

Introduction
------------
A common family of problems in satellite mission analysis involve finding *events* and *intervals* such as
finding the communication times of a satellite with a groundstation, satellite eclipse entrance and exit
times or finding the exact time a satellite crosses the Equator. They require the precise calculation
of the time where a calculated value (or its derivative) crosses a certain threshold. Usually the data
is not defined as a continuous and differentiable function, but only as a discrete data set corresponding
to the variation of a value with respect to time.

To be able to solve this problem, the :mod:`.discrete_time_intervals` module provides the
the :class:`.DiscreteTimeEvents` class. This class generates an interpolator for the data set and
computes the roots to find the start / end events that mark intervals within this data set. Similarly,
it computes the derivative of this interpolator and computes the roots again to find the max / min events
within this data set.

Using the :class:`.DiscreteTimeEvents` Class
----------------------------------------------
The :class:`.DiscreteTimeEvents` class is designed as a helper to interpret the values representing
the behaviour of a parameter. A simple use case would be to "for a given ground location, compute
the sunrise and sunset times as well as the time of the highest point of the Sun in the sky".
If the location has mountains around it, limiting the view of the Sun, a minimum elevation
can also be defined. The altitude or elevation values of the Sun as seen from a ground location would be
computed elsewhere and fed into this class.

Therefore, in this problem formulation, the `value_list` would be the Sun altitude values in degrees
over several days and the `time_list` would be the times corresponding to these altitude values. `crossing_value`
would be zero for horizon or perhaps 5, to illustrate the minimum altitude. Finally, as sunrise is defined as
the altitude crossing from negative to positive, the `neg_to_pos_is_start` parameter would be set to `True`.
If the night duration was sought, then this would be set to `False`, as night is defined as the time where the
altitude of the Sun changes from positive to negative (except things are a little bit more complicated than that,
`see here for more information <https://www.timeanddate.com/astronomy/different-types-twilight.html>`_).

To retrieve the sunrise / sunset intervals, the `start_end_intervals` can be queried, yielding a list of intervals.
Similarly, to find the times where the Sun is at the highest point, the `max_min_table` can be queried. Note that,
by definition, only the max / min events inside the intervals are shown and the remaining events are filtered.
Therefore, the max / min table will only have the Sun highest points where it is above horizon and not the lowest
points where it is well below the horizon and is technically midnight.

The class also supports events that has no definite start or end. For example, if the Sun was already up at the
first element of the `time_list`, then the first interval starts at that time. Similarly, if the Sun has not set
at the final element of the `time_list`, then the last interval is closed with this last element. If the Sun was
up all through the `time_list`, then there will be a single interval from beginning to the end.
Conversely, if the Sun was not up at all throughout the data points, then there will be no intervals.

Interpolation, Time Step and Accuracy
-------------------------------------
The interpolator :class:`scipy.interpolate.CubicSpline` is used under the hood to interpolate the values.
The accuracy of the start / end  as well as max / min events is defined by the time step of the data,
depending on the problem. For example, the sunrise / sunset problem above is very slowly varying and
would not require a time step of 60 seconds. On the other hand, finding the communication times of
a LEO satellite accurately would potentially require a stepsize of 20-40 seconds. Therefore, it is highly
recommended to check the convergence with smaller time steps and decide on a time step for
the specific problem and required level of accuracy.

Coarse data with large time steps can still yield reasonably accurate results if the slope around the
"root" is steep (i.e., transition from negative to positive values happening very quickly and clearly).
On the other hand, slowly varying events and coarse data do not play well, as it becomes extremely difficult
to find the exact point where the maximum point occurs. For example, for an observer on the ground,
the Sun moves very slowly around its highest point, with its derivative changing from positive to negative.
With a coarse data set, it is very challenging for a numerical algorithm to find the exact time with
a high accuracy.

Finally, best results can be had with continuous and smooth data. Heavy non-linearities require very small time
steps and there is no guarantee that discontiunities will be handled well.

As with all numerical algorithms, you should *"know thy problem"*.


Reference/API
-------------
.. automodule:: satmad.utils.discrete_time_events
    :members:
    :undoc-members:
