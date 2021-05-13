# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Discrete time events and intervals finding methods.

"""
import numpy as np
from astropy.timeseries import TimeSeries

from satmad.utils.interpolators import DiscreteTimeData
from satmad.utils.timeinterval import TimeInterval, TimeIntervalList


class DiscreteTimeEvents:
    """
    Discrete time events and intervals finding and storage class.

    An *time event* is described as a function value crossing a certain threshold in
    increasing (negative to positive) and decreasing (positive to negative) directions,
    defining a time interval.

    The class initialises interpolators to find where the `value_list` crosses the
    `crossing_value`, which is essentially a root finding problem.

    Parameters
    ----------
    time_list : Time
        List of time instances
    value_list : ndarray or list
        List of values corresponding to time instances (should be of the same size as
        `time_list`)
    crossing_value : float
        Crossing value for the `value_list`
    neg_to_pos_is_start : bool
        If `True` value turning from negative to positive marks a *start* event,
        otherwise marks an *end* event

    Raises
    ------
    ValueError
        `time_list` and `value_list` are of different sizes
    """

    def __init__(
        self, time_list, value_list, crossing_value=0, neg_to_pos_is_start=True
    ):

        # Init underlying Discrete Data object
        if isinstance(value_list, np.ndarray):
            self.discrete_data = DiscreteTimeData(
                time_list, value_list - crossing_value
            )
        else:
            self.discrete_data = DiscreteTimeData(
                time_list, np.array(value_list) - crossing_value
            )

        # Find start / end intervals
        self.start_end_intervals = self._find_intervals(time_list, neg_to_pos_is_start)

        # Find max / min events
        self.max_min_table = self._find_extrema_events()

        # Return the table to absolute values
        self.max_min_table["value"] = self.max_min_table["value"] + crossing_value

    def _find_intervals(self, time_list, neg_to_pos_is_start):
        """
        Finds the start / end intervals within the time list.

        Parameters
        ----------
        time_list : Time
            List of time instances
        neg_to_pos_is_start : bool
            If `True` value turning from negative to positive marks a *start* event,
            otherwise marks an *end* event

        Returns
        -------
        intervals : TimeIntervalList
            Time intervals marked with start and end events
        """

        init_time = time_list[0]
        end_time = time_list[-1]

        # *** Find start / end events (do not extrapolate to find roots) ***
        # No roots mean no events during this time
        start_end_events = self.discrete_data.roots(extrapolate=False)

        intervals = []
        if start_end_events.size != 0:
            # there should be some roots within the search range to check for intervals

            # classify start / end events and fill the timetable
            events = self._classify_start_end_events(
                start_end_events, neg_to_pos_is_start
            )

            # Put events into intervals
            event_start = init_time

            for event in events:
                # ignore events outside the search interval
                if self.search_interval.is_in_interval(event["time"]):
                    if event["type"] == "start":
                        event_start = event["time"]
                    else:
                        intervals.append(TimeInterval(event_start, event["time"]))

            # check for the last event - force end time if it is out of search bounds
            if events[-1]["type"] == "end":
                if self.search_interval.is_in_interval(events[-1]["time"]):
                    intervals.append(TimeInterval(event_start, events[-1]["time"]))
                else:
                    intervals.append(
                        TimeInterval(event_start, self.search_interval.end)
                    )
            else:
                intervals.append(TimeInterval(event_start, self.search_interval.end))

        # If interval list is empty by this stage, either the complete duration is
        # is a valid interval with no event (e.g. a continuously visible satellite)
        # or the complete duration is an invalid interval with no event (e.g. a
        # satellite continuously outside visibility)
        if not intervals:
            mid_time = self.search_interval.start + 0.5 * self.search_interval.duration
            value = self.discrete_data.interpolate(mid_time)

            if neg_to_pos_is_start:
                if value > 0:
                    # event already started
                    intervals.append(
                        TimeInterval(
                            self.search_interval.start, self.search_interval.end
                        )
                    )
            else:
                if value < 0:
                    # event already started
                    intervals.append(
                        TimeInterval(
                            self.search_interval.start, self.search_interval.end
                        )
                    )

        # generate the list of time intervals
        return TimeIntervalList(intervals, start_valid=init_time, end_valid=end_time)

    def _classify_start_end_events(self, start_end_events, neg_to_pos_is_start):
        """
        Classify the list of events into a timetable of start and end events.

        Parameters
        ----------
        start_end_events : ndarray
            Array of start and end event times in days, starting from `init_time`
        neg_to_pos_is_start : bool
            If `True` value turning from negative to positive marks a *start* event,
            otherwise marks an *end* event

        Returns
        -------
        events_table : TimeSeries
            Timetable of start and end events
        """

        # loop through each root and classify them as start / end events
        events_list = []
        for deriv in self.discrete_data.deriv_interpolate(start_end_events):

            if neg_to_pos_is_start:
                if deriv > 0:
                    events_list.append("start")
                else:
                    events_list.append("end")
            else:
                if deriv < 0:
                    events_list.append("start")
                else:
                    events_list.append("end")

        # init events table
        events_table = TimeSeries(time=start_end_events, data={"type": events_list})

        return events_table

    def _find_extrema_events(self):
        """
        Finds the extrema (max/min) events, contained within the time intervals.

        Returns
        -------
        max_min_table : TimeSeries
            Timetable of max and min events

        """

        # *** Find max / min events ***
        min_max_events = self.discrete_data.deriv_roots()

        # loop through each root and classify them as max / min events
        events_list = []
        value_list = []
        for value in self.discrete_data.interpolate(min_max_events):

            value_list.append(value)
            if value > 0:
                events_list.append("max")
            else:
                events_list.append("min")

        # init events table
        max_min_table = TimeSeries(
            time=min_max_events,
            data={"type": events_list, "value": value_list},
        )

        # Filter the events to those within the intervals
        i = 0
        remove_indexes = []
        for event_row in max_min_table:
            event_within_interval = self.start_end_intervals.is_in_interval(
                event_row["time"]
            )
            if not event_within_interval:
                remove_indexes.append(i)
            i = i + 1

        # Remove the ones outside the intervals
        max_min_table.remove_rows(remove_indexes)

        return max_min_table

    @property
    def search_interval(self):
        """Search interval for the events."""
        return self.discrete_data.data_interval
