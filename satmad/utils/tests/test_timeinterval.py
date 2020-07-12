# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Test `TimeInterval` class and associated methods and functionalities.

"""
import numpy as np
import pytest
from astropy import units as u
from astropy.time import Time, TimeDelta

from satmad.utils.timeinterval import _EPS_TIME, TimeInterval, TimeIntervalList

before = TimeInterval(
    Time("2020-04-09T00:00:00", scale="utc"),
    Time("2020-04-11T00:00:00", scale="utc"),
    end_inclusive=False,
)
within = TimeInterval(
    Time("2020-04-11T00:05:00", scale="utc"), Time("2020-04-11T00:08:00", scale="utc")
)
intersect = TimeInterval(
    Time("2020-04-10T00:00:00", scale="utc"), Time("2020-04-11T00:08:00", scale="utc")
)
exact = TimeInterval(
    Time("2020-04-11T00:00:00", scale="utc"), Time("2020-04-11T00:10:00", scale="utc")
)
after = TimeInterval(
    Time("2020-04-11T00:10:00", scale="utc"), Time("2020-04-12T00:00:00", scale="utc")
)


@pytest.fixture
def init_times():
    """Generates the initial times"""
    return Time("2020-04-10T00:00:00", scale="utc") + np.arange(1, 4) * TimeDelta(
        1, format="jd"
    )


@pytest.fixture
def durations():
    """Generates the initial durations for the initial times"""
    return np.arange(1, 4) * TimeDelta(600, format="sec")


def test_interval_init_with_copy(init_times, durations):
    """Test initialisation with copy."""
    interval_1 = TimeInterval(init_times, durations)
    interval_2 = TimeInterval(interval_1, None)
    interval_3 = TimeInterval(interval_1, None, replicate=True)

    assert interval_1.is_equal(interval_2)
    assert interval_1.is_equal(interval_3)


def test_interval_init_with_durations(init_times, durations):
    """Test initialisation with durations."""
    interval = TimeInterval(init_times, durations)

    truth_txt = "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]"
    assert truth_txt == str(interval)


def test_interval_init_with_end_times(init_times, durations):
    """Test initialisation with explicit end times."""
    end_times = init_times + durations

    interval = TimeInterval(init_times, end_times, replicate=True, end_inclusive=False)

    truth_txt = "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 )"
    assert truth_txt == str(interval)


def test_interval_init_switched_err(init_times, durations):
    """Test `init` with switched init and end times - should raise `ValueError`."""
    with pytest.raises(ValueError):
        end_times = init_times + durations
        TimeInterval(end_times, init_times)


def test_interval_init_zero_dur_err(init_times):
    """Test `init` with equal init and end times - should raise `ValueError`."""
    with pytest.raises(ValueError):
        TimeInterval(init_times, init_times + 0.1 * _EPS_TIME)


def test_interval_list_init_with_durations(init_times, durations):
    """Test initialisation with durations."""

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], durations[i]))

    intervals = TimeIntervalList(init_intervals)

    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n"
        "[ 2020-04-12T00:00:00.000  2020-04-12T00:20:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]\n"
    )

    assert truth_txt == str(intervals)


def test_interval_list_init_with_end_times(init_times, durations):
    """Test initialisation with explicit end times."""

    end_times = init_times + durations

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], end_times[i]))

    intervals = TimeIntervalList(
        init_intervals, start_valid=init_times[0], end_valid=end_times[-1]
    )

    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n"
        "[ 2020-04-12T00:00:00.000  2020-04-12T00:20:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]\n"
    )

    truth_validity = "[ 2020-04-11T00:00:00.000  2020-04-13T00:30:00.000 ]"

    assert truth_txt == str(intervals)
    assert truth_validity == str(intervals.valid_interval)


def test_interval_list_intersection(init_times, durations):
    """Test interval intersection."""

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], durations[i]))

    intervals = TimeIntervalList(init_intervals)

    interval_full_intersect = TimeInterval(Time("2020-04-13T00:00:00.000"), 5 * u.min)
    interval_part_intersect = TimeInterval(Time("2020-04-10T23:50:00.000"), 60 * u.min)
    interval_no_intersect = TimeInterval(Time("2020-04-11T12:00:00.000"), 5 * u.min)

    # *** Test TimeIntervalList vs. TimeInterval intersection ***
    assert intervals.is_intersecting(interval_full_intersect) is True
    assert "[ 2020-04-13T00:00:00.000  2020-04-13T00:05:00.000 ]" == str(
        intervals.intersect(interval_full_intersect)
    )

    assert intervals.is_intersecting(interval_part_intersect) is True
    assert "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]" == str(
        intervals.intersect(interval_part_intersect)
    )

    assert intervals.is_intersecting(interval_no_intersect) is False
    assert intervals.intersect(interval_no_intersect) is None

    # *** Test TimeIntervalList vs. TimeIntervalList intersection ***

    test_intervals = TimeIntervalList(
        [interval_full_intersect, interval_part_intersect, interval_no_intersect]
    )

    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:05:00.000 ]\n"
    )

    assert truth_txt == str(intervals.intersect_list(test_intervals))


def test_interval_list_union(init_times, durations):
    """Test interval union."""

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], durations[i]))

    intervals = TimeIntervalList(
        init_intervals,
        start_valid=Time("2020-04-11T00:00:00.000"),
        end_valid=Time("2020-04-14T00:00:00.000"),
    )

    interval_part_intersect = TimeInterval(Time("2020-04-10T23:50:00.000"), 60 * u.min)
    interval_no_intersect = TimeInterval(Time("2020-04-11T12:00:00.000"), 5 * u.min)
    interval_full_intersect = TimeInterval(Time("2020-04-13T00:00:00.000"), 5 * u.min)
    interval_outside = TimeInterval(Time("2020-04-13T23:50:00.000"), 15 * u.min)

    # *** Test TimeIntervalList vs. TimeInterval intersection ***
    assert "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]" == str(
        intervals.union(interval_full_intersect)
    )

    assert "[ 2020-04-10T23:50:00.000  2020-04-11T00:50:00.000 ]" == str(
        intervals.union(interval_part_intersect)
    )

    assert intervals.union(interval_no_intersect) is None

    # *** Test TimeIntervalList vs. TimeIntervalList intersection ***

    test_intervals = TimeIntervalList(
        [
            interval_full_intersect,
            interval_part_intersect,
            interval_no_intersect,
            interval_outside,
        ]
    )

    # print(intervals.union_list(test_intervals))
    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:50:00.000 ]\n"
        "[ 2020-04-11T12:00:00.000  2020-04-11T12:05:00.000 ]\n"
        "[ 2020-04-12T00:00:00.000  2020-04-12T00:20:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]\n"
        "[ 2020-04-13T23:50:00.000  2020-04-14T00:00:00.000 ]\n"
    )

    assert truth_txt == str(intervals.union_list(test_intervals))


def test_is_intersecting(init_times, durations):
    """Test `is_intersecting` method."""

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], durations[i]))

    interval = TimeIntervalList(init_intervals).get_interval(0)

    assert interval.is_intersecting(before) is False
    assert interval.is_intersecting(within) is True
    assert interval.is_intersecting(intersect) is True
    assert interval.is_intersecting(exact) is True
    assert interval.is_intersecting(after) is False


def test_contains(init_times, durations):
    """Test `contains` method."""

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], durations[i]))

    interval = TimeIntervalList(init_intervals).get_interval(0)

    assert interval.contains(before) is False
    assert interval.contains(within) is True
    assert interval.contains(exact) is True
    assert interval.contains(intersect) is False
    assert interval.contains(after) is False


def test_is_in_interval(init_times, durations):
    """Test `is_in_interval` method."""

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], durations[i]))

    interval = TimeIntervalList(init_intervals).get_interval(0)

    t_before = Time("2020-04-09T00:00:00", scale="utc")
    t_eqinit = Time("2020-04-11T00:00:00", scale="utc")
    t_within = Time("2020-04-11T00:05:00", scale="utc")
    t_eqend = Time("2020-04-11T00:10:00", scale="utc")
    t_after = Time("2020-04-11T12:00:00", scale="utc")

    assert interval.is_in_interval(t_before) is False
    assert interval.is_in_interval(t_eqinit) is True
    assert interval.is_in_interval(t_within) is True
    assert interval.is_in_interval(t_eqend) is True
    assert interval.is_in_interval(t_after) is False


def test_intersect(init_times, durations):
    """Test `intersect` method."""

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], durations[i]))

    interval = TimeIntervalList(init_intervals).get_interval(0)

    assert interval.intersect(before) is None
    assert interval.intersect(within).is_equal(within)
    assert interval.intersect(intersect).is_equal(
        TimeInterval(interval.start, intersect.end)
    )
    assert exact.is_equal(interval.intersect(exact))
    assert exact.is_equal(intersect) is False
    assert interval.intersect(after) is None


def test_union(init_times, durations):
    """Test `union` method."""

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], durations[i]))

    interval = TimeIntervalList(init_intervals).get_interval(0)

    assert interval.union(before) is None
    assert interval.union(within).is_equal(interval)
    assert interval.union(intersect).is_equal(
        TimeInterval(intersect.start, interval.end)
    )
    assert interval.union(exact).is_equal(interval)
    assert interval.union(after) is None


def test_expand():
    """Test `expand` method."""
    # Test expansion
    expanded = within.expand(
        start_delta=TimeDelta(5 * 60, format="sec"),
        end_delta=TimeDelta(2 * 60, format="sec"),
    )
    assert expanded.is_equal(exact)

    # Test shrinkage
    shrunk = exact.expand(
        start_delta=TimeDelta(-5 * 60, format="sec"),
        end_delta=TimeDelta(-2 * 60, format="sec"),
    )
    assert shrunk.is_equal(within)


def test_expand_shrink_zero():
    """Test `expand` method with shrink to zero - should raise `ValueError`."""
    with pytest.raises(ValueError):
        # Test shrink to zero - this raises a ValueError
        within.expand(start_delta=TimeDelta(-3 * 60, format="sec"))


def test_expand_shrink_negative():
    """Test `expand` method with shrink to negative - should raise `ValueError`."""
    with pytest.raises(ValueError):
        # Test a negative duration shrinkage - this raises a ValueError
        within.expand(
            start_delta=TimeDelta(-5 * 60, format="sec"),
            end_delta=TimeDelta(-5 * 60, format="sec"),
        )


def test_invert(init_times, durations):
    """Test `invert` method."""

    end_times = init_times + durations

    init_intervals = []
    for i, init_time in enumerate(init_times):
        init_intervals.append(TimeInterval(init_times[i], end_times[i]))

    intervals = TimeIntervalList(
        init_intervals,
        start_valid=init_times[0] - TimeDelta(1.0, format="jd"),
        end_valid=end_times[-1] + TimeDelta(1.0, format="jd"),
    )

    inverted_intervals = intervals.invert()

    truth_txt = (
        "[ 2020-04-10T00:00:00.000  2020-04-11T00:00:00.000 ]\n"
        "[ 2020-04-11T00:10:00.000  2020-04-12T00:00:00.000 ]\n"
        "[ 2020-04-12T00:20:00.000  2020-04-13T00:00:00.000 ]\n"
        "[ 2020-04-13T00:30:00.000  2020-04-14T00:30:00.000 ]\n"
    )

    truth_validity = "[ 2020-04-10T00:00:00.000  2020-04-14T00:30:00.000 ]"

    assert truth_txt == str(inverted_intervals)
    assert truth_validity == str(inverted_intervals.valid_interval)
