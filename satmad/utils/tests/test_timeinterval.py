# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Test TimeInterval class and associated methods and functionalities.

"""
import numpy as np
import pytest
from astropy.time import Time, TimeDelta

from satmad.utils.timeinterval import TimeInterval, TimeIntervalList

before = TimeInterval(
    Time("2020-04-09T00:00:00", scale="utc"), Time("2020-04-11T00:00:00", scale="utc"),
)
within = TimeInterval(
    Time("2020-04-11T00:05:00", scale="utc"), Time("2020-04-11T00:08:00", scale="utc"),
)
intersect = TimeInterval(
    Time("2020-04-10T00:00:00", scale="utc"), Time("2020-04-11T00:08:00", scale="utc"),
)
exact = TimeInterval(
    Time("2020-04-11T00:00:00", scale="utc"), Time("2020-04-11T00:10:00", scale="utc"),
)
after = TimeInterval(
    Time("2020-04-11T00:10:00", scale="utc"), Time("2020-04-12T00:00:00", scale="utc"),
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


def test_init_with_durations(init_times, durations):
    """Test initialisation with durations."""

    intervals = TimeIntervalList(init_times, durations)

    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n"
        "[ 2020-04-12T00:00:00.000  2020-04-12T00:20:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]\n"
    )

    assert truth_txt == str(intervals)


def test_init_with_end_times(init_times, durations):
    """Test initialisation with explicit end times."""

    end_times = init_times + durations

    intervals = TimeIntervalList(init_times, end_times)

    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n"
        "[ 2020-04-12T00:00:00.000  2020-04-12T00:20:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]\n"
    )

    assert truth_txt == str(intervals)


def test_is_intersecting(init_times, durations):
    """Test `is_intersecting` method."""
    interval = TimeIntervalList(init_times, durations).get_interval(0)

    assert interval.is_intersecting(before) is False
    assert interval.is_intersecting(within) is True
    assert interval.is_intersecting(intersect) is True
    assert interval.is_intersecting(exact) is True
    assert interval.is_intersecting(after) is False


def test_contains(init_times, durations):
    """Test `contains` method."""
    interval = TimeIntervalList(init_times, durations).get_interval(0)

    assert interval.contains(before) is False
    assert interval.contains(within) is True
    assert interval.contains(exact) is True
    assert interval.contains(intersect) is False
    assert interval.contains(after) is False


def test_is_within_interval(init_times, durations):
    """Test `is_within_interval` method."""
    interval = TimeIntervalList(init_times, durations).get_interval(0)

    t_before = Time("2020-04-09T00:00:00", scale="utc")
    t_eqinit = Time("2020-04-11T00:00:00", scale="utc")
    t_within = Time("2020-04-11T00:05:00", scale="utc")
    t_eqend = Time("2020-04-11T00:10:00", scale="utc")
    t_after = Time("2020-04-11T12:00:00", scale="utc")

    assert interval.is_within_interval(t_before) is False
    assert interval.is_within_interval(t_eqinit) is True
    assert interval.is_within_interval(t_within) is True
    assert interval.is_within_interval(t_eqend) is True
    assert interval.is_within_interval(t_after) is False


def test_intersect(init_times, durations):
    """Test `intersect` method."""
    interval = TimeIntervalList(init_times, durations).get_interval(0)

    assert interval.intersect(before) is None
    assert interval.intersect(within).is_equal(within)
    assert interval.intersect(intersect).is_equal(
        TimeInterval(interval.start, intersect.end)
    )
    assert exact.is_equal(interval.intersect(exact))
    assert interval.intersect(after) is None


def test_union(init_times, durations):
    """Test `union` method."""
    interval = TimeIntervalList(init_times, durations).get_interval(0)

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
        assert within.expand(start_delta=TimeDelta(-3 * 60, format="sec"))


def test_expand_shrink_negative():
    """Test `expand` method with shrink to negative - should raise `ValueError`."""
    with pytest.raises(ValueError):
        # Test a negative duration shrinkage - this raises a ValueError
        assert within.expand(
            start_delta=TimeDelta(-5 * 60, format="sec"),
            end_delta=TimeDelta(-5 * 60, format="sec"),
        )
