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

from satmad.utils.timeinterval import (
    TimeInterval,
    _contains,
    _expand,
    _intersect,
    _interval_equals,
    _is_intersecting,
    _is_within_interval,
    _union,
)

before = (
    Time("2020-04-09T00:00:00", scale="utc"),
    Time("2020-04-11T00:00:00", scale="utc"),
)
within = (
    Time("2020-04-11T00:05:00", scale="utc"),
    Time("2020-04-11T00:08:00", scale="utc"),
)
intersect = (
    Time("2020-04-10T00:00:00", scale="utc"),
    Time("2020-04-11T00:08:00", scale="utc"),
)
exact = (
    Time("2020-04-11T00:00:00", scale="utc"),
    Time("2020-04-11T00:10:00", scale="utc"),
)
after = (
    Time("2020-04-11T00:10:00", scale="utc"),
    Time("2020-04-12T00:00:00", scale="utc"),
)


@pytest.fixture
def init_times():
    return Time("2020-04-10T00:00:00", scale="utc") + np.arange(1, 4) * TimeDelta(
        1, format="jd"
    )


@pytest.fixture
def durations():
    return np.arange(1, 4) * TimeDelta(600, format="sec")


def test_init_with_durations(init_times, durations):
    """Test initialisation with durations."""

    intervals = TimeInterval(init_times, durations)

    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n"
        "[ 2020-04-12T00:00:00.000  2020-04-12T00:20:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]\n"
    )

    assert truth_txt == str(intervals)


def test_init_with_end_times(init_times, durations):
    """Test initialisation with explicit end times."""

    end_times = init_times + durations

    intervals = TimeInterval(init_times, end_times)

    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n"
        "[ 2020-04-12T00:00:00.000  2020-04-12T00:20:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]\n"
    )

    assert truth_txt == str(intervals)


def test_is_intersecting(init_times, durations):
    """Test `is_intersecting` method."""
    interval = TimeInterval(init_times, durations).get_interval(0)

    assert _is_intersecting(before, interval) is False
    assert _is_intersecting(within, interval) is True
    assert _is_intersecting(intersect, interval) is True
    assert _is_intersecting(exact, interval) is True
    assert _is_intersecting(after, interval) is False


def test_contains(init_times, durations):
    """Test `_contains` method."""
    interval = TimeInterval(init_times, durations).get_interval(0)

    assert _contains(before, interval) is False
    assert _contains(within, interval) is True
    assert _contains(exact, interval) is True
    assert _contains(intersect, interval) is False
    assert _contains(after, interval) is False


def test_is_within_interval(init_times, durations):
    """Test `is_within_interval` method."""
    interval = TimeInterval(init_times, durations).get_interval(0)

    t_before = Time("2020-04-09T00:00:00", scale="utc")
    t_eqinit = Time("2020-04-11T00:00:00", scale="utc")
    t_within = Time("2020-04-11T00:05:00", scale="utc")
    t_eqend = Time("2020-04-11T00:10:00", scale="utc")
    t_after = Time("2020-04-11T12:00:00", scale="utc")

    assert _is_within_interval(t_before, interval) is False
    assert _is_within_interval(t_eqinit, interval) is True
    assert _is_within_interval(t_within, interval) is True
    assert _is_within_interval(t_eqend, interval) is True
    assert _is_within_interval(t_after, interval) is False


def test_intersect(init_times, durations):
    """Test `_intersect` method."""
    interval = TimeInterval(init_times, durations).get_interval(0)

    assert _intersect(before, interval) is None
    assert _interval_equals(_intersect(within, interval), within)
    assert _interval_equals(
        _intersect(intersect, interval), (interval[0], intersect[1])
    )
    assert _interval_equals(_intersect(exact, interval), exact)
    assert _intersect(after, interval) is None


def test_union(init_times, durations):
    """Test `_union` method."""
    interval = TimeInterval(init_times, durations).get_interval(0)

    assert _union(before, interval) is None
    assert _interval_equals(_union(within, interval), interval)
    assert _interval_equals(_union(intersect, interval), (intersect[0], interval[1]))
    assert _interval_equals(_union(exact, interval), interval)
    assert _union(after, interval) is None


def test_expand():
    """Test `_expand` method."""
    with pytest.raises(ValueError):
        # Test expansion
        expanded = _expand(
            within,
            delta_begin=TimeDelta(5 * 60, format="sec"),
            delta_end=TimeDelta(2 * 60, format="sec"),
        )
        assert _interval_equals(expanded, exact)

        # Test shrinkage
        shrunk = _expand(
            exact,
            delta_begin=TimeDelta(-5 * 60, format="sec"),
            delta_end=TimeDelta(-2 * 60, format="sec"),
        )
        assert _interval_equals(shrunk, within)

        # Test shrink to zero - this raises a ValueError
        assert _expand(within, delta_begin=TimeDelta(-3 * 60, format="sec"))

        # Test a negative duration shrinkage - this raises a ValueError
        assert _expand(
            within,
            delta_begin=TimeDelta(-5 * 60, format="sec"),
            delta_end=TimeDelta(-5 * 60, format="sec"),
        )
