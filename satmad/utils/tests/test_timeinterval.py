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

from satmad.utils.timeinterval import TimeInterval


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

    # combo = Time(
    #     np.concatenate((init_times.isot, end_times.isot), axis=0)
    #     .reshape(2, 3)
    #     .transpose()
    # )

    intervals = TimeInterval(init_times, end_times)

    truth_txt = (
        "[ 2020-04-11T00:00:00.000  2020-04-11T00:10:00.000 ]\n"
        "[ 2020-04-12T00:00:00.000  2020-04-12T00:20:00.000 ]\n"
        "[ 2020-04-13T00:00:00.000  2020-04-13T00:30:00.000 ]\n"
    )

    assert truth_txt == str(intervals)
