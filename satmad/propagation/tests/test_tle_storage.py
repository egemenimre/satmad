# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
TLE Storage class tests

"""
from pathlib import Path

from astropy import units as u

from satmad.propagation.tle_storage import TleFilterParams, TleStorage

extra_path = Path("satmad", "propagation", "tests")

mixed_tle_file_path_1 = Path("data", "tle_mixed_1.txt")


def test_parse_storage_file():
    """Test parsing of the storage file with mixed TLE input."""
    file_path = process_paths(mixed_tle_file_path_1)

    tle_storage = TleStorage.from_path(file_path)

    truth_no_of_elements = 19
    truth_pos = 3
    truth_name = "APRIZESAT 1"
    truth_line1 = (
        "1 28372U 04025G   21088.15052337  .00000024  00000-0  18112-4 0  9991"
    )
    truth_line2 = (
        "2 28372  98.4324  42.3607 0046090 190.1598 169.8668 14.49665428885585"
    )

    # check number of elements
    assert len(tle_storage.tle_list) == truth_no_of_elements

    # check specific element at `truth_pos`
    assert str(tle_storage.tle_list[truth_pos]).split("\n")[0:3] == [
        truth_name,
        truth_line1,
        truth_line2,
    ]


def test_filter_value():
    """Tests filtering for a value equivalence."""
    file_path = process_paths(mixed_tle_file_path_1)

    tle_storage = TleStorage.from_path(file_path)

    filtered_list_sat_nr = tle_storage.filter_by_value(
        TleFilterParams.SAT_NUMBER, 46495
    )
    filtered_list_int_deg = tle_storage.filter_by_value(
        TleFilterParams.INTL_DESIGNATOR, "18014H"
    )
    empty_filtered_list = tle_storage.filter_by_value(TleFilterParams.SAT_NUMBER, 56495)

    assert tle_storage.tle_list[0].name == "APRIZESAT 2"
    assert len(filtered_list_sat_nr.tle_list) == 2
    assert len(filtered_list_int_deg.tle_list) == 1
    assert len(empty_filtered_list.tle_list) == 0


def test_filter_func():
    """Tests filtering for a value range with a function."""
    file_path = process_paths(mixed_tle_file_path_1)

    tle_storage = TleStorage.from_path(file_path)

    def sma_filter_1(a):
        """Semimajor axis filter min."""
        return True if a > 7000 * u.km else False

    def sma_filter_2(a):
        """Semimajor axis filter max."""
        return True if 7000 * u.km > a else False

    def sma_filter_3(a):
        """Semimajor axis filter min/max."""
        return True if 7100 * u.km > a > 7000 * u.km else False

    filtered_list_sma_1 = tle_storage.filter_by_func(
        TleFilterParams.SM_AXIS, sma_filter_1
    )
    filtered_list_sma_2 = tle_storage.filter_by_func(
        TleFilterParams.SM_AXIS, sma_filter_2
    )
    filtered_list_sma_3 = tle_storage.filter_by_func(
        TleFilterParams.SM_AXIS, sma_filter_3
    )

    assert len(filtered_list_sma_1.tle_list) == 11
    assert len(filtered_list_sma_2.tle_list) == 8
    assert len(filtered_list_sma_3.tle_list) == 8


def test_filter_range():
    """Tests filtering for a value range with min/max parameters."""
    file_path = process_paths(mixed_tle_file_path_1)

    tle_storage = TleStorage.from_path(file_path)

    filtered_list_sma_1 = tle_storage.filter_by_range(
        TleFilterParams.SM_AXIS, min_value=7000 * u.km
    )

    filtered_list_sma_2 = tle_storage.filter_by_range(
        TleFilterParams.SM_AXIS, max_value=7000 * u.km
    )

    filtered_list_sma_3 = tle_storage.filter_by_range(
        TleFilterParams.SM_AXIS, min_value=7000 * u.km, max_value=7100 * u.km
    )

    filtered_list_sma_4 = tle_storage.filter_by_range(TleFilterParams.SM_AXIS)

    assert len(filtered_list_sma_1.tle_list) == 11
    assert len(filtered_list_sma_2.tle_list) == 8
    assert len(filtered_list_sma_3.tle_list) == 8
    assert len(filtered_list_sma_4.tle_list) == 0


def process_paths(path):
    """
    Processes the path depending on the run environment.
    """
    working_dir = Path.cwd()

    file_path = working_dir.joinpath(path)
    if not working_dir.joinpath(file_path).exists():
        file_path = working_dir.joinpath(extra_path).joinpath(path)

    return file_path
