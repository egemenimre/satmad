# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
TLE Storage class tests

"""
from pathlib import Path

from satmad.propagation.tle_storage import TleStorage

extra_path = Path("satmad", "propagation", "tests")


def test_parse_storage_file():
    """Test parsing of the storage file with mixed TLE input."""
    file_path = Path("data", "tle_mixed_1.txt")

    working_dir = Path.cwd()
    file_path = process_paths(working_dir, file_path)

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


def process_paths(working_dir, path):
    """
    Processes the path depending on the run environment.
    """
    file_path = working_dir.joinpath(path)
    if not working_dir.joinpath(file_path).exists():
        file_path = working_dir.joinpath(extra_path).joinpath(path)

    return file_path
