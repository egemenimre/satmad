# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
TLE storage classes.

"""
from satmad.propagation.tle import TLE


class TleStorage:
    """TLE storage class that keeps a list of TLE data from multiple satellites,
    at multiple times ant without any ordering.

    This class is the entry point for reading a TLE file, from which various sublists
    (e.g. single satellite, all LEO sats etc.) can be derived."""

    @classmethod
    def from_path(cls, tle_file_path):
        """
        Read a set of TLE data from file. Tries to extract satellite names from the
        list, if no name is found, an empty string (not `None`) is assigned.

        Parameters
        ----------
        tle_file_path : Path
            Path of the text file to be read

        Returns
        -------
        tle_storage
            A `TleStorage` object that contains the list of TLE data
        """
        with open(tle_file_path, "r") as f:
            tle_source_str = f.read()

        return cls.from_string(tle_source_str)

    @classmethod
    def from_string(cls, tle_string):
        """
        Read a set of TLE data from a text. Tries to extract satellite names from the
        list, if no name is found, an empty string (not `None`) is assigned.

        Parameters
        ----------
        tle_string : str
            String containing successive TLE data

        Returns
        -------
        tle_storage
            A `TleStorage` object that contains the list of TLE data
        """

        tle_source_str = tle_string.split("\n")

        # create object without calling `__init__`
        tle_storage = cls.__new__(cls)

        # parse the string and load parsed items into list
        tle_storage.tle_list = _parse_tle_list(tle_source_str)

        return tle_storage


def _parse_tle_list(tle_source_str_list):
    """
    Parses the TLE list.

    Parameters
    ----------
    tle_source_str_list : list[str]
        TLE data as a list of strings.

    Returns
    -------
    tle_list : list[TLE]
        List of TLE data
    """

    tle_list = []

    name = line1 = line2 = None
    for i, line in enumerate(tle_source_str_list):

        # strip spaces and EOF around the line
        line = line.strip()
        # skip empty lines
        if not line.strip():
            continue

        if __is_tle_line(line, 1):
            line1 = line
            if __is_tle_line(tle_source_str_list[i + 1], 2):
                line2 = tle_source_str_list[i + 1]
                if i > 0 and (
                    not __is_tle_line(tle_source_str_list[i - 1], 1)
                    and not __is_tle_line(tle_source_str_list[i - 1], 2)
                ):
                    name = tle_source_str_list[i - 1].strip("\n ")

        if line1 and line2:
            tle = TLE.from_tle(line1, line2, name if name else "")
            tle_list.append(tle)
            # reset temp fields
            name = line1 = line2 = None

    return tle_list


def __is_tle_line(line, line_nr):
    """Checks whether line is of type Line 1 or Line 2 standard TLE line."""
    if line.strip().startswith(f"{line_nr} "):
        return True
    else:
        return False
