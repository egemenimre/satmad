# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
TLE storage classes.

"""
from abc import ABC
from enum import Enum
from typing import List

from satmad.propagation.tle import TLE


class TleFilterParams(Enum):
    """TLE Filtering Parameters."""

    PERIOD = "period"
    SM_AXIS = "sm_axis"
    NODE_ROTATION_RATE = "node_rotation_rate"
    ARGP_ROTATION_RATE = "argp_rotation_rate"
    EPOCH = "epoch"
    BSTAR = "bstar"
    N_DOT = "n_dot"
    N_DOTDOT = "n_dotdot"
    INCLINATION = "inclination"
    RAAN = "raan"
    ECCENTRICITY = "eccentricity"
    ARG_PERIGEE = "arg_perigee"
    MEAN_ANOMALY = "mean_anomaly"
    MEAN_MOTION = "mean_motion"
    NAME = "name"
    SAT_NUMBER = "sat_number"
    INTL_DESIGNATOR = "intl_designator"
    CLASSIFICATION = "classification"
    REV_NR = "rev_nr"


class _TleList(ABC):
    """Abstract Base Class for TLE lists."""

    tle_list: List[TLE] = []

    def filter_by_value(self, param, value):
        """
        Filters the TLE list for equivalence to a given value.

        For example `param` can be equal to `TleFilterParams.SAT_NUMBER` and `value`
        can be equal to `46945`, which filters all TLEs with a sat number 46945.

        Note that this filter is not appropriate for float values such as eccentricity
        where equivalence is very brittle. For these applications, use
        `filter_by_range()` or `filter_by_func()` instead.

        This method returns a `TleStorage` object even if the filtering result is empty.
        In this case, `tle_list` parameter of the `TleStorage` object will be an empty
        list.

        Also note that, the returned TLE objects in the list are just shallow copies
        of the objects in the master TLE list. Any change to them will change the
        relevant item in the backing TLE list.

        Parameters
        ----------
        param : TleFilterParams
            Filter parameter (such as name or satellite number)
        value
            Value associated with the parameter

        Returns
        -------
        TleStorage
            A `TleStorage` object that contains the filtered list of TLE data
        """
        filtered_list = [
            tle for tle in self.tle_list if getattr(tle, param.value) == value
        ]

        return TleStorage(filtered_list)

    def filter_by_func(self, param, filter_func):
        """
        Filters the TLE list for compliance to a given filter function.

        For example `param` can be equal to `TleFilterParams.SM_AXIS` and
        `filter_func` can be a filtering function that tests this parameter and
        returns `True` or `False` accordingly. The following function filters for
        semi-major axis values above 7000 km. Note that units should be defined and
        compatible with the value to be compared against. For semimajor axis, distance
        units such as meters are acceptable but using no dimensions or using
        wrong units (such as degrees) will throw an error.

        >>> from astropy import units as u
        >>> def sma_filter(a):
        >>>     return True if a > 7000 * u.km else False

        For exact equivalences (such as satellite names or ID numbers),
        using `filter_by_value` method will be easier and more appropriate.

        This method returns a `TleStorage` object even if the filtering result is empty.
        In this case, `tle_list` parameter of the `TleStorage` object will be an empty
        list.

        Also note that, the returned TLE objects in the list are just shallow copies
        of the objects in the master TLE list. Any change to them will change the
        relevant item in the backing TLE list.

        Parameters
        ----------
        param : TleFilterParams
            Filter parameter (such as name or satellite number)
        filter_func
            Function to test the parameter against

        Returns
        -------
        TleStorage
            A `TleStorage` object that contains the filtered list of TLE data
        """
        filtered_list = [
            tle for tle in self.tle_list if filter_func(getattr(tle, param.value))
        ]

        return TleStorage(filtered_list)

    def filter_by_range(self, param, min_value=None, max_value=None):
        """
        Filters the TLE list for compliance to a given min/max values.

        The test is simply:

        `max_value > param > min_value`

        For example `param` can be equal to `TleFilterParams.SM_AXIS`, then
        this parameter is tested against the minimum and maximum semimajor axis values
        supplied. If `None` is supplied for `min_value` or `max_value`, then there
        is no range or range check defined for this parameter. For example,
        if `min_value` is `None`, the parameter check reduces to `max_value > param`.

        Note that units should be defined and
        compatible with the value to be compared against. For semimajor axis, distance
        units such as meters are acceptable but using no dimensions or using
        wrong units (such as degrees) will throw an error.

        For exact equivalences (such as satellite names or ID numbers),
        using `filter_by_value` method will be easier and more appropriate.

        This method returns a `TleStorage` object even if the filtering result is empty.
        In this case, `tle_list` parameter of the `TleStorage` object will be an empty
        list.

        Also note that, the returned TLE objects in the list are just shallow copies
        of the objects in the master TLE list. Any change to them will change the
        relevant item in the backing TLE list.

        Parameters
        ----------
        param : TleFilterParams
            Filter parameter (such as name or satellite number)
        min_value
            Minimum value to test the parameter against
        max_value
            Maximum value to test the parameter against

        Returns
        -------
        TleStorage
            A `TleStorage` object that contains the filtered list of TLE data
        """
        # `min_value` and `max_value` are quantities and should be checked explicitly
        # for `None`, otherwise can be interpreted as `True` or `False`.
        if min_value is not None and max_value is not None:
            filtered_list = [
                tle
                for tle in self.tle_list
                if max_value > getattr(tle, param.value) > min_value
            ]
        elif min_value is not None:
            filtered_list = [
                tle for tle in self.tle_list if getattr(tle, param.value) > min_value
            ]
        elif max_value is not None:
            filtered_list = [
                tle for tle in self.tle_list if max_value > getattr(tle, param.value)
            ]
        else:
            filtered_list = []

        return TleStorage(filtered_list)


class TleTimeSeries(_TleList):
    """TLE storage class that keeps a list of TLE data from a single satellite,
    at multiple times and with time order.

    The entry point is ideally the `TleStorage` class, where a TLE file is usually read,
    and a single satellite is filtered. Once this class is initialised, various sublists
    (e.g. specific time range) can be derived.

    Parameters
    ----------
    tle_list : list[TLE]
        initial list of TLE objects (shallow copied into object)
    """

    def __init__(self, tle_list, sat_number):

        # init a TLE Storage and filter for the sat number
        self.tle_list = (
            TleStorage(tle_list)
            .filter_by_value(TleFilterParams.SAT_NUMBER, sat_number)
            .tle_list
        )

        # order the internal TLE list with respect to epoch
        self.tle_list.sort(key=lambda tle: tle.epoch)


class TleStorage(_TleList):
    """TLE storage class that keeps a list of TLE data from multiple satellites,
    at multiple times and without any ordering.

    This class is the entry point for reading a TLE file, from which various sublists
    (e.g. single satellite, all LEO sats etc.) can be derived.

    Parameters
    ----------
    tle_list : list[TLE]
        initial list of TLE objects (shallow copied into object)
    """

    def __init__(self, tle_list):
        self.tle_list = tle_list

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
        TleStorage
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
        TleStorage
            A `TleStorage` object that contains the list of TLE data
        """

        tle_source_str = tle_string.split("\n")

        # create object without calling `__init__`
        tle_storage = cls.__new__(cls)

        # parse the string and load parsed items into list
        tle_storage.tle_list = _parse_tle_list(tle_source_str)

        return tle_storage

    def to_tle_timeseries(self, sat_number):
        """
        Filters the TLE list for a single satellite to initialise a `TleTimeSeries`.

        Parameters
        ----------
        sat_number
            Satellite Catalog Number

        Returns
        -------
        TleTimeSeries
            A `TleTimeSeries` object that contains the list of TLE data of a single
            satellite

        """
        return TleTimeSeries(self.tle_list, sat_number)


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
                    if name.startswith("0 "):
                        name = name[2:]

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
