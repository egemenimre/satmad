# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Two-Line-Elements to represent satellites in Earth orbit.

"""
from astropy.time import Time
from sgp4.exporter import export_tle
from sgp4.model import WGS72, Satrec


class TLE:
    """
    A two-line element set (TLE) is a data format encoding a list of TEME
    (True Equator, Mean Equinox) mean orbital elements
    of an Earth-orbiting object for a given point in time, called the Epoch Time.

    These orbital elements are solely for use with the SGP4 propagator due to the
    analytical orbit theory used in its derivation.

    See the `TLE page in Wikipedia
    <https://en.wikipedia.org/wiki/Two-line_element_set>`_ or `NASA definition
    <https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/SSOP_Help/tle_def.html>`_
    for more information.

    An example TLE is given as::

        ISS (ZARYA)
        1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
        2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537

    A TLE object is usually initialised from these two lines of strings:

        >>> line1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
        >>> line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"
        >>> tle = TLE.from_tle(line1, line2, "ISS (ZARYA)")

   """

    # # Gravitational constants(defaults to WGS72)
    # _mu = 398600.8 * u.km ** 3 / u.s ** 2  # in km3 / s2
    # _earth_radius = 6378.135 * u.km  # km
    # _j2 = 0.001082616

    # Hardcoded defaults
    _grav_model = WGS72
    _opsmode = "i"

    def __init__(
        self,
        epoch: Time,
        inclination,
        raan,
        eccentricity,
        arg_perigee,
        mean_anomaly,
        mean_motion,
        bstar,
        n_dot,
        n_dotdot=0.0,
        name="NONAME",
        intl_designator="12345A",
        sat_num=99999,
        classification="U",
        rev_nr=0,
        el_number=1,
    ):
        """
        Parameters
        ----------
        epoch : Time
            Epoch Time corresponding to the orbital elements (nominally very near the
            time of true ascending node passage)
        inclination : float or Quantity
            TEME mean inclination of the orbit [rad]
        raan : float or Quantity
            TEME mean right ascension of ascending node (RAAN) of the orbit [rad]
        eccentricity : float
            mean eccentricity of the orbit
        arg_perigee : float or Quantity
            TEME mean argument of perigee [rad]
        mean_anomaly : float or Quantity
            mean anomaly of the orbit [rad]
        mean_motion : float or Quantity
            mean motion of the orbit [orbits/day]
        bstar : float or Quantity
            sgp4 type drag coefficient [1 / earth radius] (see TLE class documentation)
        n_dot : float or Quantity
            First time derivative of the mean motion or Ballistic Coefficient
            [revs/day] (see TLE class documentation)
        n_dotdot : float or Quantity
            Second time derivative of the mean motion (see TLE class documentation)
        name : str
            Common name of the satellite
        intl_designator : str
            international designator on card 1 (up to 8 characters) (see class
            definition)
        sat_num : int
            satellite catalog number (see TLE class documentation)
        classification : str
            Classification (`U`=Unclassified, `C`=Classified, `S`=Secret)
        rev_nr : int
            Revolution number of the object at Epoch Time [revolutions]
        el_number : int
            Element set number. Incremented when a new TLE is generated for this object.
        """
        # init internal satrec object
        self._satrec = Satrec()

        # recreate the epoch composite
        yydd_str = epoch.utc.to_value("yday", subfmt="date").split(":")
        self._satrec.epochyr = int(yydd_str[0][-2:])
        self._satrec.epochdays = int(yydd_str[1]) + epoch.mjd % 1
        epoch_yydd = self._satrec.epochyr * 1000 + self._satrec.epochdays

        # fill the Satrec object
        self._satrec.sgp4init(
            TLE._grav_model,
            TLE._opsmode,
            sat_num,
            epoch_yydd,
            bstar,
            n_dot,
            n_dotdot,
            eccentricity,
            arg_perigee,
            inclination,
            mean_anomaly,
            mean_motion,
            raan,
        )

        # fill time with precise Time information
        self._satrec.jdsatepoch = epoch.utc.jd1
        self._satrec.jdsatepochF = epoch.utc.jd2

        # load other TLE info
        self._satrec.classification = classification
        self._satrec.intldesg = intl_designator
        self._satrec.revnum = rev_nr
        self._satrec.ephtype = "0"  # always 0 in distributed ephemeris
        self._satrec.elnum = el_number

        # load TLE name
        self._name = name

    @classmethod
    def from_tle(
        cls, line1, line2, name="No Name",
    ):
        """
        Initialises the TLE from two strings.

        Parameters
        ----------
        line1 : str
            First line of the TLE
        line2 : str
            Second line of the TLE
        name : str
            Name of the object / satellite

        Returns
        -------
        TLE
            `TLE` object initialised with the two line input.
        """
        # create object without calling `__init__`
        tle = cls.__new__(cls)

        # load TLE name
        tle._name = name

        # init Satrec object with TLE strings
        tle._satrec = Satrec.twoline2rv(line1, line2, whichconst=TLE._grav_model)

        return tle

    def __str__(self):
        """
        Exports the TLE as a Three-Line string, with name, line1 and line2.

        A sample output looks like this::

            ISS (ZARYA)
            1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
            2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537

        Returns
        -------
        str
            string containing line 1, line 2 and name of the satellite
        """
        line1, line2 = export_tle(self._satrec)

        txt = ""
        if self._name is not None and self._name.strip():
            txt += self._name + "\n"

        txt += line1 + "\n" + line2 + "\n"

        return txt

    # def _set_default_info(self):
    #     self._satrec.intldesg = "12345A"
    #     self._satrec.sat_num = 99999
    #     self._satrec.classification = "U"
