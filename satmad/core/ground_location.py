# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Ground and Earth Location classes. These are based on
:class:`astropy.coordinates.EarthLocation` class, but generalised
to non-Earth cases.

"""
import collections

import erfa
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, Latitude, Longitude

from satmad.core.celestial_bodies import EARTH_ELLIPSOID_GRS80

GeodeticLocation = collections.namedtuple("GeodeticLocation", ["lon", "lat", "height"])


class GroundLocation(u.Quantity):
    """
    Location on a Central Body (e.g. Earth or Moon).

    Initialization is first attempted assuming geocentric (x, y, z) coordinates
    are given; if that fails, another attempt is made assuming geodetic
    coordinates (longitude, latitude, height above a reference ellipsoid).
    When using the geodetic forms, Longitudes are measured increasing to the
    east, so west longitudes are negative. Internally, the coordinates are
    stored as geocentric.

    To ensure a specific type of coordinates is used, use the corresponding
    class methods (`from_geocentric` and `from_geodetic`) or initialize the
    arguments with names (``x``, ``y``, ``z`` for geocentric; ``lon``, ``lat``,
    ``height`` for geodetic).  See the class methods for details.


    Notes
    -----
    This class fits into the coordinates transformation framework in that it
    encodes a position on the `~astropy.coordinates.ITRS` frame.  To get a
    proper `~astropy.coordinates.ITRS` object from this object, use the ``itrs``
    property.
    """

    _location_dtype = np.dtype({"names": ["x", "y", "z"], "formats": [np.float64] * 3})
    _array_dtype = np.dtype((np.float64, (3,)))

    def __new__(cls, *args, **kwargs):

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], GroundLocation):
            # input is a GroundLocation object, copy it
            return args[0].copy()

        if len(args) == 1 and isinstance(args[0], GeodeticLocation):
            # try parsing as geodetic coordinates (lat, lon, alt)
            self = cls.from_geodetic(
                lon=args[0].lon, lat=args[0].lat, height=args[0].height, **kwargs
            )
            return self

        try:
            # try parsing as geocentric cartesian
            self = cls.from_geocentric(*args, **kwargs)
        except (u.UnitsError, TypeError) as exc_geocentric:
            try:
                # try parsing as geodetic coordinates (lat, lon, alt)
                self = cls.from_geodetic(*args, **kwargs)
            except Exception as exc_geodetic:
                raise TypeError(
                    "Coordinates could not be parsed as either "
                    "geocentric or geodetic, with respective "
                    'exceptions "{}" and "{}"'.format(exc_geocentric, exc_geodetic)
                )
        return self

    @classmethod
    def from_geocentric(cls, x, y, z, unit=None, ellipsoid=EARTH_ELLIPSOID_GRS80):
        """
        Location on Central Body (e.g. Earth or Moon), initialized from
        geocentric cartesian coordinates.

        Parameters
        ----------
        x, y, z : `~astropy.units.Quantity` or array_like
            Cartesian coordinates.  If not quantities, ``unit`` should be given.
        unit : `~astropy.units.UnitBase` object or None
            Physical unit of the coordinate values.  If ``x``, ``y``, and/or
            ``z`` are quantities, they will be converted to this unit.
        ellipsoid : CelestialBodyEllipsoid, optional
            Definition of the reference ellipsoid to use
            (default: 'EARTH_ELLIPSOID_GRS80').

        Raises
        ------
        astropy.units.UnitsError
            If the units on ``x``, ``y``, and ``z`` do not match or an invalid
            unit is given.
        ValueError
            If the shapes of ``x``, ``y``, and ``z`` do not match.
        TypeError
            If ``x`` is not a `~astropy.units.Quantity` and no unit is given.
        """
        if unit is None:
            try:
                unit = x.unit
            except AttributeError:
                raise TypeError(
                    "Geocentric coordinates should be Quantities "
                    "unless an explicit unit is given."
                )
        else:
            unit = u.Unit(unit)

        if unit.physical_type != "length":
            raise u.UnitsError("Geocentric coordinates should be in units of length.")

        try:
            x = u.Quantity(x, unit, copy=False)
            y = u.Quantity(y, unit, copy=False)
            z = u.Quantity(z, unit, copy=False)
        except u.UnitsError:
            raise u.UnitsError("Geocentric coordinate units should all be consistent.")

        x, y, z = np.broadcast_arrays(x, y, z)
        struc = np.empty(x.shape, cls._location_dtype)
        struc["x"], struc["y"], struc["z"] = x, y, z
        self = super().__new__(cls, struc, unit, copy=False)
        self._ellipsoid = ellipsoid
        return self

    @classmethod
    def from_geodetic(cls, lon, lat, height=0.0, ellipsoid=EARTH_ELLIPSOID_GRS80):
        """
        Location on Central Body (e.g. Earth or Moon), initialized from
        geodetic coordinates.

        Parameters
        ----------
        lon : `~astropy.coordinates.Longitude` or float
            Earth East longitude.  Can be anything that initialises an
            `~astropy.coordinates.Angle` object (if float, in degrees).
        lat : `~astropy.coordinates.Latitude` or float
            Earth latitude.  Can be anything that initialises an
            `~astropy.coordinates.Latitude` object (if float, in degrees).
        height : `~astropy.units.Quantity` or float, optional
            Height above reference ellipsoid (if float, in meters; default: 0).
        ellipsoid : CelestialBodyEllipsoid, optional
            Definition of the reference ellipsoid to use
            (default: 'EARTH_ELLIPSOID_GRS80').


        Raises
        ------
        astropy.units.UnitsError
            If the units on ``lon`` and ``lat`` are inconsistent with angular
            ones, or that on ``height`` with a length.
        ValueError
            If ``lon``, ``lat``, and ``height`` do not have the same shape, or
            if ``ellipsoid`` is not recognized as among the ones implemented.

        Notes
        -----
        For the conversion to geocentric coordinates, the ERFA routine
        ``gd2gce`` is used.  See https://github.com/liberfa/erfa
        """
        # We use Angle here since there is no need to wrap the longitude -
        # gd2gce will just take cos/sin anyway.  And wrapping might fail
        # on readonly input.
        lon = Angle(lon, u.degree, copy=False)
        lat = Latitude(lat, u.degree, copy=False)
        # don't convert to m by default, so we can use the height unit below.
        if not isinstance(height, u.Quantity):
            height = u.Quantity(height, u.m, copy=False)
        # get geocentric coordinates. Have to give one-dimensional array.
        xyz = erfa.gd2gce(
            ellipsoid.re.to_value(u.m),
            1 / ellipsoid.inv_f.to_value(),
            lon.to_value(u.radian),
            lat.to_value(u.radian),
            height.to_value(u.m),
        )
        self = xyz.ravel().view(cls._location_dtype, cls).reshape(xyz.shape[:-1])
        self._unit = u.m
        self._ellipsoid = ellipsoid
        return self

    @property
    def ellipsoid(self):
        """The default ellipsoid used to convert to geodetic coordinates."""
        return self._ellipsoid

    @property
    def geodetic(self):
        """Convert to geodetic coordinates for the default ellipsoid."""
        return self.to_geodetic()

    def to_geodetic(self, ellipsoid=EARTH_ELLIPSOID_GRS80):
        """Convert to geodetic coordinates.

        Parameters
        ----------
        ellipsoid : CelestialBodyEllipsoid, optional
            Definition of the reference ellipsoid to use
            (default: 'EARTH_ELLIPSOID_GRS80').

        Returns
        -------
        (lon, lat, height) : tuple
            The tuple contains instances of `~astropy.coordinates.Longitude`,
            `~astropy.coordinates.Latitude`, and `~astropy.units.Quantity`


        Notes
        -----
        For the conversion to geodetic coordinates, the ERFA routine
        ``gc2gde`` is used.  See https://github.com/liberfa/erfa
        """

        self_array = self.to(u.meter).view(self._array_dtype, np.ndarray)
        lon, lat, height = erfa.gc2gde(
            ellipsoid.re.to_value(u.m), 1 / ellipsoid.inv_f.to_value(), self_array
        )
        return GeodeticLocation(
            Longitude(
                lon * u.radian, u.degree, wrap_angle=180.0 * u.degree, copy=False
            ),
            Latitude(lat * u.radian, u.degree, copy=False),
            u.Quantity(height * u.meter, self.unit, copy=False),
        )

    @property
    def lon(self):
        """Longitude of the location, for the default ellipsoid."""
        return self.geodetic[0]

    @property
    def lat(self):
        """Longitude of the location, for the default ellipsoid."""
        return self.geodetic[1]

    @property
    def height(self):
        """Height of the location, for the default ellipsoid."""
        return self.geodetic[2]

    # mostly for symmetry with geodetic and to_geodetic.
    @property
    def geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return self.to_geocentric()

    def to_geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return self.x, self.y, self.z

    @property
    def x(self):
        """The X component of the geocentric coordinates."""
        return self["x"]

    @property
    def y(self):
        """The Y component of the geocentric coordinates."""
        return self["y"]

    @property
    def z(self):
        """The Z component of the geocentric coordinates."""
        return self["z"]

    def __getitem__(self, item):
        result = super().__getitem__(item)
        if result.dtype is self.dtype:
            return result.view(self.__class__)
        else:
            return result.view(u.Quantity)

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)
        if hasattr(obj, "_ellipsoid"):
            self._ellipsoid = obj.ellipsoid

    def __len__(self):
        if self.shape == ():
            raise IndexError("0-d GroundLocation arrays cannot be indexed")
        else:
            return super().__len__()

    def _to_value(self, unit, equivalencies=None):
        """Helper method for to and to_value."""
        # Conversion to another unit in both ``to`` and ``to_value`` goes
        # via this routine. To make the regular quantity routines work, we
        # temporarily turn the structured array into a regular one.
        if equivalencies is None:
            equivalencies = []
        array_view = self.view(self._array_dtype, np.ndarray)
        if not equivalencies:
            equivalencies = self._equivalencies
        new_array = self.unit.to(unit, array_view, equivalencies=equivalencies)
        return new_array.view(self.dtype).reshape(self.shape)
