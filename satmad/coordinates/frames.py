# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Coordinate systems and frames defined by `satmad`.

"""
import numpy as np
from astropy import _erfa as erfa
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ICRS,
    ITRS,
    AffineTransform,
    BaseCoordinateFrame,
    CartesianDifferential,
    DynamicMatrixTransform,
    StaticMatrixTransform,
    TimeAttribute,
    frame_transform_graph,
    get_body_barycentric,
    get_body_barycentric_posvel,
)
from astropy.coordinates import representation as r
from astropy.coordinates.builtin_frames.utils import (
    DEFAULT_OBSTIME,
    get_jd12,
    get_polar_motion,
)

from satmad.core.celestial_bodies import MOON

_w = np.array([0, 0, 7.292115e-5]) / u.s
"""Nominal mean angular velocity of the Earth [rad/s] as per GRS 80.
(see IERS TN 36)"""

_FRAME_BIAS_MATRIX = np.array(
    [
        [0.9999999999999942, 7.078279477857338e-8, -8.056217380986972e-8],
        [-7.078279744199198e-8, 0.9999999999999969, -3.306040883980552e-8],
        [8.056217146976134e-8, 3.3060414542221364e-8, 0.9999999999999962],
    ]
)
"""This is the fixed 3x3 frame bias matrix that does the conversion between the
 J2000 and GCRS frames.
"""


class CelestialBodyCRS(BaseCoordinateFrame):
    """
    A coordinate frame in the generic Celestial Reference System (CRS). This CRS is
    derived from ICRS by simply carrying the origin to the origin of the celestial body.
    A specific example is GCRS, where the origin is at the centre of the Earth.

    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    def __new__(cls, *args, **kwargs):
        frame_transform_graph.transform(AffineTransform, ICRS, cls)(icrs_to_cb_crs)
        frame_transform_graph.transform(AffineTransform, cls, ICRS)(cb_crs_to_icrs)
        return super().__new__(cls)


_NEED_ORIGIN_HINT = (
    "The input {0} coordinates do not have length units. This "
    "probably means you created coordinates with lat/lon but "
    "no distance.  Heliocentric<->ICRS transforms cannot "
    "function in this case because there is an origin shift."
)


def cb_crs_to_icrs(cb_crs_coord, icrs_frame):
    """Conversion from Celestial Reference System of a Central Body to ICRS."""

    if not u.m.is_equivalent(cb_crs_coord.cartesian.x.unit):
        raise u.UnitsError(_NEED_ORIGIN_HINT.format(cb_crs_coord.__class__.__name__))

    if cb_crs_coord.data.differentials:
        # Calculate the barycentric position and velocity (of a solar system body).
        # Uses default ephemeris
        r_icrs, v_icrs = get_body_barycentric_posvel(
            cb_crs_coord.body.name, cb_crs_coord.obstime, ephemeris="builtin"
        )

        v_icrs = CartesianDifferential.from_cartesian(v_icrs)

        # Prepare final coord vector with velocity
        icrs_coord = r_icrs.with_differentials(v_icrs)
    else:
        # Calculate the barycentric position ONLY (of a solar system body).
        # Uses default ephemeris. This is faster than the one above with velocities for
        # JPL ephemerides.
        icrs_coord = get_body_barycentric(
            cb_crs_coord.body.name, cb_crs_coord.obstime, ephemeris="builtin"
        )

    # Return transformation matrix (None) and translation vector (with velocities)
    return None, icrs_coord


def icrs_to_cb_crs(icrs_coord, cb_crs_frame):
    """Conversion from ICRS to Celestial Reference System of a Central Body."""

    if not u.m.is_equivalent(icrs_coord.cartesian.x.unit):
        raise u.UnitsError(_NEED_ORIGIN_HINT.format(icrs_coord.__class__.__name__))

    if icrs_coord.data.differentials:
        # Calculate the barycentric position and velocity (of a solar system body).
        # Uses default ephemeris
        r_icrs, v_icrs = get_body_barycentric_posvel(
            cb_crs_frame.body.name, cb_crs_frame.obstime, ephemeris="builtin"
        )

        v_icrs = CartesianDifferential.from_cartesian(v_icrs)

        # Prepare final coord vector with velocity
        cb_crs_coord = (-r_icrs).with_differentials(-v_icrs)
    else:
        # Calculate the barycentric position ONLY (of a solar system body).
        # Uses default ephemeris. This is faster than the one above with velocities for
        # JPL ephemerides.
        cb_crs_coord = -get_body_barycentric(
            cb_crs_frame.body.name, cb_crs_frame.obstime, ephemeris="builtin"
        )

    # Return transformation matrix (None) and translation vector (with velocities)
    return None, cb_crs_coord


class MoonCRS(CelestialBodyCRS):
    """Moon Celestial Reference System. This is simply the ICRS shifted to the
    centre of the Moon."""

    body = MOON


# ******  Mean Pole and Equinox at J2000.0 Reference System (J2000) ******


class J2000(BaseCoordinateFrame):
    """
    A coordinate or frame in the Mean Pole and Equinox at J2000.0 Reference
    System (J2000).

    This coordinate frame is similar to GCRS but rotated by the frame bias.
    This rotation is applicable only to the
    equinox based approach, and is only an approximation. The difference
    betweenGCRS and J2000 is less than 1m for
    the Low Earth Orbit, therefore these two can be used interchangeably with
    a small error.

    This frame is for legacy applications that store data in J2000. GCRS
    should be used wherever possible.

    References
    ----------
    The definitions and conversions are from IERS Conventions 2010, Chapter 5
    [TCF1]_.

    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential


# *********** Terrestrial Intermediate Reference System (TIRS) definition
# and conversions. ***********


class TIRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the Terrestrial Intermediate Reference System
    (TIRS).

    References
    ----------
    The definitions and conversions are from IERS Conventions 2010, Chapter 5
    [TCF1]_.

    """

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


def _gmst82_angle(obstime):
    """
    Universal Time to Greenwich mean sidereal time (IAU 1982 model).

    Parameters
    ----------
    obstime : Time
        time at which the polar motion should be calculated.
    Returns
    -------
    float
        Greenwich mean sidereal time (radians)
    """
    # Get GMST82 angle in rad
    # noinspection PyArgumentList
    gmst82 = erfa.gmst82(*get_jd12(obstime, "ut1")) * u.rad

    return gmst82


def _polar_mot_matrix(obstime):
    """
    Form the matrix of polar motion for a given date, IAU 2000.

    The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning that
    it is the final rotation when computing the
    pointing direction to a celestial source.

    Parameters
    ----------
    obstime : Time
        time at which the polar motion should be calculated.
    Returns
    -------
        3x3 rotation matrix due to polar motion
    """
    # compute the polar motion p-matrix
    xp, yp = get_polar_motion(obstime)
    # noinspection PyArgumentList
    sp = erfa.sp00(*get_jd12(obstime, "tt"))
    polar_mot_mat = erfa.pom00(xp, yp, sp)

    return polar_mot_mat


@frame_transform_graph.transform(DynamicMatrixTransform, TIRS, ITRS)
def tirs_to_itrs(tirs_coord, itrs_frame):
    """Dynamic conversion matrix (Polar Motion) from TIRS to ITRS."""
    tirs_to_itrs_mat = _polar_mot_matrix(tirs_coord.obstime)

    return tirs_to_itrs_mat


@frame_transform_graph.transform(DynamicMatrixTransform, ITRS, TIRS)
def itrs_to_tirs(itrs_coord, tirs_frame):
    """Dynamic conversion matrix (Polar Motion) from ITRS to TIRS."""
    itrs_to_tirs_mat = _polar_mot_matrix(itrs_coord.obstime).transpose()

    return itrs_to_tirs_mat


@frame_transform_graph.transform(StaticMatrixTransform, J2000, GCRS)
def j2000_to_gcrs():
    """Constant conversion matrix (Frame Bias) from J2000 to GCRS."""
    return _FRAME_BIAS_MATRIX


@frame_transform_graph.transform(StaticMatrixTransform, GCRS, J2000)
def gcrs_to_j2000():
    """Constant conversion matrix (Frame Bias) from GCRS to J2000."""
    return _FRAME_BIAS_MATRIX.transpose()
