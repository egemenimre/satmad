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
    ITRS,
    BaseCoordinateFrame,
    CartesianDifferential,
    CartesianRepresentation,
    DynamicMatrixTransform,
    FunctionTransform,
    StaticMatrixTransform,
    TimeAttribute,
    frame_transform_graph,
)
from astropy.coordinates import representation as r
from astropy.coordinates.builtin_frames.utils import (
    DEFAULT_OBSTIME,
    get_jd12,
    get_polar_motion,
)
from astropy.coordinates.matrix_utilities import rotation_matrix

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


# *********** True Equator Mean Equinox Reference System definition and
# conversions. ***********


class TEME(BaseCoordinateFrame):
    """
    A coordinate or frame in the True Equator Mean Equinox Reference System
    (TEME).

    This frame is used as the output of the SGP4 Satellite Propagator.
    This should not be used for any other purpose.
    GCRS should be used wherever possible.

    References
    ----------
    The definitions and conversions are from Fundamentals of Astrodynamics and
    Applications 4th Ed. Section 3.7, pg 231
    [OM1]_.

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


@frame_transform_graph.transform(FunctionTransform, TEME, TIRS)
def teme_to_tirs(teme_coord, tirs_frame):
    # TEME to TIRS basic rotation matrix
    teme_to_pef_mat = rotation_matrix(_gmst82_angle(teme_coord.obstime), axis="z")

    # rotate position vector: TEME to TIRS
    r_tirs = teme_coord.cartesian.transform(teme_to_pef_mat)

    # Check for velocity - skip velocity transform if not present
    if teme_coord.data.differentials:
        # prepare rotation offset: w x r_TIRS
        wxr = CartesianRepresentation(_w).cross(r_tirs)

        # do the velocity rotation and then add rotation offset
        v_tirs = teme_coord.velocity.to_cartesian().transform(teme_to_pef_mat) - wxr
        v_tirs = CartesianDifferential.from_cartesian(v_tirs)

        # Prepare final coord vector with velocity
        tirs_coord = r_tirs.with_differentials(v_tirs)
    else:
        # Prepare final coord vector without velocity
        tirs_coord = r_tirs

    # Add coord data to the existing frame
    return tirs_frame.realize_frame(tirs_coord)


@frame_transform_graph.transform(FunctionTransform, TIRS, TEME)
def tirs_to_teme(tirs_coord, teme_frame):
    # TIRS to TEME basic rotation matrix
    teme_to_pef_mat = rotation_matrix(_gmst82_angle(tirs_coord.obstime), axis="z")
    pef_to_teme_mat = teme_to_pef_mat.transpose()

    # rotate position vector: TIRS to TEME
    r_teme = tirs_coord.cartesian.transform(pef_to_teme_mat)

    # Check for velocity - skip velocity transform if not present
    if tirs_coord.data.differentials:
        # prepare rotation offset: w x r_TIRS
        wxr = CartesianRepresentation(_w).cross(tirs_coord.cartesian)

        # add rotation offset and then do the velocity rotation
        v_teme = (tirs_coord.velocity.to_cartesian() + wxr).transform(pef_to_teme_mat)
        v_teme = CartesianDifferential.from_cartesian(v_teme)

        # Prepare final coord vector with velocity
        teme_coord = r_teme.with_differentials(v_teme)
    else:
        # Prepare final coord vector without velocity
        teme_coord = r_teme

    # Add coord data to the existing frame
    return teme_frame.realize_frame(teme_coord)


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
