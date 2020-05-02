"""
True Equator Mean Equinox Reference System definition and conversions.

Copyright (C) 2020 Egemen Imre

Licensed under GNU GPL v3.0. See LICENSE.rst for more info.

"""
import astropy.coordinates.representation as r
import numpy as np
from astropy import _erfa as erfa
from astropy import units as u
from astropy.coordinates import (BaseCoordinateFrame, TimeAttribute, frame_transform_graph, CartesianRepresentation,
                                 CartesianDifferential, FunctionTransform)
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME, get_jd12
from astropy.coordinates.matrix_utilities import rotation_matrix

from .tirs import TIRS

# __all__ = ['TEME']

# Nominal mean angular velocity of the Earth [rad/s] as per GRS 80. (see IERS TN 36)
_w = np.array([0, 0, 7.292115E-5]) / u.s


class TEME(BaseCoordinateFrame):
    """
    A coordinate or frame in the True Equator Mean Equinox Reference System (TEME).

    This frame is used as the output of the SGP4 Satellite Propagator. This should not be used for any other purpose.
    GCRS should be used wherever possible.

    References
    ----------
    The definitions and conversions are from Fundamentals of Astrodynamics and Applications 4th Ed. Section 3.7, pg 231
    [1]_.

    .. [1] Fundamentals of Astrodynamics and Applications 4th Ed. David A. Vallado. Microcosm Press, 2013.,
        ISBN: 978-11881883203


    """
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential


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
    gmst82 = erfa.gmst82(*get_jd12(obstime, 'ut1')) * u.rad

    return gmst82


@frame_transform_graph.transform(FunctionTransform, TEME, TIRS)
def teme_to_tirs(teme_coord, tirs_frame):
    # TEME to TIRS basic rotation matrix
    teme_to_pef_mat = rotation_matrix(_gmst82_angle(teme_coord.obstime), axis='z')

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
    teme_to_pef_mat = rotation_matrix(_gmst82_angle(tirs_coord.obstime), axis='z')
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
