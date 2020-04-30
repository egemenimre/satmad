"""
Terrestrial Intermediate Reference System (TIRS) definition and conversions.

Copyright (C) 2020 Egemen Imre

Licensed under GNU GPL v3.0. See LICENSE.rst for more info.

"""
import astropy.coordinates.representation as r
from astropy import _erfa as erfa
from astropy.coordinates import (BaseCoordinateFrame, TimeAttribute, frame_transform_graph, ITRS,
                                 DynamicMatrixTransform)
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME, get_polar_motion, get_jd12

# __all__ = ['TIRS']


class TIRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the Terrestrial Intermediate Reference System (TIRS).

    References
    ----------
    The definitions and conversions are from IERS Conventions 2010, Chapter 5 [1]_.

    .. [1] IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). (IERS Technical Note ; 36) Frankfurt am Main;
        Verlag des Bundesamts für Kartographie und Geodäsie, 2010. 179 pp., ISBN 3-89888-989-6

    """

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


def _polar_mot_matrix(obstime):
    """
    Form the matrix of polar motion for a given date, IAU 2000.

    The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning that it is the final rotation when computing the
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
    sp = erfa.sp00(*get_jd12(obstime, 'tt'))
    polar_mot_mat = erfa.pom00(xp, yp, sp)

    return polar_mot_mat


@frame_transform_graph.transform(DynamicMatrixTransform, TIRS, ITRS)
def tirs_to_itrs(tirs_coord, itrs_frame):
    tirs_to_itrs_mat = _polar_mot_matrix(tirs_coord.obstime)

    return tirs_to_itrs_mat


@frame_transform_graph.transform(DynamicMatrixTransform, ITRS, TIRS)
def itrs_to_tirs(itrs_coord, tirs_frame):
    itrs_to_tirs_mat = _polar_mot_matrix(itrs_coord.obstime).transpose()

    return itrs_to_tirs_mat
