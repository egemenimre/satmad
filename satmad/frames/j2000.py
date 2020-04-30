"""
Mean Pole and Equinox at J2000.0 Reference System definition and conversions.

Copyright (C) 2020 Egemen Imre

Licensed under GNU GPL v3.0. See LICENSE.rst for more info.

"""
import astropy.coordinates.representation as r
import numpy as np
from astropy.coordinates import (BaseCoordinateFrame, TimeAttribute, frame_transform_graph, StaticMatrixTransform, GCRS)

# __all__ = ['J2000']

from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME

_FRAME_BIAS_MATRIX = np.array(
    [[0.9999999999999942, 7.078279477857338E-8, -8.056217380986972E-8],
     [-7.078279744199198E-8, 0.9999999999999969, -3.306040883980552E-8],
     [8.056217146976134E-8, 3.3060414542221364E-8, 0.9999999999999962]])

"""This is the fixed 3x3 frame bias matrix that does the conversion between the J2000 and GCRS frames. 
"""


class J2000(BaseCoordinateFrame):
    """
    A coordinate or frame in the Mean Pole and Equinox at J2000.0 Reference System (J2000).

    This coordinate frame is similar to GCRS but rotated by the frame bias. This rotation is applicable only to the
    equinox based approach, and is only an approximation. The difference between GCRS and J2000 is less than 1m for
    the Low Earth Orbit, therefore these two can be used interchangeably with a small error.

    This frame is for legacy applications that store data in J2000. GCRS should be used wherever possible.

    References
    ----------
    The definitions and conversions are from IERS Conventions 2010, Chapter 5 [1]_.

    .. [1] IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). (IERS Technical Note ; 36)
        Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 2010. 179 pp., ISBN 3-89888-989-6


    """
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential


@frame_transform_graph.transform(StaticMatrixTransform, J2000, GCRS)
def j2000_to_gcrs():
    return _FRAME_BIAS_MATRIX


@frame_transform_graph.transform(StaticMatrixTransform, GCRS, J2000)
def gcrs_to_j2000():
    return _FRAME_BIAS_MATRIX.transpose()
