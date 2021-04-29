# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Methods related to occultations, shadows and basic illumination.

"""
from enum import Enum, auto

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation, SkyCoord, get_body_barycentric

from satmad.core.celestial_bodies import EARTH, SUN


class IlluminationStatus(Enum):
    """Illumination status: UMBRA, PENUMBRA and ILLUMINATED"""

    UMBRA = auto()
    PENUMBRA = auto()
    ILLUMINATED = auto()


def compute_occultation(
    rv_obj, occulting_body=EARTH, illum_body=SUN, ephemeris="builtin"
):
    """
    Computes the instantaneous occultation status and related geometry.

    Based on the paper NASA Technical Paper 3547, "Method for the Calculation of
    Spacecraft Umbra and Penumbra Shadow Terminator Points", Carlos R. Ortiz Longo,
    Steven L. Rickman. Apr 1995.
    <https://ntrs.nasa.gov/api/citations/19950023025/downloads/19950023025.pdf?attachment=true>

    In addition to the paper, this function takes into account the oblateness of the
    occulting body, assuming that the oblateness is along the local Z axis (i.e.,
    equatorial radius is the maximum radius and polar radius is the minimum radius).


    a) Shadow terminators may only be encountered when `proj_distance >= 0`
    b) However,the object will still be illuminated if `penumbra_param > 0`
    c) object in penumbra if  `penumbra_param < 0 < umbra_param`
    d) object in umbra if `0 < umbra_param`
    e) `penumbra_param = 0` is penumbra terminator
    f) `umbra_param = 0` is umbra terminator

    Parameters
    ----------
    rv_obj : SkyCoord
        Position of the object to be tested
    occulting_body : CelestialBody
        Occulting body (defaults to `EARTH`)
    illum_body : CelestialBody
        Illuminating body (defaults to `SUN`)
    ephemeris : str, optional
        Ephemeris to use.  By default uses "builtin".
        Can also be "jpl" (default JPL ephemeris e.g. DE430).


    Returns
    -------

    """

    # TODO little point in calculating where the sun is every 30 secs.
    #      generate a trajectory for the sun and interpolate.
    #      Can even be cached in the Sun object in ICRS

    # TODO check with moon, most probably will have to project to the
    #  inertial frame of the occulting body

    # TODO illumination ratio is probably incorrect

    # compute shadow geometry params
    ksi, kappa, proj_distance, delta = _compute_shadow_geometry(
        rv_obj, occulting_body, illum_body, ephemeris
    )

    # set the penumbra and umbra params
    penumbra_param = delta - kappa
    umbra_param = delta - np.abs(ksi)

    # compute illumination status
    if proj_distance >= 0:
        # Case a: object is on the illuminated side of the occulting body
        illum_ratio = 1.0
        illum_status = IlluminationStatus.ILLUMINATED
        penumbra_param = np.abs(penumbra_param)
        umbra_param = np.abs(umbra_param)

    else:
        # object not on the illuminated side of the occulting body
        if penumbra_param > 0:
            # Case b: object outside penumbra cone
            illum_ratio = 1.0
            illum_status = IlluminationStatus.ILLUMINATED

        else:
            # object is inside penumbra or umbra cones
            if penumbra_param <= 0 < umbra_param:
                # object is inside penumbra
                illum_ratio = (delta - ksi) / (kappa - ksi)
                illum_status = IlluminationStatus.PENUMBRA
            else:
                # object is inside umbra
                illum_ratio = 0.0
                illum_status = IlluminationStatus.UMBRA

    return (
        illum_status,
        illum_ratio,
        penumbra_param.to_value(u.km),
        umbra_param.to_value(u.km),
    )


def _compute_shadow_geometry(
    rv_obj, occulting_body, illum_body=SUN, ephemeris="builtin"
):
    """
    Computes the shadow geometry and related parameters.

    Note the following::

    a) Shadow terminators may only be encountered when `(r_obj . r_illum) < 0`
    b) However,the object will still be illuminated if `delta > kappa`
    c) object in penumbra if  `ksi < delta < kappa`
    d) object in umbra if `delta < ksi`
    e) `delta = kappa` is penumbra terminator
    f) `delta = ksi` is umbra terminator

    Parameters
    ----------
    rv_obj : SkyCoord
        Position of the object to be tested
    occulting_body : CelestialBody
        Occulting body (e.g. Earth)
    illum_body : CelestialBody
        Illuminating body (defaults to `SUN`)
    ephemeris : str, optional
        Ephemeris to use.  By default uses "builtin".
        Can also be "jpl" (default JPL ephemeris e.g. DE430).

    Returns
    -------
    ksi, kappa, proj_distance, delta
        umbra radius, penumbra radius, distance of the object in the shadow axis,
        "flattened" distance (or radius) of the object to the shadow axis in the
        shadow plane

    """

    # define the inertial coord frame of the occulting body,
    # we will execute all computations in this frame
    frame = occulting_body.inert_coord_frame

    # shorthand for object position vector
    r_obj = rv_obj.cartesian.without_differentials()

    # `get_body_barycentric` is in ICRS, convert to the frame needed
    r_occult_body = SkyCoord(
        get_body_barycentric(occulting_body.name, rv_obj.obstime, ephemeris=ephemeris),
        obstime=rv_obj.obstime,
        frame="icrs",
    ).transform_to(frame)
    r_illum_body = SkyCoord(
        get_body_barycentric(illum_body.name, rv_obj.obstime, ephemeris=ephemeris),
        obstime=rv_obj.obstime,
        frame="icrs",
    ).transform_to(frame)
    r_occ_body_to_illum = (
        r_illum_body.cartesian - r_occult_body.cartesian
    )  # s_dot in GMAT

    # unit vector body to illum (e.g. Earth to Sun)
    r_occ_body_to_illum_unit = r_occ_body_to_illum / r_occ_body_to_illum.norm()

    # projection of object position vector on the occ body to illum body
    # (e.g. Earth to Sun).
    # if this is positive, it means that the object is over the illuminated half of
    # the occulting body, therefore is guaranteed to be illuminated
    proj_distance = r_obj.dot(r_occ_body_to_illum_unit)  # -l in GMAT

    # this is practically the normal of the instantaneous "shadow plane"
    rs = r_occ_body_to_illum_unit * proj_distance
    # position vector of the object in the instantaneous "shadow plane"
    r_diff = r_obj - rs

    # compute flattening scale for z axis
    f = 1.0 / occulting_body.ellipsoid.inv_f
    flattening_scale = 1.0 / np.sqrt(1.0 - (2 * f - f * f))

    # scale the z axis in the shadow plane to account for the ellipsoid shape of
    # the shadow - in the local frame z axis is the "squashed" axis
    r_diff = CartesianRepresentation(
        r_diff.x,
        r_diff.y,
        r_diff.z * flattening_scale,
    )

    # distance of position vector on the shadow plane (or radius of the object)
    delta = r_diff.norm()  # d in GMAT (though no flattening there)

    # Illumination depends on the diameters of the Illuminating and Occulting Bodies
    diam_occult = 2 * occulting_body.ellipsoid.re
    diam_illum = 2 * illum_body.ellipsoid.re

    # compute umbra params
    x_u = (diam_occult * r_occ_body_to_illum.norm()) / (diam_illum - diam_occult)
    alpha_u = np.arcsin(diam_occult / (2 * x_u))

    # compute penumbra params
    x_p = (diam_occult * r_occ_body_to_illum.norm()) / (diam_illum + diam_occult)
    alpha_p = np.arcsin(diam_occult / (2 * x_p))

    # compute radii of the shadow terminator points at the shadow plane
    ksi = (x_u - rs.norm()) * np.tan(alpha_u)  # radius of umbra cone - r_U in GMAT
    kappa = (x_p + rs.norm()) * np.tan(alpha_p)  # radius of penumbra cone - r_P in GMAT

    return ksi, kappa, proj_distance, delta
