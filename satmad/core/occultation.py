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

from satmad.core.celestial_bodies import EARTH, SUN
from satmad.utils.discrete_time_events import DiscreteTimeEvents


class IlluminationStatus(Enum):
    """Illumination status: UMBRA, PENUMBRA and ILLUMINATED"""

    UMBRA = auto()
    PENUMBRA = auto()
    ILLUMINATED = auto()


def compute_single_body_occultation_times(
    traj_obj, traj_occult_body, traj_illum_body, occulting_body=EARTH, illum_body=SUN
):

    # TODO there can be multiple occult bodies

    time_list = traj_obj.coord_list.obstime
    illum_pos_list = traj_illum_body(time_list)

    # Do the following for each occulting object

    # init interpolated planet positions
    # this is 10-15% faster than the list comprehension
    occult_pos_list = traj_occult_body(time_list)

    occultation_results = [
        compute_occultation(
            coord,
            occult_pos_list[i],
            illum_pos_list[i],
            occulting_body=occulting_body,
            illum_body=illum_body,
        )
        for i, coord in enumerate(traj_obj.coord_list)
    ]

    umbra_params = np.asarray(
        [result[3].to_value(u.km) for result in occultation_results]
    )
    penumbra_params = np.asarray(
        [result[2].to_value(u.km) for result in occultation_results]
    )

    # ------------------- find umbra times -------------
    umbra_intervals = DiscreteTimeEvents(
        time_list, umbra_params, 0.0, neg_to_pos_is_start=False
    ).start_end_intervals

    # ------------------- find penumbra times -------------
    penumbra_events = DiscreteTimeEvents(
        time_list, penumbra_params, 0.0, neg_to_pos_is_start=False
    )
    # this nominally includes penumbra and umbra times. Subtract the umbra times.
    penumbra_intervals = umbra_intervals.invert().intersect_list(
        penumbra_events.start_end_intervals
    )

    # For each occulting body, return the umbra and penumbra intervals

    return occulting_body.name, umbra_intervals, penumbra_intervals


def compute_occultation(
    rv_obj, r_occult_body, r_illum_body, occulting_body=EARTH, illum_body=SUN
):
    """
    Computes the instantaneous occultation status and related geometry.

    Note that the input position vectors can be in any frame - they will be converted
    to the local inertial frame of the occulting body (e.g. GCRS for the Earth) to
    carry out the calculations. However, it seems to be much faster to convert to the
    correct coordinates outside this method.

    The times of the vectors are not checked to increase speed. It is the user's
    responsibility to make sure that they match.

    Based on the paper NASA Technical Paper 3547, "Method for the Calculation of
    Spacecraft Umbra and Penumbra Shadow Terminator Points", Carlos R. Ortiz Longo,
    Steven L. Rickman. Apr 1995.
    <https://ntrs.nasa.gov/api/citations/19950023025/downloads/19950023025.pdf?attachment=true>

    In addition to the paper, this function takes into account the oblateness of the
    occulting body, assuming that the oblateness is along the local Z axis (i.e.,
    equatorial radius is the maximum radius and polar radius is the minimum radius).

    Also note the following::

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
    r_occult_body : SkyCoord
        Position vector of the Occulting Body
    r_illum_body : SkyCoord
        Position vector of the Illuminating Body
    occulting_body : CelestialBody
        Occulting body (defaults to `EARTH`)
    illum_body : CelestialBody
        Illuminating body (defaults to `SUN`)


    Returns
    -------
    illum_status, illum_ratio, penumbra_param, umbra_param
        Illumination Status, Illumination Ratio (between 0 and 1), penumbra parameter,
        umbra parameter

    Raises
    ------
    ValueError
        Input times of the position vectors do not match
    """

    # TODO check with moon

    # TODO illumination ratio is probably incorrect

    # define the inertial coord frame of the occulting body,
    # we will execute all computations in this frame
    frame = occulting_body.inert_coord_frame

    # convert to occulting body inertial frame if and as needed
    r_obj = __check_frame(rv_obj, frame)
    pos_occult_body = __check_frame(r_occult_body, frame)
    pos_illum_body = __check_frame(r_illum_body, frame)

    # compute pos vector of the illum body
    r_occ_body_to_illum = pos_illum_body - pos_occult_body  # s_dot in GMAT

    # compute shadow geometry params
    ksi, kappa, proj_distance, delta = _compute_shadow_geometry(
        r_obj.xyz.to_value(u.km),
        r_occ_body_to_illum.xyz.to_value(u.km),
        occulting_body.ellipsoid,
        illum_body.ellipsoid,
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
        penumbra_param.to(u.km),
        umbra_param.to(u.km),
    )


def __check_frame(r, tgt_frame):
    """Checks the frame and converts to the target frame if necessary."""
    if r.frame.name != tgt_frame.name:
        pos = r.transform_to(tgt_frame).cartesian.without_differentials()
    else:
        pos = r.cartesian.without_differentials()

    return pos


def _compute_shadow_geometry(
    r_obj, r_occ_body_to_illum, occulting_body_ellipsoid, illum_body_ellipsoid
):
    """
    Computes the shadow geometry and related parameters.

    The method does not use `Quantity` objects to increase processing speed.

    Note the following::

    a) Shadow terminators may only be encountered when `(r_obj . r_illum) < 0`
    b) However,the object will still be illuminated if `delta > kappa`
    c) object in penumbra if  `ksi < delta < kappa`
    d) object in umbra if `delta < ksi`
    e) `delta = kappa` is penumbra terminator
    f) `delta = ksi` is umbra terminator

    Parameters
    ----------
    r_obj : ndarray
        Position vector of the object to be tested in the local inertial frame of the
        Occulting Body [km]
    r_occ_body_to_illum : ndarray
        Position vector of the Illuminating Body in the local inertial frame of the
        Occulting Body
    occulting_body_ellipsoid : CelestialBodyEllipsoid
        Ellipsoid parameters of the Occulting Body
    illum_body_ellipsoid : CelestialBodyEllipsoid
        Ellipsoid parameters of the Illuminating Body

    Returns
    -------
    ksi, kappa, proj_distance, delta
        umbra radius, penumbra radius, distance of the object in the shadow axis,
        "flattened" distance (or radius) of the object to the shadow axis in the
        shadow plane

    """

    # unit vector body to illum (e.g. Earth to Sun)
    r_occ_body_to_illum_norm = np.linalg.norm(r_occ_body_to_illum)
    r_occ_body_to_illum_unit = r_occ_body_to_illum / r_occ_body_to_illum_norm

    # projection of object position vector on the occ body to illum body
    # (e.g. Earth to Sun).
    # if this is positive, it means that the object is over the illuminated half of
    # the occulting body, therefore is guaranteed to be illuminated
    proj_distance = r_obj.dot(r_occ_body_to_illum_unit)  # -l in GMAT

    # this is practically the normal of the instantaneous "shadow plane"
    rs = r_occ_body_to_illum_unit * proj_distance
    rs_norm = np.linalg.norm(rs)
    # position vector of the object in the instantaneous "shadow plane"
    r_diff = r_obj - rs

    # compute flattening scale for z axis
    f = 1.0 / occulting_body_ellipsoid.inv_f.value
    flattening_scale = 1.0 / np.sqrt(1.0 - (2 * f - f * f))

    # scale the z axis in the shadow plane to account for the ellipsoid shape of
    # the shadow - in the local frame z axis is the "squashed" axis
    r_diff[2] = r_diff[2] * flattening_scale

    # distance of position vector on the shadow plane (or radius of the object)
    delta = np.linalg.norm(r_diff)  # d in GMAT (though no flattening there)

    # Illumination depends on the diameters of the Illuminating and Occulting Bodies
    diam_occult = 2 * occulting_body_ellipsoid.re.to_value(u.km)
    diam_illum = 2 * illum_body_ellipsoid.re.to_value(u.km)

    # compute umbra params
    x_u = (diam_occult * r_occ_body_to_illum_norm) / (diam_illum - diam_occult)
    alpha_u = np.arcsin(diam_occult / (2 * x_u))

    # compute penumbra params
    x_p = (diam_occult * r_occ_body_to_illum_norm) / (diam_illum + diam_occult)
    alpha_p = np.arcsin(diam_occult / (2 * x_p))

    # compute radii of the shadow terminator points at the shadow plane
    ksi = (x_u - rs_norm) * np.tan(alpha_u)  # radius of umbra cone - r_U in GMAT
    kappa = (x_p + rs_norm) * np.tan(alpha_p)  # radius of penumbra cone - r_P in GMAT

    return ksi * u.km, kappa * u.km, proj_distance * u.km, delta * u.km
