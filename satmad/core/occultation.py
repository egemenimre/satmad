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

from satmad.coordinates.trajectory import Trajectory
from satmad.core.celestial_bodies import EARTH, SUN
from satmad.utils.discrete_time_events import DiscreteTimeEvents
from satmad.utils.timeinterval import TimeIntervalList


class IlluminationStatus(Enum):
    """Illumination status: UMBRA, PENUMBRA and ILLUMINATED"""

    UMBRA = auto()
    PENUMBRA = auto()
    ILLUMINATED = auto()


def multi_body_occultation_intervals(
    traj_obj, occult_bodies, traj_illum_body, illum_body=SUN
):
    """
    Computes the occultation intervals of an object for multiple occulting bodies
    (e.g. Earth and Moon).

    The method is initialised with a dictionary containing occulting bodies and their
    trajectories::

        occult_traj_Earth = Trajectory(
            SkyCoord(
                get_body_barycentric(EARTH.name, sparse_time_list, ephemeris="jpl"),
                obstime=sparse_time_list,
                frame="icrs",
                representation_type="cartesian",
                differential_type="cartesian",
            ).transform_to("gcrs")
        )

        occult_traj_Moon = Trajectory(
            SkyCoord(
                get_body_barycentric(MOON.name, sparse_time_list, ephemeris="jpl"),
                obstime=sparse_time_list,
                frame="icrs",
                representation_type="cartesian",
                differential_type="cartesian",
            ).transform_to("MoonCRS")
        )

        occult_bodies = {EARTH: occult_traj_Earth, MOON: occult_traj_Moon}

    Note that the input trajectories can be in any frame - they will be converted
    to the local inertial frame of the occulting body (e.g. GCRS for the Earth) to
    carry out the calculations. However, it seems to be much faster to convert to the
    correct coordinates outside this method.

    Parameters
    ----------
    traj_obj : Trajectory
        Trajectory of the object to be checked (e.g. satellite)
    occult_bodies : dict
        Dictionary containing occulting bodies and trajectories
    traj_illum_body : Trajectory
        Trajectory of the illuminating body (e.g. Sun)
    illum_body : CelestialBody
        Properties of the illuminating body (its name and ellipsoid properties are used)

    Returns
    -------

    """
    # This is the complete set of points where the intervals will be searched
    time_list = traj_obj.coord_list.obstime

    # TODO check intervals match? Occult and Illum Traj must be equal or larger than
    #  object trajectory

    output_dict = {}

    # Compute intervals for each occulting object
    for occulting_body, traj_occult_body in occult_bodies.items():

        # regenerate interpolated planet positions in the coordinate frame necessary
        # make sure no unnecessary coord transformations are taking place
        occult_pos_list = __generate_pos_list(
            traj_occult_body, time_list, occulting_body
        )
        illum_pos_list = __generate_pos_list(traj_illum_body, time_list, occulting_body)
        if traj_obj.coord_list.frame.name != occulting_body.inert_coord_frame.name:
            traj_pos_list = traj_obj.coord_list.transform_to(
                occulting_body.inert_coord_frame
            )
        else:
            traj_pos_list = traj_obj.coord_list

        (
            occulting_body_name,
            umbra_intervals,
            penumbra_intervals,
        ) = _single_body_occultation_intervals(
            traj_pos_list,
            occult_pos_list,
            illum_pos_list,
            occulting_body=occulting_body,
            illum_body=illum_body,
        )

        # check the intervals for the spurious "positive projected distance" cases
        # check the umbra intervals
        umbra_intervals = _generate_corrected_intervals(
            umbra_intervals,
            traj_obj,
            traj_occult_body,
            occulting_body,
            traj_illum_body,
            illum_body,
        )
        # check the penumbra intervals
        penumbra_intervals = _generate_corrected_intervals(
            penumbra_intervals,
            traj_obj,
            traj_occult_body,
            occulting_body,
            traj_illum_body,
            illum_body,
        )

        output_dict[occulting_body_name] = (umbra_intervals, penumbra_intervals)

    return output_dict


def __generate_pos_list(traj_body, time_list, occulting_body):
    """Checks whether the coordinates are in the inertial frame of the occulting body,
    transform them if necessary."""
    if traj_body.coord_list.frame.name != occulting_body.inert_coord_frame.name:
        pos_list = Trajectory(
            traj_body.coord_list.transform_to(occulting_body.inert_coord_frame)
        )(time_list)
        return pos_list
    else:
        return traj_body(time_list)


def _generate_corrected_intervals(
    event_intervals,
    traj_obj,
    traj_occult_body,
    occulting_body,
    traj_illum_body,
    illum_body,
):
    """Checks the intervals to make sure they don't contain false events with zero
    umbra/penumbra param but on the illuminated side of the occulting body."""
    checked_intervals = [
        interval
        for interval in event_intervals.intervals
        if not _is_obj_on_illum_side(
            interval,
            traj_obj,
            traj_occult_body,
            occulting_body,
            traj_illum_body,
            illum_body,
        )
    ]

    return TimeIntervalList(
        checked_intervals, start_valid=event_intervals.valid_interval, replicate=False
    )


def _is_obj_on_illum_side(
    interval, traj_obj, traj_occult_body, occulting_body, traj_illum_body, illum_body
):
    """Checks the interval to make sure they it is on the illuminated side of
    the occulting body."""
    obstime = interval.start + interval.duration * 0.5
    (
        illum_status,
        illum_ratio,
        penumbra_param,
        umbra_param,
        proj_distance,
    ) = compute_occultation(
        traj_obj(obstime),
        traj_occult_body(obstime),
        traj_illum_body(obstime),
        occulting_body=occulting_body,
        illum_body=illum_body,
    )

    if proj_distance >= 0:
        return True
    else:
        return False


def occultation_intervals(
    traj_obj, traj_occult_body, traj_illum_body, occulting_body=EARTH, illum_body=SUN
):
    """
    Computes the occultation intervals of an object for a single occulting body.

    Note that the input trajectories can be in any frame - they will be converted
    to the local inertial frame of the occulting body (e.g. GCRS for the Earth) to
    carry out the calculations. However, it seems to be much faster to convert to the
    correct coordinates outside this method.

    This is a thin wrapper around the `multi_body_occultation_intervals()` method for
    convenience.

    Parameters
    ----------
    traj_obj : Trajectory
        Trajectory of the object to be checked (e.g. satellite)
    traj_occult_body : Trajectory
        Trajectory of the occulting body (e.g. Earth)
    traj_illum_body : Trajectory
        Trajectory of the illuminating body (e.g. Sun)
    occulting_body : CelestialBody
        Properties of the occulting body (its name and ellipsoid properties are used)
    illum_body : CelestialBody
        Properties of the illuminating body (its name and ellipsoid properties are used)

    Returns
    -------
    TimeIntervalList, TimeIntervalList
        Umbra intervals, penumbra intervals
    """

    occult_bodies = {occulting_body: traj_occult_body}

    output_dict = multi_body_occultation_intervals(
        traj_obj, occult_bodies, traj_illum_body, illum_body=SUN
    )

    return output_dict[occulting_body.name]


def _single_body_occultation_intervals(
    obj_pos_list,
    occult_pos_list,
    illum_pos_list,
    occulting_body=EARTH,
    illum_body=SUN,
):
    """
    Computes the occultation intervals for a single body. The times of the coordinate
    lists of the object, occulting body and illuminating body should match - this is
    not checked within the method.

    Parameters
    ----------
    obj_pos_list : SkyCoord
        List of coordinates of the object to be checked (e.g. satellite)
    occult_pos_list: SkyCoord
        List of coordinates of the occulting body (e.g. Earth)
    illum_pos_list: SkyCoord
        List of coordinates of the illuminating body (e.g. Sun)
    occulting_body : CelestialBody
        Properties of the occulting body (its name and ellipsoid properties are used)
    illum_body : CelestialBody
        Properties of the illuminating body (its name and ellipsoid properties are used)

    Returns
    -------
    str, TimeIntervalList, TimeIntervalList
        Name of the occulting body, umbra intervals, penumbra intervals
    """

    time_list = obj_pos_list.obstime

    # compute occultation results for each point in time
    occultation_results = [
        compute_occultation(
            coord,
            occult_pos_list[i],
            illum_pos_list[i],
            occulting_body=occulting_body,
            illum_body=illum_body,
        )
        for i, coord in enumerate(obj_pos_list)
    ]

    # ------------------- find umbra times -------------
    umbra_params = np.asarray(
        [result[3].to_value(u.km) for result in occultation_results]
    )
    umbra_intervals = DiscreteTimeEvents(
        time_list, umbra_params, 0.0, neg_to_pos_is_start=False
    ).start_end_intervals

    # plot umbra params if required
    # from satmad.plots.basic_plots import plot_time_param
    #
    # plot_time_param(time_list, umbra_params, x_rotation=True)
    # ------------------- find penumbra times -------------
    penumbra_params = np.asarray(
        [result[2].to_value(u.km) for result in occultation_results]
    )
    penumbra_events = DiscreteTimeEvents(
        time_list, penumbra_params, 0.0, neg_to_pos_is_start=False
    )

    # this nominally includes penumbra and umbra times. Subtract the umbra times.
    penumbra_intervals = umbra_intervals.invert().intersect_list(
        penumbra_events.start_end_intervals
    )

    # Return the umbra and penumbra intervals
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
    illum_status, illum_ratio, penumbra_param, umbra_param, proj_distance
        Illumination Status, Illumination Ratio (between 0 and 1), penumbra parameter,
        umbra parameter, projected distance

    Raises
    ------
    ValueError
        Input times of the position vectors do not match
    """

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
        proj_distance.to(u.km),
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
