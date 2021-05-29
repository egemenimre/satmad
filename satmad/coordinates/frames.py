# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Coordinate systems and frames defined by `satmad`.

"""
from abc import ABC

import erfa
import numpy as np
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ICRS,
    ITRS,
    AffineTransform,
    BaseCoordinateFrame,
    CartesianDifferential,
    DynamicMatrixTransform,
    FunctionTransformWithFiniteDifference,
    SkyCoord,
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
from astropy.coordinates.matrix_utilities import matrix_product, rotation_matrix
from astropy.time import Time

from satmad.coordinates.planetary_rot_elems import celestial_body_rot_params

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


class CelestialBodyJ2000Equatorial(BaseCoordinateFrame, ABC):
    """
    A coordinate frame representing the Equatorial System of the
    Celestial Body at J2000 Epoch. This is not defined for the Earth as it has its own
    J2000 Equatorial frame.

    The axis orientations are:

    - x-axis: Along the line formed by the intersection of the body equator and
      the x-y plane of the FK5 system, at the J2000 epoch
    - y-axis: Completes the right-handed set.
    - z-axis: Along the instantaneous body spin axis direction at the J2000 epoch

    This corresponds to "Body Inertial" in GMAT.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    def __new__(cls, *args, **kwargs):
        frame_transform_graph.transform(DynamicMatrixTransform, cls.cb_crs, cls)(
            cb_crs_to_cb_j2000_eq
        )
        frame_transform_graph.transform(DynamicMatrixTransform, cls, cls.cb_crs)(
            cb_j2000_eq_to_cb_crs
        )
        return super().__new__(cls)


def cb_j2000_eq_to_cb_crs(cb_eq_coord, _):
    """Conversion from Local Equatorial at J2000 Epoch Reference System
    of a Central Body to its Celestial Reference System."""

    ra, dec, w = celestial_body_rot_params(None, cb_eq_coord.body_name)

    body_eq_to_body_crs = matrix_product(
        rotation_matrix(90 * u.deg - dec, "x"), rotation_matrix(90 * u.deg + ra, "z")
    ).transpose()

    return body_eq_to_body_crs


def cb_crs_to_cb_j2000_eq(cb_crs_coord, _):
    """Conversion from Celestial Reference System of a Central Body
    to its Local Equatorial at J2000 Epoch Reference System."""

    ra, dec, w = celestial_body_rot_params(None, cb_crs_coord.body_name)

    body_crs_to_body_eq = matrix_product(
        rotation_matrix(90 * u.deg - dec, "x"), rotation_matrix(90 * u.deg + ra, "z")
    )

    return body_crs_to_body_eq


class CelestialBodyTODEquatorial(BaseCoordinateFrame, ABC):
    """
    A coordinate frame representing the True-of-Date Equatorial System of the
    Celestial Body. This is not defined for the Earth as it has its own TOD

    The axis orientations are:

    - x-axis: Along the line formed by the intersection of the body equator
      and the ecliptic planes.
    - y-axis: Completes the right-handed set.
    - z-axis: Normal to the equatorial plane.

    This corresponds to "Equator System" in GMAT.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    def __new__(cls, *args, **kwargs):
        frame_transform_graph.transform(DynamicMatrixTransform, cls.cb_crs, cls)(
            cb_crs_to_cb_tod_eq
        )
        frame_transform_graph.transform(DynamicMatrixTransform, cls, cls.cb_crs)(
            cb_tod_eq_to_cb_crs
        )
        return super().__new__(cls)


def cb_tod_eq_to_cb_crs(cb_eq_coord, _):
    """Conversion from Local True-of-Date Equatorial Reference System of a Central Body
    to its Celestial Reference System."""

    body_eq_to_body_crs = body_crs_to_body_eq_matrix(
        cb_eq_coord.obstime, cb_eq_coord.body_name
    ).transpose()

    return body_eq_to_body_crs


def cb_crs_to_cb_tod_eq(cb_crs_coord, _):
    """Conversion from Celestial Reference System of a Central Body
    to its Local Equatorial True-of-Date Reference System."""

    body_crs_to_body_eq = body_crs_to_body_eq_matrix(
        cb_crs_coord.obstime, cb_crs_coord.body_name
    )

    return body_crs_to_body_eq


def body_crs_to_body_eq_matrix(obstime, body_name):
    """Computes the rotation matrix from Celestial Body CRS to
    Local Equatorial True-of-Date Reference System."""
    ra, dec, w = celestial_body_rot_params(obstime, body_name)

    body_crs_to_body_eq = matrix_product(
        rotation_matrix(90 * u.deg - dec, "x"),
        rotation_matrix(90 * u.deg + ra, "z"),
    )

    return body_crs_to_body_eq


def body_crs_to_body_fixed_matrix(obstime, body_name):
    """Computes the rotation matrix from Celestial Body CRS to
    Body Fixed Reference System."""
    ra, dec, w = celestial_body_rot_params(obstime, body_name)

    body_crs_to_body_fixed = matrix_product(
        rotation_matrix(w, "z"),
        rotation_matrix(90 * u.deg - dec, "x"),
        rotation_matrix(90 * u.deg + ra, "z"),
    )

    return body_crs_to_body_fixed


class CelestialBodyFixed(BaseCoordinateFrame, ABC):
    """
    A coordinate frame representing the Body Fixed System of the
    Celestial Body. This is not defined for the Earth as it has its own ITRS.

    The axis orientations are:

    - x-axis: Time-dependent instantaneous direction
    - y-axis: Completes the right-handed set.
    - z-axis: Normal to the equatorial plane.

    This corresponds to "Celestial Body Fixed System" in GMAT.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    def __new__(cls, *args, **kwargs):
        frame_transform_graph.transform(
            FunctionTransformWithFiniteDifference, cls.cb_crs, cls
        )(cb_crs_to_cb_fixed)
        frame_transform_graph.transform(
            FunctionTransformWithFiniteDifference, cls, cls.cb_crs
        )(cb_fixed_to_cb_crs)
        return super().__new__(cls)


def cb_crs_to_cb_fixed(cb_crs_coord, cb_fixed_coord):
    """Conversion from Local Celestial Body CRS
    to its Body Fixed Reference System."""
    # compute the rotation matrix and rotate the position vector
    cart_repr = cb_crs_coord.cartesian.transform(
        body_crs_to_body_fixed_matrix(cb_crs_coord.obstime, cb_fixed_coord.body_name)
    )
    return cb_fixed_coord.realize_frame(cart_repr)


def cb_fixed_to_cb_crs(cb_fixed_coord, cb_crs_coord):
    """Conversion from Body Fixed Reference System of a Central Body
    to its Local Celestial Body CRS."""
    # compute the rotation matrix and rotate the position vector
    cart_repr = cb_fixed_coord.cartesian.transform(
        body_crs_to_body_fixed_matrix(
            cb_fixed_coord.obstime, cb_fixed_coord.body_name
        ).transpose()
    )

    cb_crs = cb_crs_coord.cb_crs(cart_repr, obstime=cb_fixed_coord.obstime)
    return cb_crs.transform_to(cb_crs_coord)


# @frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, TETE)
# def itrs_to_tete(itrs_coo, tete_frame):
#     # compute the pmatrix, and then multiply by its transpose
#     pmat = tete_to_itrs_mat(itrs_coo.obstime)
#     newrepr = itrs_coo.cartesian.transform(matrix_transpose(pmat))
#     tete = TETE(newrepr, obstime=itrs_coo.obstime)
#
#     # now do any needed offsets (no-op if same obstime)
#     return tete.transform_to(tete_frame)


class CelestialBodyCRS(BaseCoordinateFrame, ABC):
    """
    A coordinate frame in the generic Celestial Reference System (CRS). This CRS is
    derived from ICRS by simply carrying the origin to the origin of the celestial body.

    Uses the `builtin` ephemeris to compute planetary positions in ICRS, unless the
    coordinate definition explicitly specifies an ephemeris type.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    ephemeris_type = "builtin"

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


def cb_crs_to_icrs(cb_crs_coord, _):
    """Conversion from Celestial Reference System of a Central Body to ICRS."""

    if not u.m.is_equivalent(cb_crs_coord.cartesian.x.unit):
        raise u.UnitsError(_NEED_ORIGIN_HINT.format(cb_crs_coord.__class__.__name__))

    if cb_crs_coord.data.differentials:
        # Calculate the barycentric position and velocity (of a solar system body).
        # Uses default ephemeris
        r_icrs, v_icrs = get_body_barycentric_posvel(
            cb_crs_coord.body_name,
            cb_crs_coord.obstime,
            ephemeris=cb_crs_coord.ephemeris_type,
        )

        v_icrs = CartesianDifferential.from_cartesian(v_icrs)

        # Prepare final coord vector with velocity
        icrs_coord = r_icrs.with_differentials(v_icrs)
    else:
        # Calculate the barycentric position ONLY (of a solar system body).
        # Uses default ephemeris. This is faster than the one above with velocities for
        # JPL ephemerides.
        icrs_coord = get_body_barycentric(
            cb_crs_coord.body_name,
            cb_crs_coord.obstime,
            ephemeris=cb_crs_coord.ephemeris_type,
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
            cb_crs_frame.body_name,
            cb_crs_frame.obstime,
            ephemeris=cb_crs_frame.ephemeris_type,
        )

        v_icrs = CartesianDifferential.from_cartesian(v_icrs)

        # Prepare final coord vector with velocity
        cb_crs_coord = (-r_icrs).with_differentials(-v_icrs)
    else:
        # Calculate the barycentric position ONLY (of a solar system body).
        # Uses default ephemeris. This is faster than the one above with velocities for
        # JPL ephemerides.
        cb_crs_coord = -get_body_barycentric(
            cb_crs_frame.body_name,
            cb_crs_frame.obstime,
            ephemeris=cb_crs_frame.ephemeris_type,
        )

    # Return transformation matrix (None) and translation vector (with velocities)
    return None, cb_crs_coord


class MoonCRS(CelestialBodyCRS):
    """Moon Celestial Reference System. This is simply the ICRS shifted to the
    centre of the Moon with the velocity adjusted with respect to the Moon.

    This uses the ephemeris type `jpl` to ensure that Moon velocity is computed
    (`builtin` cannot compute velocity). This is critical for coordinate
    transformations that involve velocity.
    """

    body_name = "Moon"
    ephemeris_type = "jpl"


class MarsCRS(CelestialBodyCRS):
    """Mars Celestial Reference System. This is simply the ICRS shifted to the
    centre of the Mars with the velocity adjusted with respect to the Mars.
    """

    body_name = "Mars"
    ephemeris_type = "jpl"


class MarsTODEquatorial(CelestialBodyTODEquatorial):
    """
    A coordinate frame representing the True-of-Date Equatorial System of Mars.
    """

    body_name = "Mars"
    cb_crs = MarsCRS


class MarsJ2000Equatorial(CelestialBodyJ2000Equatorial):
    """
    A coordinate frame representing the Equatorial System of Mars at J2000 Epoch.
    """

    body_name = "Mars"
    cb_crs = MarsCRS


class MarsBodyFixed(CelestialBodyTODEquatorial):
    """
    A coordinate frame representing the Body Fixed System of Mars.
    """

    body_name = "Mars"
    cb_crs = MarsCRS


class MoonJ2000Equatorial(CelestialBodyJ2000Equatorial):
    """
    A coordinate frame representing the Equatorial System of Moon at J2000 Epoch.

    This corresponds to the Moon Principal Axis (PA) system. This has a fixed and small
    rotational offset with respect to the Mean Earth/rotation (ME) system.
    """

    body_name = "Moon"
    cb_crs = MoonCRS


class MoonTODEquatorial(CelestialBodyTODEquatorial):
    """
    A coordinate frame representing the True-of-Date Equatorial System of Moon.
    """

    body_name = "Moon"
    cb_crs = MoonCRS


class MoonBodyFixed(CelestialBodyTODEquatorial):
    """
    A coordinate frame representing the Body Fixed System of Moon.
    """

    body_name = "Moon"
    cb_crs = MoonCRS


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
    [TCF1] in :doc:`References <../references>`.

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
    [TCF1] in :doc:`References <../references>`.

    """

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


@frame_transform_graph.transform(DynamicMatrixTransform, TIRS, ITRS)
def tirs_to_itrs(tirs_coord, _):
    """Dynamic conversion matrix (Polar Motion) from TIRS to ITRS."""
    tirs_to_itrs_mat = _polar_mot_matrix(tirs_coord.obstime)

    return tirs_to_itrs_mat


@frame_transform_graph.transform(DynamicMatrixTransform, ITRS, TIRS)
def itrs_to_tirs(itrs_coord, _):
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


def init_pvt(frame, time, pos, vel=None, copy=False):
    """
    Convenience method to initialise a `SkyCoord` object using the inputs.

    If velocity input is not given (or is `None`), then velocity input is not set.
    If a velocity value of (0,0,0) is desired, these values should be explicitly
    assigned and the corresponding velocity vector should be given as input.

    Note that, the position vector input may be including the velocity inputs. If
    velocity input is not given, then these velocity values are not overwritten.

    Also note that, `Time`, `CartesianRepresentation` and `CartesianDifferential`
    objects as well as the resulting `SkyCoord` object support multiple value inputs,
    as long as the input sizes match. As an example, it is possible to input
    30 time, position and velocity values, representing discrete points on
    a trajectory.

    Parameters
    ----------
    frame : str or Type[BaseCoordinateFrame]
        Frame name (string) or frame object itself (e.g.GCRS)
    time : Time
        time(s) associated with the coordinates
    pos : CartesianRepresentation
        Position vector(s)
    vel : CartesianDifferential or None
        Velocity vector(s)
    copy : bool, optional
        If `True` (default), a copy of any coordinate data is made.

    Returns
    -------
    SkyCoord
        `SkyCoord` object representing the cartesian input
    """

    if isinstance(frame, str):
        # frame is defined with its name, try to identify it
        rv_frame = frame_transform_graph.lookup_name(frame.lower())
        if not rv_frame:
            raise ValueError(
                f"Frame name {frame} not recognised as a valid frame name."
            )
    else:
        # frame is (hopefully) defined with an object
        rv_frame = frame

    if not isinstance(time, Time):
        raise TypeError(f"Time input ({time}) not recognised as a valid Time object.")

    # All is ready, init the frame object
    try:
        # Certain coordinate frames (notably ICRS) don't like obstime, check for this
        if vel:
            rv_skycoord = rv_frame(
                pos.with_differentials(vel),
                obstime=time,
                representation_type="cartesian",
                differential_type="cartesian",
            )
        else:
            rv_skycoord = rv_frame(
                pos,
                obstime=time,
                representation_type="cartesian",
                differential_type="cartesian",
                copy=copy,
            )

        return SkyCoord(rv_skycoord, copy=copy)
    except TypeError:
        if vel:
            rv_skycoord = rv_frame(
                pos.with_differentials(vel),
                representation_type="cartesian",
                differential_type="cartesian",
            )
        else:
            rv_skycoord = rv_frame(
                pos,
                representation_type="cartesian",
                differential_type="cartesian",
                copy=copy,
            )

        return SkyCoord(
            rv_skycoord,
            obstime=time,
            representation_type="cartesian",
            differential_type="cartesian",
            copy=copy,
        )
