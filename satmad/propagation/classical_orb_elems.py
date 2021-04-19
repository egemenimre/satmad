# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Classical Orbital Elements definitions.

"""
from abc import ABC, abstractmethod

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from astropy.time import Time
from astropy.units import Quantity

from satmad.coordinates.frames import init_pvt
from satmad.core.celestial_bodies import EARTH


class AbstractKeplerianOrbitElements(ABC):
    """
    Classical (Keplerian) Orbital Elements that should be overridden with the specific
    type (Mean, Osculating etc.).

    Parameters
    ----------
    central_body : CelestialBody
        celestial body (e.g. Earth) around which the orbit is defined
    epoch : Time
        Epoch Time corresponding to the orbital elements
    sm_axis: Quantity
        semimajor axis of the orbit [km]
    inclination : Quantity
        inclination of the orbit [rad]
    raan : Quantity
        right ascension of ascending node (RAAN) of the orbit [rad]
    eccentricity : Quantity or float
        eccentricity of the orbit [dimensionless]
    arg_periapsis : Quantity
        argument of periapsis [rad]
    true_anomaly : Quantity
        true anomaly of the orbit [rad]

    Raises
    ------
    ValueError
        Parabolic orbits or singularity
    """

    _arg_periapsis: u.rad
    _raan: u.rad
    _true_anomaly: u.rad

    def __init__(
        self,
        epoch,
        sm_axis,
        eccentricity,
        inclination,
        raan,
        arg_periapsis,
        true_anomaly,
        central_body=EARTH,
    ):

        # Fill init values with setters
        # Setters carry out the type and value checks
        # Ignore codestyle checks as mypy thinks these are read only
        self.central_body = central_body
        self.epoch = epoch
        self.sm_axis = sm_axis  # type: ignore
        self.eccentricity = eccentricity
        self.inclination = inclination  # type: ignore
        self.raan = raan  # type: ignore
        self.arg_periapsis = arg_periapsis  # type: ignore
        self.true_anomaly = true_anomaly  # type: ignore

        # check initialisation
        self._self_check()

    def _self_check(self):
        """Check the elements for errors. """
        # If conic section is too close to singular, raise error.
        if np.abs(self.sm_axis * (1 - self.eccentricity)) < 1e-3 * u.km:
            raise ValueError(
                "Orbit periapsis too close to zero. Singular conic sections are not allowed."
            )

        # If magnitude of the position vector is infinite, raise error.
        if 1 + self.eccentricity * np.cos(self.true_anomaly) < 1e-30:
            raise ValueError("Magnitude of the position vector is near infinite.")

        # If ecc_mag too close to 1 then the orbit is parabolic.
        if np.abs(1 - self.eccentricity) < 1e-7:
            raise ValueError(
                "Orbit eccentricity is close to 1. Parabolic orbits are not allowed."
            )

        # Check true anomaly in the hyperbolic case
        nu_rad = self.true_anomaly.to_value(u.rad)
        nu_rad = nu_rad if nu_rad < np.pi else nu_rad - 2 * np.pi
        if self.eccentricity > 1 and np.abs(nu_rad) >= np.pi - np.arccos(
            1 / self.eccentricity
        ):
            raise ValueError(
                "True anomaly value physically impossible for the sm_axis and ecc values."
            )

    @classmethod
    @abstractmethod
    def from_cartesian(cls, init_coords, central_body=EARTH):
        """
        Generates Keplerian Elements from cartesian coordinates in the
        inertial frame of the celestial body.

        If the initial coordinate is in a different frame than the
        inertial frame of the celestial body, it is automatically converted to
        the proper frame for conversion.

        Parameters
        ----------
        init_coords : SkyCoord
            Initial coordinates (the first value is used)
        central_body : CelestialBody
            celestial body (e.g. Earth) around which the orbit is defined

        Returns
        -------
        OsculatingKeplerianOrbElems
            Osculating classical (or Keplerian) orbital elements

        Raises
        ------
        ValueError
            Parabolic orbits or singularity
        """

        pass

    @abstractmethod
    def to_cartesian(self):
        """
        Converts the orbital elements to the cartesian coordinates in the
        local inertial frame of the central body.

        Returns
        -------
        SkyCoord
            cartesian coordinates in the local inertial frame of the central body

        Raises
        ------
        ValueError
            Parabolic orbits or singularity
        """
        pass

    @property
    def epoch(self) -> Time:
        """Returns the epoch time associated with the orbital parameters.
        Setter replicates the time with `copy=False`."""
        return self._epoch

    @epoch.setter
    def epoch(self, epoch: Time):
        self._epoch = epoch.replicate()

    @property  # type: ignore
    @u.quantity_input
    def sm_axis(self) -> u.km:
        """Gets the semimajor axis [km]."""
        return self._sm_axis

    @sm_axis.setter  # type: ignore
    @u.quantity_input(sm_axis="length")
    def sm_axis(self, sm_axis):
        self._sm_axis = sm_axis

    @property
    def eccentricity(self) -> float:
        """Gets and sets the eccentricity of the orbit.

        Input can be a Quantity (dimensionless) or float. Getter value is a float.

        Eccentricity should be in range 0 <= e < 1.0.
        Raises a `ValueError` otherwise."""
        return self._eccentricity

    @eccentricity.setter
    def eccentricity(self, e):

        # check eccentricity as Quantity and reduce to value
        if isinstance(e, Quantity):
            e = e.value

        # If ecc_mag too close to 1 then the orbit is parabolic.
        if np.abs(1 - e) < 1e-7:
            raise ValueError(
                "Orbit eccentricity is close to 1. Parabolic orbits are not allowed."
            )

        if 0 <= e:
            self._eccentricity = e
        else:
            raise ValueError(
                f"Given argument ({e}) is an invalid "
                f"Eccentricity value, "
                f"only values in range 0 <= e < 1.0 are allowed."
            )

    @property  # type: ignore
    @u.quantity_input
    def inclination(self) -> u.rad:
        """Gets and sets the inclination of the orbit [rad].

        Inclination should be in range 0 <= om < PI.
        Raises a `ValueError` otherwise.
        """
        return self._incl

    @inclination.setter  # type: ignore
    @u.quantity_input(incl="angle")
    def inclination(self, incl):

        if 0 <= incl.to_value(u.rad) < np.pi:
            self._incl = incl
        else:
            raise ValueError(
                f"Given argument ({incl}) is an invalid "
                f"Inclination value, "
                f"only values in range 0 <= i < PI are allowed."
            )

    @property  # type: ignore
    @u.quantity_input
    def raan(self) -> u.rad:
        """Gets and sets the right ascension of ascending node (RAAN)
        of the orbit [rad].

        RAAN should be in range 0 <= om < 2*PI. Any input outside this range will
        be forced into this range."""
        return self._raan

    @raan.setter  # type: ignore
    @u.quantity_input(om="angle")
    def raan(self, om):
        self._raan = _force_angle_range(om)

    @property  # type: ignore
    @u.quantity_input
    def arg_periapsis(self) -> u.rad:
        """Gets and sets the argument of periapsis [rad].

        Argument of periapsis should be in range 0 <= argp < 2*PI.
        Any input outside this range will
        be forced into this range."""
        return self._arg_periapsis

    @arg_periapsis.setter  # type: ignore
    @u.quantity_input(argp="angle")
    def arg_periapsis(self, argp):
        self._arg_periapsis = _force_angle_range(argp)

    @property  # type: ignore
    @u.quantity_input
    def true_anomaly(self) -> u.rad:
        """Gets and sets the true anomaly of the orbit [rad].

        True Anomaly should be in range 0 <= true_anomaly < 2*PI.
        Any input outside this range will
        be forced into this range."""
        return self._true_anomaly

    @true_anomaly.setter  # type: ignore
    @u.quantity_input(true_anomaly="angle")
    def true_anomaly(self, true_anomaly):
        self._true_anomaly = _force_angle_range(true_anomaly)

    @property
    def mean_motion(self):
        """Computes the mean motion of the orbit [rad/sec].

        Returns
        -------
        Quantity
            mean motion in u.rad / u.s
        """
        return np.sqrt(self.central_body.mu / np.abs(self.sm_axis ** 3)).to(
            u.rad / u.s, equivalencies=u.dimensionless_angles()
        )

    @property  # type: ignore
    @u.quantity_input
    def period(self) -> u.s:
        """Computes the satellite period [sec]."""
        return 2 * np.pi * u.rad / self.mean_motion

    @property  # type: ignore
    @u.quantity_input
    def periapsis(self) -> u.km:
        """Computes the periapsis distance [km]."""
        return self.sm_axis * (1 - self.eccentricity)

    @property  # type: ignore
    @u.quantity_input
    def apoapsis(self) -> u.km:
        """Computes the apoapsis distance [km]."""
        return self.sm_axis * (1 + self.eccentricity)


class OsculatingKeplerianOrbElems(AbstractKeplerianOrbitElements):
    """
    Osculating Classical (Keplerian) Orbital Elements in the local inertial frame.

    By definition, this uses a two-body potential for computations. Therefore, this is
    not a "mean elements" model such as a TLE. Over a trajectory generated with a
    two-body force model, the orbital elements (apart from true anomaly) should stay
    constant, limited by the accuracy of the trajectory generation algorithm.

    Osculating orbital elements should be used with care,
    particularly in contexts where instantaneous orbital elements are computed on a
    trajectory generated by a non-two-body model (e.g. including geopotentials). The
    orbital elements will *not* stay constant along the trajectory, simply because the
    force model over successive points are not strictly two-body.

    Parameters
    ----------
    central_body : CelestialBody
        celestial body (e.g. Earth) around which the orbit is defined
    epoch : Time
        Epoch Time corresponding to the orbital elements
    sm_axis: Quantity
        semimajor axis of the orbit [km]
    inclination : Quantity
        inclination of the orbit [rad]
    raan : Quantity
        right ascension of ascending node (RAAN) of the orbit [rad]
    eccentricity : Quantity or float
        eccentricity of the orbit [dimensionless]
    arg_periapsis : Quantity
        argument of periapsis [rad]
    true_anomaly : Quantity
        true anomaly of the orbit [rad]

    Raises
    ------
    ValueError
        Parabolic orbits or singularity
    """

    @u.quantity_input(
        sm_axis="length",
        inclination="angle",
        raan="angle",
        arg_periapsis="angle",
        true_anomaly="angle",
    )
    def __init__(
        self,
        epoch,
        sm_axis,
        eccentricity,
        inclination,
        raan,
        arg_periapsis,
        true_anomaly,
        central_body=EARTH,
    ):
        super().__init__(
            epoch,
            sm_axis,
            eccentricity,
            inclination,
            raan,
            arg_periapsis,
            true_anomaly,
            central_body=central_body,
        )

    @classmethod
    def from_cartesian(cls, init_coords, central_body=EARTH):
        """
        Generates Osculating Keplerian Elements from cartesian coordinates in the
        inertial frame of the celestial body.

        If the initial coordinate is in a different frame than the
        inertial frame of the celestial body, it is automatically converted to
        the proper frame for conversion.

        Parameters
        ----------
        init_coords : SkyCoord
            Initial coordinates (the first value is used)
        central_body : CelestialBody
            celestial body (e.g. Earth) around which the orbit is defined

        Returns
        -------
        OsculatingKeplerianOrbElems
            Osculating classical (or Keplerian) orbital elements

        Raises
        ------
        ValueError
            Parabolic orbits or singularity
        """

        # Check initial conditions
        if init_coords.isscalar:
            # Single SkyCoord element available
            init_coords_converted = init_coords
        else:
            # Multiple SkyCoord elements available, use first
            init_coords_converted = init_coords[0]

        # Check for the initial coordinate system, convert if necessary
        if central_body.inert_coord != init_coords_converted.frame.name:
            init_coords_converted = init_coords_converted.transform_to(
                central_body.inert_coord
            )

        # Compute orbital elements
        orb_elems = _rv_to_keplerian_elems(
            init_coords_converted.obstime,
            init_coords_converted.cartesian.without_differentials(),
            init_coords_converted.velocity,
            central_body,
        )

        return orb_elems

    def to_cartesian(self):
        """
        Converts the orbital elements to the cartesian coordinates in the
        local inertial frame of the central body.

        Returns
        -------
        SkyCoord
            cartesian coordinates in the local inertial frame of the central body

        Raises
        ------
        ValueError
            Parabolic orbits or singularity
        """

        # shorthands
        a = self.sm_axis
        i = self.inclination
        e = self.eccentricity
        argp = self.arg_periapsis
        raan = self.raan
        nu = self.true_anomaly
        node = argp + nu
        mu = self.central_body.mu

        p = a * (1 - e ** 2)

        r = p / (1 + e * np.cos(nu))

        # compute position components
        x = r * (np.cos(node) * np.cos(raan) - np.cos(i) * np.sin(node) * np.sin(raan))
        y = r * (np.cos(node) * np.sin(raan) + np.cos(i) * np.sin(node) * np.cos(raan))
        z = r * (np.sin(node) * np.sin(i))

        # If orbit is too close to parabolic, raise error.
        if np.abs(p) < 1e-30 * u.km:
            raise ValueError(
                "Orbit eccentricity is close to 1. Parabolic orbits are not allowed."
            )

        sqrt_mup = np.sqrt(mu / p)

        # compute velocity components
        v_x = sqrt_mup * (np.cos(nu) + e) * (
            -np.sin(argp) * np.cos(raan) - np.cos(i) * np.sin(raan) * np.cos(argp)
        ) - sqrt_mup * np.sin(nu) * (
            np.cos(argp) * np.cos(raan) - np.cos(i) * np.sin(raan) * np.sin(argp)
        )

        v_y = sqrt_mup * (np.cos(nu) + e) * (
            -np.sin(argp) * np.sin(raan) + np.cos(i) * np.cos(raan) * np.cos(argp)
        ) - sqrt_mup * np.sin(nu) * (
            np.cos(argp) * np.sin(raan) + np.cos(i) * np.cos(raan) * np.sin(argp)
        )
        v_z = sqrt_mup * (
            (np.cos(nu) + e) * np.sin(i) * np.cos(argp)
            - np.sin(i) * np.sin(nu) * np.sin(argp)
        )

        # collate the SkyCoord object
        time = self.epoch.replicate()
        r = CartesianRepresentation([x, y, z], unit=u.km)
        v = CartesianDifferential(
            [v_x, v_y, v_z],
            unit=u.km / u.s,
        )

        rv_init = init_pvt(self.central_body.inert_coord, time, r.with_differentials(v))

        return rv_init


def _rv_to_keplerian_elems(epoch, r, v, central_body):
    """
    Converts position and velocity vectors to an Osculating Keplerian Orbital Elems
    object.

    This section is adapted from GMAT Math Specifications Draft for R2020a Section 3.1.2
    (6 Mar 2020) http://gmat.sourceforge.net/doc/R2020a/GMATMathSpec.pdf.

    Parameters
    ----------
    epoch : Time
        Epoch Time corresponding to the orbital elements
    r: CartesianRepresentation
        Position vector
    v: CartesianDifferential
        Velocity vector
    central_body : CelestialBody
        celestial body (e.g. Earth) around which the orbit is defined

    Returns
    -------
    OsculatingKeplerianOrbElems
        Osculating classical (or Keplerian) orbital elements

    Raises
    ------
    ValueError
        Parabolic orbits or singularity
    """

    # shorthands
    mu = central_body.mu
    r_mag = r.norm()
    v_mag = v.norm()

    # ang mom
    h = r.cross(v)
    h_mag = h.norm()

    # k-axis unit vector (Z-axis of the coord frame)
    k = CartesianRepresentation([0, 0, 1])

    # orbit normal
    n = k.cross(h)
    n_mag = n.norm()

    # find eccentricity vector and magnitude
    ecc = ((v_mag ** 2 - mu / r_mag) * r - r.dot(v) * v) / mu
    ecc_mag = ecc.norm().decompose()

    # Find two-body energy
    energy = 0.5 * v_mag ** 2 - mu / r_mag

    # If ecc_mag too close to 1  then the orbit is parabolic. Abort the conversion.
    if np.abs(1 - ecc_mag) < 1e-7:
        raise ValueError(
            "Orbit eccentricity is equal to 1. Parabolic orbits are not allowed."
        )

    # compute sm_axis
    sm_axis = -mu / (2 * energy)

    # If conic section is too close to singular, then conversion is to be aborted.
    if np.abs(sm_axis * (1 - ecc_mag)) < 1e-3 * u.km:
        raise ValueError(
            "Orbit periapsis too close to zero. Singular conic sections are not allowed."
        )

    # compute the inclination
    i = np.arccos(h.z / h_mag)

    # The rest depends on ellipticity and inclination
    is_circular = True if ecc_mag < 1e-11 else False
    is_equatorial = False if (1e-11 * u.rad <= i <= (np.pi - 1e-11) * u.rad) else True

    # *** Case 1: Non-circular, Inclined Orbit ***
    if not is_circular and not is_equatorial:

        # Compute RAAN and fix quadrant
        raan = np.arccos(n.x / n_mag)
        if n.y < 0:
            raan = 2 * np.pi * u.rad - raan

        # Compute argument of periapsis and fix quadrant
        arg_p = np.arccos(n.dot(ecc) / (n_mag * ecc_mag))
        if ecc.z < 0:
            arg_p = 2 * np.pi * u.rad - arg_p

        # Compute true anomaly and fix quadrant
        true_an = np.arccos(ecc.dot(r) / (ecc_mag * r_mag))
        if r.dot(v) < 0:
            true_an = 2 * np.pi * u.rad - true_an

    # *** Case 2: Non-circular, Equatorial Orbit ***
    elif not is_circular and is_equatorial:
        # Set RAAN to zero
        raan = 0 * u.rad

        # Compute argument of periapsis and fix quadrant
        arg_p = np.arccos(ecc.x / ecc_mag)
        if ecc.y < 0:
            arg_p = 2 * np.pi * u.rad - arg_p

        # Compute true anomaly and fix quadrant
        true_an = np.arccos(ecc.dot(r) / (ecc_mag * r_mag))
        if r.dot(v) < 0:
            true_an = 2 * np.pi * u.rad - true_an

    # *** Case 3: Circular, Inclined Orbit ***
    elif is_circular and not is_equatorial:

        # Compute RAAN and fix quadrant
        raan = np.arccos(n.x / n_mag)
        if n.y < 0:
            raan = 2 * np.pi * u.rad - raan

        # Set argument of periapsis to zero
        arg_p = 0 * u.rad

        # Compute true anomaly and fix quadrant
        true_an = np.arccos(n.dot(r) / (n_mag * r_mag))
        if r.z < 0:
            true_an = 2 * np.pi * u.rad - true_an

    # *** Case 4: Circular, Equatorial Orbit ***
    else:
        # Set RAAN to zero
        raan = 0 * u.rad

        # Set argument of periapsis to zero
        arg_p = 0 * u.rad

        # Compute true anomaly and fix quadrant
        true_an = np.arccos(r.x / r_mag)
        if r.y < 0:
            true_an = 2 * np.pi * u.rad - true_an

    return OsculatingKeplerianOrbElems(
        epoch,
        sm_axis,
        ecc_mag,
        i,
        raan,
        arg_p,
        true_an,
        central_body,
    )


@u.quantity_input(angle_rad="angle", min_range="angle", max_range="angle")
def _force_angle_range(
    angle, min_range=0 * u.rad, max_range=2 * np.pi * u.rad
) -> u.rad:
    """
    Forces the angles into a prescribed range.

    Parameters
    ----------
    angle : Quantity
        Angle [rad]
    min_range : Quantity
        Minimum value allowable [rad]
    max_range : Quantity
        Max value allowable  [rad]

    """

    while not (min_range <= angle < max_range):
        if max_range < angle:
            angle = angle - 2 * np.pi * u.rad

        if angle < min_range:
            angle = angle + 2 * np.pi * u.rad

    return angle
