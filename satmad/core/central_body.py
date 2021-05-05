# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Definition of a Central Body, which is a celestial body around which
an object can rotate.

"""
from astropy import units as u


class CelestialBodyEllipsoid:
    """Defines an ellipsoid, the shape of most celestial bodies.

    Parameters
    ----------
    name : str
        Name of the Ellipsoid
    re : Quantity
        Equatorial Radius or Semimajor Axis of the Ellipsoid. (km)
    inv_f : Quantity
        Inverse flattening (:math:`1/f`). (dimensionless)
    j2 : Quantity
        Dynamical form factor
    om : Quantity
        Nominal mean angular velocity of the celestial body

    """

    j2 = None

    #     mu : Quantity
    #         `GM` value (Gravitational constant) of the Celestial body
    # mu = None

    om = None

    @u.quantity_input(
        re=u.km,
        inv_f=u.dimensionless_unscaled,
        j2=u.dimensionless_unscaled,
        mu=u.m ** 3 / u.s ** 2,
        om=u.rad / u.s,
    )
    def __init__(self, name, re, inv_f, **kwargs):
        self.name = name

        self.re = re
        self.inv_f = inv_f

        self.j2 = kwargs.get("j2")
        # self.mu = kwargs.get("mu")

        self.om = kwargs.get("om")


class CelestialBody:
    """
    Container class for Celestial Bodies and related information.

    Parameters
    ----------
    name : str
        Name of the body
    info : str
        Info on the constants etc.
    mu : Quantity
        `GM` value (Gravitational constant) of the Celestial body
    inert_coord : `~astropy.coordinates.BaseRepresentation`
        Default Inertial Coordinate to run the propagations (e.g. `GCRS` for Earth)
    body_fixed_coord : `~astropy.coordinates.BaseRepresentation`
        Default Central Body Fixed Coordinate to run the propagations (e.g. `ITRS` for
        Earth)
    """

    @u.quantity_input(mu=u.km ** 3 / u.s ** 2)
    def __init__(self, name, info, mu, **kwargs):
        self.name = name
        self.info = info

        self.mu = mu

        self.ellipsoid = kwargs.get("ellipsoid")

        self.inert_coord = kwargs.get("inert_coord")
        self.body_fixed_coord = kwargs.get("body_fixed_coord")

    def __str__(self) -> str:
        out = f"Name: {self.name}\n"

        out += f"{self.info}\n"

        out += f"mu: {self.mu}\n"

        out += f"Default inertial coord: {self.inert_coord.name}\n"
        out += f"Default body fixed coord: {self.body_fixed_coord.name}\n"

        return out

    @property
    def inert_coord_frame(self):
        """Gets the underlying inertial coordinate frame class (e.g. GCRS for Earth)
        or None if no such class exists."""
        return self.inert_coord

    @property
    def body_fixed_coord_frame(self):
        """Gets the underlying body fixed coordinate frame class (e.g. ITRS for Earth)
        or None if no such class exists."""
        return self.body_fixed_coord
