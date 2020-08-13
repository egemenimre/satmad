# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Defines a Celestial Body, which acts as the central repository of information for
the planets or similar Central Bodies.

"""
from astropy import constants as const
from astropy import units as u


class CelestialBody:
    """
    Container for the various Celestial Bodies.

    Parameters
    ----------
    name : str
        Name of the body
    info : str
        Info on the constants etc.
    mu : Quantity
        `GM` value of the Celestial body
    inert_coord : str
        Default Inertial Coordinate to run the propagations (e.g. "gcrs" for Earth)
    body_fixed_coord : str
        Default Central Body Fixed Coordinate to run the propagations (e.g. "itrs" for
        Earth)
    """

    @u.quantity_input(mu=u.km ** 3 / u.s ** 2)
    def __init__(self, name, info, mu, **kwargs):
        self.name = name
        self.info = info

        self.mu = mu

        self.inert_coord = kwargs.get("inert_coord")
        self.body_fixed_coord = kwargs.get("body_fixed_coord")

    def __str__(self) -> str:
        out = f"Name: {self.name}\n"

        out += f"{self.info}\n"

        out += f"mu: {self.mu}\n"

        out += f"Default inertial coord: {self.inert_coord}\n"
        out += f"Default body fixed coord: {self.body_fixed_coord}\n"

        return out


EARTH = CelestialBody(
    "Earth",
    "mu: Astropy Default",
    const.GM_earth.to(u.km ** 3 / u.s ** 2),
    inert_coord="gcrs",
    body_fixed_coord="itrs",
)
SUN = CelestialBody("Sun", "mu: Astropy Default", const.GM_sun.to(u.km ** 3 / u.s ** 2))
