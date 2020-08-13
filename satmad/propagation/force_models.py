# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Force models generally for use with the numerical propagators.

"""
import numpy as np


def two_body_accel(r, mu):
    r"""
    Two-body acceleration.

    The input vectors as well as the output acceleration is in in the inertial frame
    of the central body. For the Earth this is GCRS.

    The acceleration value is given by:

    .. math::

        \ddot{\vec{r}} = - \dfrac{\mu}{r^3} \vec{r}

    where :math:`\ddot{\vec{r}}` is the inertial acceleration,
    :math:`\vec{r}` is the position vector (with :math:`r` its norm).
    :math:`\mu` is equal to the constant :math:`GM` (Gravitational Constant times
    the Mass of the main attracting body.)

    .. note::
        This method is used in the Scipy ODE solver, therefore it does not accept
        `Quantity` inputs. Care must be taken with the units.


    Parameters
    ----------
    r : ndarray
        inertial position vector with 3 elements (km)
    mu : float
        GM value (:math:`km^3 / s^2`)

    Returns
    -------
    float
        inertial acceleration (:math:`km / s^2`)

    """
    r_norm = np.linalg.norm(r)
    a = -mu / (r_norm ** 3) * r
    return a


def two_body_energy(r, v, mu):
    r"""Two-body specific energy.

    The input vectors are in in the inertial frame of the central body.
    For the Earth this is GCRS.

    The specific energy value is given by:

    .. math::

        \varepsilon = \dfrac{v^2}{2} - \dfrac{\mu}{r}

    with :math:`v` is the norm of the velocity vector, :math:`r` is the norm of the
    position vector and :math:`\mu` is equal to the constant :math:`GM`
    (Gravitational Constant times the Mass of the main attracting body).

    This is the formula for the specific energy, or the sum of
    potential and kinetic energy *per unit mass*.

    Parameters
    ----------
    r : ndarray or Quantity
        inertial position vector with 3 elements (km)
    v : ndarray or Quantity
        inertial position vector with 3 elements (km/s)
    mu : float or Quantity
        GM value (:math:`km^3 / s^2`)

    Returns
    -------
    sp_energy : float or Quantity
        specific energy (:math:`km^2 / s^2`)
    """
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    sp_energy = 0.5 * v_norm ** 2 - mu / r_norm

    return sp_energy
