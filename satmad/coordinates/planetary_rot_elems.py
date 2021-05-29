# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Planetary rotational parameters, defined relative to their mean axis of rotation
and body-dependent definitions of longitude.

All formulas (with the exception of the Moon) are from:
Archinal, B.A., Acton, C.H., A’Hearn, M.F. et al. Report of the IAU Working Group on
Cartographic Coordinates and Rotational Elements: 2015.
Celest Mech Dyn Astr 130, 22 (2018). https://doi.org/10.1007/s10569-017-9805-5

Moon formulas are from the 2009 edition as the 2015 version no longer lists the
explicit lunar formulae - they refer to DE430 or the older DE421 values.
Archinal, B.A., A’Hearn, M.F., Bowell, E. et al. Report of the IAU Working Group on
Cartographic Coordinates and Rotational Elements: 2009.
Celest Mech Dyn Astr 109, 101–135 (2011). https://doi.org/10.1007/s10569-010-9320-4
"""
import numpy as np
from astropy import units as u
from astropy.time import Time


def celestial_body_rot_params(epoch, body_name):
    """Retrieves the rotational parameters for planets and satellites for the
    requested epoch.

    These parameters link the Celestial Body ICRS to the equatorial inertial or
    body-fixed coordinates.

    Parameters
    ----------
    epoch : ~astropy.time.time
        Epoch at which rotational parameters are evaluated, defaults to J2000 epoch
    body_name : str
        Name of the celestial body (planet or satellite)

    Returns
    -------
    ra, dec, w : Tuple[~astropy.units.Quantity]
        Right ascension and declination of celestial body north pole, and prime meridian angle

    Raises
    ------
    NotImplementedError
        Rotational parameters not available for the requested body name
    """
    if epoch:
        # Interval from the standard epoch, in days
        d = (epoch.tdb - Time("J2000", scale="tdb")).to(u.day).value
    else:
        d = 0
    # Interval in Julian centuries (36,525 days) from the standard epoch
    t = d / 36525.0

    # generate function name
    func_name = body_name.lower() + "_rot_params"

    try:
        # run function
        ra, dec, w = eval(func_name)(t, d)
    except NameError:
        raise NotImplementedError(
            f"Requested Planet/Body Rotational parameters are not defined for {body_name}"
        )

    return ra * u.deg, dec * u.deg, w * u.deg


def moon_rot_params(t, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for the Moon.

    These parameters are useful to convert local ICRS frame to the Principal Axis (PA)
    frame of the Moon. The Mean Earth/rotation axis (ME) has an additional fixed offset
    with respect to PA.

    Parameters
    ----------
    t : float
        Julian Centuries since J2000 epoch
    d : float
        Days since J2000 epoch
    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    e_1 = 125.045 - 0.0529921 * d
    e_2 = 250.089 - 0.1059842 * d
    e_3 = 260.008 + 13.0120009 * d
    e_4 = 176.625 + 13.3407154 * d
    e_5 = 357.529 + 0.9856003 * d
    e_6 = 311.589 + 26.4057084 * d
    e_7 = 134.963 + 13.0649930 * d
    e_8 = 276.617 + 0.3287146 * d
    e_9 = 34.226 + 1.7484877 * d
    e_10 = 15.134 - 0.1589763 * d
    e_11 = 119.743 + 0.0036096 * d
    e_12 = 239.961 + 0.1643573 * d
    e_13 = 25.053 + 12.9590088 * d

    ra = (
        269.9949
        + 0.0031 * t
        - 3.8787 * np.sin(e_1)
        - 0.1204 * np.sin(e_2)
        + 0.0700 * np.sin(e_3)
        - 0.0172 * np.sin(e_4)
        + 0.0072 * np.sin(e_6)
        - 0.0052 * np.sin(e_10)
        + 0.0043 * np.sin(e_13)
    )
    dec = (
        66.5392
        + 0.0130 * t
        + 1.5419 * np.cos(e_1)
        + 0.0239 * np.cos(e_2)
        - 0.0278 * np.cos(e_3)
        + 0.0068 * np.cos(e_4)
        - 0.0029 * np.cos(e_6)
        + 0.0009 * np.cos(e_7)
        + 0.0008 * np.cos(e_10)
        - 0.0009 * np.cos(e_13)
    )
    w = (
        38.3213
        + 13.17635815 * d
        - 1.4e-12 * d ** 2
        + 3.5610 * np.sin(e_1)
        + 0.1208 * np.sin(e_2)
        - 0.0642 * np.sin(e_3)
        + 0.0158 * np.sin(e_4)
        + 0.0252 * np.sin(e_5)
        - 0.0066 * np.sin(e_6)
        - 0.0047 * np.sin(e_7)
        - 0.0046 * np.sin(e_8)
        + 0.0028 * np.sin(e_9)
        + 0.0052 * np.sin(e_10)
        + 0.0040 * np.sin(e_11)
        + 0.0019 * np.sin(e_12)
        - 0.0044 * np.sin(e_13)
    )

    return ra, dec, w


def sun_rot_params(_, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for the Sun.

    Parameters
    ----------
    _
        Empty parameter
    d : float
        Days since J2000 epoch

    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    ra = 286.13
    dec = 63.87
    w = 84.176 + 14.1844000 * d

    return ra, dec, w


def mercury_rot_params(t, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for Mercury.

    Parameters
    ----------
    t : float
        Julian Centuries since J2000 epoch
    d : float
        Days since J2000 epoch
    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    m_1 = 174.7910857 + 4.092335 * d
    m_2 = 349.5821714 + 8.184670 * d
    m_3 = 164.3732571 + 12.277005 * d
    m_4 = 339.1643429 + 16.369340 * d
    m_5 = 153.9554286 + 20.461675 * d

    ra = 281.0103 - 0.0328 * t
    dec = 61.45 - 0.005 * t
    w = (329.5988 + 6.1385108 * d) + (
        0.01067257 * np.sin(m_1)
        - 0.00112309 * np.sin(m_2)
        - 0.00011040 * np.sin(m_3)
        - 0.00002539 * np.sin(m_4)
        - 0.00000571 * np.sin(m_5)
    )

    return ra, dec, w


def venus_rot_params(_, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for Venus.

    Parameters
    ----------
    _
        Empty parameter
    d : float
        Days since J2000 epoch
    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    ra = 272.76
    dec = 67.16
    w = 160.20 - 1.4813688 * d

    return ra, dec, w


def mars_rot_params(t, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for Mars.

    Parameters
    ----------
    t : float
        Julian Centuries since J2000 epoch
    d : float
        Days since J2000 epoch

    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    ra = (
        317.269202
        - 0.10927547 * t
        + 0.000068 * np.sin(198.991226 + 19139.4819985 * t)
        + 0.000238 * np.sin(226.292679 + 38280.8511281 * t)
        + 0.000052 * np.sin(249.663391 + 57420.7251593 * t)
        + 0.000009 * np.sin(266.183510 + 76560.6367950 * t)
        + 0.419057 * np.sin(79.398797 + 0.5042615 * t)
    )

    dec = (
        54.432516
        - 0.05827105 * t
        + 0.000051 * np.cos(122.433576 + 19139.9407476 * t)
        + 0.000141 * np.cos(43.058401 + 38280.8753272 * t)
        + 0.000031 * np.cos(57.663379 + 57420.7517205 * t)
        + 0.000005 * np.cos(79.476401 + 76560.6495004 * t)
        + 1.591274 * np.cos(166.325722 + 0.5042615 * t)
    )

    w = (
        176.049863
        + 350.891982443297 * d
        + 0.000145 * np.sin(129.071773 + 19140.0328244 * t)
        + 0.000157 * np.sin(36.352167 + 38281.0473591 * t)
        + 0.000040 * np.sin(56.668646 + 57420.9295360 * t)
        + 0.000001 * np.sin(67.364003 + 76560.2552215 * t)
        + 0.000001 * np.sin(104.792680 + 95700.4387578 * t)
        + 0.584542 * np.sin(95.391654 + 0.5042615 * t)
    )

    return ra, dec, w


def mars_iau_2000_rot_params(t, d):
    """This is used for test purposes only. Parameters are taken from the
    2000 version of the IAU Report."""
    ra = 317.68143 - 0.1061 * t
    dec = 52.88650 - 0.0609 * t
    w = 176.630 + 350.89198226 * d
    return ra, dec, w


def jupiter_rot_params(t, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for Jupiter.

    Parameters
    ----------
    t : float
        Julian Centuries since J2000 epoch
    d : float
        Days since J2000 epoch
    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    j_a = 99.360714 + 4850.4046 * t
    j_b = 175.895369 + 1191.9605 * t
    j_c = 300.323162 + 262.5475 * t
    j_d = 114.012305 + 6070.2476 * t
    j_e = 49.511251 + 64.3000 * t

    ra = (
        268.056595
        - 0.006499 * t
        + 0.000117 * np.sin(j_a)
        + 0.000938 * np.sin(j_b)
        + 0.001432 * np.sin(j_c)
        + 0.000030 * np.sin(j_d)
        + 0.002150 * np.sin(j_e)
    )
    dec = (
        64.495303
        + 0.002413 * t
        + 0.000050 * np.cos(j_a)
        + 0.000404 * np.cos(j_b)
        + 0.000617 * np.cos(j_c)
        - 0.000013 * np.cos(j_d)
        + 0.000926 * np.cos(j_e)
    )
    w = 284.95 + 870.536 * d

    return ra, dec, w


def saturn_rot_params(t, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for Saturn.

    Parameters
    ----------
    t : float
        Julian Centuries since J2000 epoch
    d : float
        Days since J2000 epoch

    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    ra = 40.589 - 0.036 * t
    dec = 83.537 - 0.004 * t
    w = 38.90 + 810.7939024 * d

    return ra, dec, w


def uranus_rot_params(_, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for Uranus.

    Parameters
    ----------
    _
        Empty parameter
    d : float
        Days since J2000 epoch

    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    ra = 257.311
    dec = -15.175
    w = 203.81 - 501.1600928 * d

    return ra, dec, w


def neptune_rot_params(t, d):
    """Computes the North Pole rotational parameters (ra, dec, prime meridian angle)
    for Neptune.

    Parameters
    ----------
    t : float
        Julian Centuries since J2000 epoch
    d : float
        Days since J2000 epoch

    Returns
    -------
    ra, dec, w : Tuple[float]
        Right ascension and declination of celestial body north pole, and prime meridian angle
    """
    n = 357.85 + 52.316 * t

    ra = 299.36 + 0.70 * np.sin(n)
    dec = 43.46 - 0.51 * np.cos(n)
    w = 249.978 + 541.1397757 * d - 0.48 * np.sin(n)

    return ra, dec, w
