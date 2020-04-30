"""
Test coordinate transformations for the TEME, J2000 and TIRS coordinate frames.

Copyright (C) 2020 Egemen Imre

Licensed under GNU GPL v3.0. See LICENSE.rst for more info.

"""
from pytest import approx
from astropy import units as u
from astropy.coordinates import CartesianRepresentation, ITRS, CIRS, CartesianDifferential, GCRS
from astropy.time import Time

from satmad.frames.j2000 import J2000
from satmad.frames.teme import TEME
from satmad.frames.tirs import TIRS

time: Time = Time("2004-04-06T07:51:28.386009", scale="utc")

# Vallado IAU 2000 - Table 3-6
v_gcrs_true = CartesianDifferential(d_x=[-4.7432201610, 0.7905364950, 5.5337557240], unit=u.km / u.s, copy=True)
r_gcrs_true = CartesianRepresentation(x=[5102.50895290, 6123.01139910, 6378.13693380], unit=u.km, copy=True)
rv_gcrs_true = GCRS(r_gcrs_true.with_differentials(v_gcrs_true), obstime=time, representation_type="cartesian",
                    differential_type="cartesian")

# Vallado IAU 2000 - Table 3-6 J2000
v_j2000_true = CartesianDifferential(d_x=[-4.7432196000, 0.7905366000, 5.5337561900], unit=u.km / u.s, copy=True)
r_j2000_true = CartesianRepresentation(x=[5102.50960000, 6123.01152000, 6378.13630000], unit=u.km, copy=True)
rv_j2000_true = J2000(r_j2000_true.with_differentials(v_j2000_true), obstime=time, representation_type="cartesian",
                      differential_type="cartesian")
# Vallado IAU 2000 - Table 3-6 CIRS
v_cirs_true = CartesianDifferential(d_x=[-4.7453803300, 0.7903414530, 5.5319312880], unit=u.km / u.s, copy=True)
r_cirs_true = CartesianRepresentation(x=[5100.01840470, 6122.78636480, 6380.34453270], unit=u.km, copy=True)
rv_cirs_true = CIRS(r_cirs_true.with_differentials(v_cirs_true), obstime=time, representation_type="cartesian",
                    differential_type="cartesian")

# Vallado IAU 2000 - Table 3-6 ITRS
v_itrs_true = CartesianDifferential(d_x=[-3.2256365200, -2.8724514500, 5.5319244460], unit=u.km / u.s, copy=True)
r_itrs_true = CartesianRepresentation(x=[-1033.4793830, 7901.29527540, 6380.35659580], unit=u.km, copy=True)
rv_itrs_true = ITRS(r_itrs_true.with_differentials(v_itrs_true), obstime=time, representation_type="cartesian",
                    differential_type="cartesian")

# Vallado IAU 2000 - Table 3-6 TIRS
v_tirs_true = CartesianDifferential(d_x=[-3.2256327470, -2.8724425110, 5.5319312880], unit=u.km / u.s, copy=True)
r_tirs_true = CartesianRepresentation(x=[-1033.47503120, 7901.30558560, 6380.34453270], unit=u.km, copy=True)
rv_tirs_true = TIRS(r_tirs_true.with_differentials(v_tirs_true), obstime=time, representation_type="cartesian",
                    differential_type="cartesian")

# Vallado IAU 2000 - Table 3-6 TEME
v_teme_true = CartesianDifferential(d_x=[-4.7461314870, 0.7858180410, 5.5319312880], unit=u.km / u.s, copy=True)
r_teme_true = CartesianRepresentation(x=[5094.18016210, 6127.64465950, 6380.34453270], unit=u.km, copy=True)
rv_teme_true = TEME(r_teme_true.with_differentials(v_teme_true), obstime=time, representation_type="cartesian",
                    differential_type="cartesian")


def pos_err(rv_test, rv_true):
    """
    Computes positional error between two vectors defined in a frame.

    Parameters
    ----------
    rv_test : State with position vector under test
    rv_true : State with true position vector

    Returns
    -------
    Quantity
        3D position difference

    """
    r_diff = rv_test.cartesian.without_differentials() - rv_true.cartesian.without_differentials()
    return r_diff.norm().to(u.mm)


def vel_err(rv_test, rv_true):
    """
    Computes velocity error between two vectors defined in a frame.

    Parameters
    ----------
    rv_test : State with velocity vector under test
    rv_true : State with true velocity vector

    Returns
    -------
    Quantity
        3D velocity difference

    """
    v_diff = rv_test.velocity - rv_true.velocity
    return v_diff.norm().to(u.mm / u.s)


def test_itrs_roundtrip():
    # test_frame = "ITRS"
    allowable_pos_diff = 1.5e-6 * u.mm
    allowable_vel_diff = 2.2e-9 * u.mm / u.s

    rv_teme = rv_itrs_true.transform_to(TEME(obstime=time))
    rv_itrs_from_teme = rv_teme.transform_to(ITRS(obstime=time))

    r_diff = pos_err(rv_itrs_from_teme, rv_itrs_true)
    v_diff = vel_err(rv_itrs_from_teme, rv_itrs_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0


def test_j2000_roundtrip():
    # test_frame = "J2000"
    allowable_pos_diff = 1.5e-6 * u.mm
    allowable_vel_diff = 1.0e-9 * u.mm / u.s

    rv_gcrs = rv_j2000_true.transform_to(GCRS(obstime=time))
    rv_j2000_from_gcrs = rv_gcrs.transform_to(J2000(obstime=time))

    r_diff = pos_err(rv_j2000_from_gcrs, rv_j2000_true)
    v_diff = vel_err(rv_j2000_from_gcrs, rv_j2000_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0


def test_j2000_to_gcrs():
    # test_frame = "GCRS"
    allowable_pos_diff = 800 * u.mm
    allowable_vel_diff = 0.36 * u.mm / u.s

    rv_gcrs = rv_j2000_true.transform_to(GCRS(obstime=time))

    r_diff = pos_err(rv_gcrs, rv_gcrs_true)
    v_diff = vel_err(rv_gcrs, rv_gcrs_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0


def test_teme_to_tirs():
    # test_frame = "TIRS"
    allowable_pos_diff = 300 * u.mm
    allowable_vel_diff = 0.18 * u.mm / u.s

    rv_tirs = rv_teme_true.transform_to(TIRS(obstime=time))

    r_diff = pos_err(rv_tirs, rv_tirs_true)
    v_diff = vel_err(rv_tirs, rv_tirs_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0


def test_itrs_to_teme():
    # test_frame = "TEME"
    allowable_pos_diff = 300 * u.mm
    allowable_vel_diff = 0.21 * u.mm / u.s

    rv_teme = rv_itrs_true.transform_to(TEME(obstime=time))

    r_diff = pos_err(rv_teme, rv_teme_true)
    v_diff = vel_err(rv_teme, rv_teme_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0


def test_itrs_to_tirs():
    # test_frame = "TIRS"
    allowable_pos_diff = 60 * u.mm
    allowable_vel_diff = 0.05 * u.mm / u.s

    rv_tirs = rv_itrs_true.transform_to(TIRS(obstime=time))

    r_diff = pos_err(rv_tirs, rv_tirs_true)
    v_diff = vel_err(rv_tirs, rv_tirs_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0


def test_tirs_to_itrs():
    # test_frame = "ITRS"
    allowable_pos_diff = 60 * u.mm
    allowable_vel_diff = 0.05 * u.mm / u.s

    rv_itrs = rv_tirs_true.transform_to(ITRS(obstime=time))

    r_diff = pos_err(rv_itrs, rv_itrs_true)
    v_diff = vel_err(rv_itrs, rv_itrs_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert approx(r_diff.value, abs=allowable_pos_diff.value) == 0.0
    assert approx(v_diff.value, abs=allowable_vel_diff.value) == 0.0
