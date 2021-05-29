# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Test coordinate transformations for the coordinate frames defined by SatMAD.

"""
import pytest
from astropy import units as u
from astropy.coordinates import (
    CIRS,
    GCRS,
    HCRS,
    ICRS,
    ITRS,
    TEME,
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord,
)
from astropy.time import Time

from satmad.coordinates.frames import (
    J2000,
    TIRS,
    CelestialBodyCRS,
    CelestialBodyFixed,
    CelestialBodyJ2000Equatorial,
    CelestialBodyTODEquatorial,
    MarsCRS,
    MarsJ2000Equatorial,
    MarsTODEquatorial,
    init_pvt,
)
from satmad.tests.common_test_funcs import pos_err, vel_err

time: Time = Time("2004-04-06T07:51:28.386009", scale="utc")

# Vallado IAU 2000 - Table 3-6
v_gcrs_true = CartesianDifferential(
    [-4.7432201610, 0.7905364950, 5.5337557240], unit=u.km / u.s
)
r_gcrs_true = CartesianRepresentation(
    [5102.50895290, 6123.01139910, 6378.13693380], unit=u.km
)
rv_gcrs_true = GCRS(
    r_gcrs_true.with_differentials(v_gcrs_true),
    obstime=time,
    representation_type="cartesian",
    differential_type="cartesian",
)

# Vallado IAU 2000 - Table 3-6 J2000
v_j2000_true = CartesianDifferential(
    [-4.7432196000, 0.7905366000, 5.5337561900], unit=u.km / u.s
)
r_j2000_true = CartesianRepresentation(
    [5102.50960000, 6123.01152000, 6378.13630000], unit=u.km
)
rv_j2000_true = J2000(
    r_j2000_true.with_differentials(v_j2000_true),
    obstime=time,
    representation_type="cartesian",
    differential_type="cartesian",
)
# Vallado IAU 2000 - Table 3-6 CIRS
v_cirs_true = CartesianDifferential(
    [-4.7453803300, 0.7903414530, 5.5319312880], unit=u.km / u.s
)
r_cirs_true = CartesianRepresentation(
    [5100.01840470, 6122.78636480, 6380.34453270], unit=u.km
)
rv_cirs_true = CIRS(
    r_cirs_true.with_differentials(v_cirs_true),
    obstime=time,
    representation_type="cartesian",
    differential_type="cartesian",
)

# Vallado IAU 2000 - Table 3-6 ITRS
v_itrs_true = CartesianDifferential(
    [-3.2256365200, -2.8724514500, 5.5319244460], unit=u.km / u.s
)
r_itrs_true = CartesianRepresentation(
    [-1033.4793830, 7901.29527540, 6380.35659580], unit=u.km
)
rv_itrs_true = ITRS(
    r_itrs_true.with_differentials(v_itrs_true),
    obstime=time,
    representation_type="cartesian",
    differential_type="cartesian",
)

# Vallado IAU 2000 - Table 3-6 TIRS
v_tirs_true = CartesianDifferential(
    [-3.2256327470, -2.8724425110, 5.5319312880], unit=u.km / u.s
)
r_tirs_true = CartesianRepresentation(
    [-1033.47503120, 7901.30558560, 6380.34453270], unit=u.km
)
rv_tirs_true = TIRS(
    r_tirs_true.with_differentials(v_tirs_true),
    obstime=time,
    representation_type="cartesian",
    differential_type="cartesian",
)

# Vallado IAU 2000 - Table 3-6 TEME
v_teme_true = CartesianDifferential(
    [-4.7461314870, 0.7858180410, 5.5319312880], unit=u.km / u.s
)
r_teme_true = CartesianRepresentation(
    [5094.18016210, 6127.64465950, 6380.34453270], unit=u.km
)
rv_teme_true = TEME(
    r_teme_true.with_differentials(v_teme_true),
    obstime=time,
    representation_type="cartesian",
    differential_type="cartesian",
)


def _init_orbit_mars():
    """Initialises the test orbit in a Martian orbit."""
    obstime = Time("2020-01-10T11:30:00.0003", scale="tdb")

    v_mars = CartesianDifferential(
        [3.057934230169251e-01, -3.068219328642159e00, -1.397389382391015e00],
        unit=u.km / u.s,
    )
    r_mars = CartesianRepresentation(
        [1.786125452441044e03, -1.191541086863880e03, 3.003639539248147e03], unit=u.km
    )

    rv_mars = SkyCoord(
        r_mars.with_differentials(v_mars),
        obstime=obstime,
        frame=MarsCRS,
        representation_type="cartesian",
        differential_type="cartesian",
    )

    return rv_mars


# ********** Cartesian init testing **********


def test_init_cart():
    """Tests the cartesian initialisation convenience method."""

    rv_gcrs_true_sky = SkyCoord(rv_gcrs_true)

    # easy init
    rv_gcrs_0 = init_pvt(GCRS, time, r_gcrs_true, v_gcrs_true)

    # init with str frame name (should be lowercase letters)
    rv_gcrs_1 = init_pvt("GCRS", time, r_gcrs_true, v_gcrs_true)

    # init without velocity
    rv_gcrs_2 = init_pvt(GCRS, time, r_gcrs_true, copy=False)

    assert rv_gcrs_0 == rv_gcrs_true_sky
    assert rv_gcrs_1 == rv_gcrs_true_sky
    assert rv_gcrs_2 == SkyCoord(
        GCRS(
            r_gcrs_true,
            obstime=time,
            representation_type="cartesian",
            differential_type="cartesian",
        )
    )


def test_init_cart_wrong_frame():
    """Tests cartesian init with incorrect frame name."""
    with pytest.raises(ValueError):
        init_pvt("GCRSx", time, r_gcrs_true, v_gcrs_true)


def test_init_cart_wrong_time():
    """Tests cartesian init with incorrect frame name."""
    with pytest.raises(TypeError):
        # noinspection PyTypeChecker
        init_pvt("GCRS", "2004-04-06T07:51:28.386009", r_gcrs_true, v_gcrs_true)


# ********** Functional testing **********


def test_tirs_to_teme_no_vel():
    """Check whether coord transform without velocity is possible."""
    rv_tirs_no_vel = TIRS(r_tirs_true, obstime=time, representation_type="cartesian")
    rv_tirs_no_vel.transform_to(TEME(obstime=time))


def test_teme_to_tirs_no_vel():
    """Check whether coord transform without velocity is possible."""
    rv_teme_no_vel = TEME(r_teme_true, obstime=time, representation_type="cartesian")
    rv_teme_no_vel.transform_to(TIRS(obstime=time))


def test_itrs_roundtrip():
    """Check whether transforming to a coord and then
    transforming back yields the same output."""
    # test_frame = "ITRS"
    allowable_pos_diff = 1.5e-6 * u.mm
    allowable_vel_diff = 2.8e-6 * u.mm / u.s

    rv_teme = rv_itrs_true.transform_to(TEME(obstime=time))
    rv_itrs_from_teme = rv_teme.transform_to(ITRS(obstime=time))

    r_diff = pos_err(rv_itrs_from_teme, rv_itrs_true)
    v_diff = vel_err(rv_itrs_from_teme, rv_itrs_true)

    # print(f"r {rv_itrs_true.name} diff      :  {r_diff}")
    # print(f"v {rv_itrs_true.name} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


def test_j2000_roundtrip():
    """Check whether transforming to a coord and then
    transforming back yields the same output."""
    # test_frame = "J2000"
    allowable_pos_diff = 1.5e-6 * u.mm
    allowable_vel_diff = 1.0e-9 * u.mm / u.s

    rv_gcrs = rv_j2000_true.transform_to(GCRS(obstime=time))
    rv_j2000_from_gcrs = rv_gcrs.transform_to(J2000(obstime=time))

    r_diff = pos_err(rv_j2000_from_gcrs, rv_j2000_true)
    v_diff = vel_err(rv_j2000_from_gcrs, rv_j2000_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


# ********** Performance testing **********


def test_j2000_to_gcrs():
    """Check the coordinate transform accuracy."""
    # test_frame = "GCRS"
    allowable_pos_diff = 800 * u.mm
    allowable_vel_diff = 0.36 * u.mm / u.s

    rv_gcrs = rv_j2000_true.transform_to(GCRS(obstime=time))

    r_diff = pos_err(rv_gcrs, rv_gcrs_true)
    v_diff = vel_err(rv_gcrs, rv_gcrs_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


def test_teme_to_tirs():
    """Check the coordinate transform accuracy."""
    # test_frame = "TIRS"
    allowable_pos_diff = 300 * u.mm
    allowable_vel_diff = 0.20 * u.mm / u.s

    rv_tirs = rv_teme_true.transform_to(TIRS(obstime=time))

    r_diff = pos_err(rv_tirs, rv_tirs_true)
    v_diff = vel_err(rv_tirs, rv_tirs_true)

    # print(f"r {rv_tirs.name} diff      :  {r_diff}")
    # print(f"v {rv_tirs.name} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


def test_itrs_to_teme():
    """Check the coordinate transform accuracy."""
    # test_frame = "TEME"
    allowable_pos_diff = 300 * u.mm
    allowable_vel_diff = 0.23 * u.mm / u.s

    rv_teme = rv_itrs_true.transform_to(TEME(obstime=time))

    r_diff = pos_err(rv_teme, rv_teme_true)
    v_diff = vel_err(rv_teme, rv_teme_true)

    # print(f"r {rv_teme.name} diff      :  {r_diff}")
    # print(f"v {rv_teme.name} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


def test_itrs_to_tirs():
    """Check the coordinate transform accuracy."""
    # test_frame = "TIRS"
    allowable_pos_diff = 60 * u.mm
    allowable_vel_diff = 0.05 * u.mm / u.s

    rv_tirs = rv_itrs_true.transform_to(TIRS(obstime=time))

    r_diff = pos_err(rv_tirs, rv_tirs_true)
    v_diff = vel_err(rv_tirs, rv_tirs_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


def test_tirs_to_itrs():
    """Check the coordinate transform accuracy."""
    # test_frame = "ITRS"
    allowable_pos_diff = 60 * u.mm
    allowable_vel_diff = 0.05 * u.mm / u.s

    rv_itrs = rv_tirs_true.transform_to(ITRS(obstime=time))

    r_diff = pos_err(rv_itrs, rv_itrs_true)
    v_diff = vel_err(rv_itrs, rv_itrs_true)

    # print(f"r {test_frame} diff      :  {r_diff}")
    # print(f"v {test_frame} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


# ********************** Test CRS *************************


class _EarthCRS(CelestialBodyCRS):
    body_name = "EARTH"
    ephemeris_type = "builtin"


class _SunCRS(CelestialBodyCRS):
    body_name = "SUN"
    ephemeris_type = "builtin"


def test_cb_crs_to_icrs_earth():
    """Check the basic quasi-round-trip accuracy for the generic Central Body Celestial
    Reference System, using Earth (equal to GCRS). Converts GCRS to ICRS, then ICRS to
    Earth CRS. However there is a difference, possibly due to the difference between
    the two transformations (`FunctionTransformWithFiniteDifference` vs
    `AffineTransform`)"""
    allowable_pos_diff = 800 * u.m
    allowable_vel_diff = 0.510 * u.m / u.s

    rv_icrs = SkyCoord(rv_gcrs_true).transform_to(ICRS)

    rv_gcrs = rv_icrs.transform_to(_EarthCRS)

    r_diff = pos_err(rv_gcrs, rv_gcrs_true)
    v_diff = vel_err(rv_gcrs, rv_gcrs_true)

    # print(f"r {rv_gcrs.name} diff      :  {r_diff}")
    # print(f"v {rv_gcrs.name} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


def test_cb_crs_to_icrs_sun():
    """Check the basic quasi-round-trip accuracy for the generic Central Body Celestial
    Reference System, using Sun (equal to HCRS). Converts Sun CRS to ICRS, then ICRS to
    HCRS."""
    allowable_pos_diff = 5e-4 * u.mm
    allowable_vel_diff = 0.0001 * u.mm / u.s

    rv_suncrs_true = init_pvt(
        _SunCRS, time, r_gcrs_true.with_differentials(v_gcrs_true)
    )

    rv_icrs = rv_suncrs_true.transform_to(ICRS)

    rv_hcrs = rv_icrs.transform_to(HCRS)

    r_diff = pos_err(rv_hcrs, rv_suncrs_true)
    v_diff = vel_err(rv_hcrs, rv_suncrs_true)

    # print(f"r {rv_hcrs.name} diff      :  {r_diff}")
    # print(f"v {rv_hcrs.name} diff      :  {v_diff}")

    assert r_diff < allowable_pos_diff
    assert v_diff < allowable_vel_diff


# ********************** Test TOD and J2000 Equatorial *************************


def test_equatorial_round_trip_mars():
    """Round trip testing between CRS and TOD & J2000 Equatorial."""

    allowable_pos_diff = 2e-6 * u.mm
    allowable_vel_diff = 8e-10 * u.mm / u.s

    # init CRS
    pvt_crs = _init_orbit_mars()

    # Convert to TOD Equatorial
    pvt_tod_eq = pvt_crs.transform_to(MarsTODEquatorial)
    # Convert to J2000 Equatorial
    pvt_j2000_eq = pvt_crs.transform_to(MarsJ2000Equatorial)

    # Convert J2000 back to CRS
    assert pos_err(pvt_j2000_eq.transform_to(MarsCRS), pvt_crs) < allowable_pos_diff
    assert vel_err(pvt_j2000_eq.transform_to(MarsCRS), pvt_crs) < allowable_vel_diff

    # print(pos_err(pvt_tod_eq.transform_to(MarsCRS), pvt_crs))
    # print(vel_err(pvt_tod_eq.transform_to(MarsCRS), pvt_crs))

    # Convert TOD back to CRS
    assert pos_err(pvt_tod_eq.transform_to(MarsCRS), pvt_crs) < allowable_pos_diff
    assert vel_err(pvt_tod_eq.transform_to(MarsCRS), pvt_crs) < allowable_vel_diff


class _MarsIAU2000TODEquatorial(CelestialBodyTODEquatorial):
    """
    A coordinate frame representing the True-of-Date Equatorial System of Mars.
    """

    body_name = "Mars_IAU_2000"
    cb_crs = MarsCRS


class _MarsIAU2000J2000Equatorial(CelestialBodyJ2000Equatorial):
    """
    A coordinate frame representing the Equatorial System of Mars at J2000 Epoch.
    """

    body_name = "Mars_IAU_2000"
    cb_crs = MarsCRS


class _MarsIAU2000BodyFixed(CelestialBodyFixed):
    """
    A coordinate frame representing the Body Fixed System of Mars.
    """

    body_name = "Mars_IAU_2000"
    cb_crs = MarsCRS


def test_equatorial_mars():
    """NASA HORIZONS/GMAT testing between CRS and TOD & J2000 Equatorial
    as well as Body Fixed (GMAT).

    GMAT and HORIZONS use a different model in the IAU polar rotation parameters
    (IAU 2000). Therefore test setup uses this version but not IAU 2015 versions.

    Satellite orbit is taken from Mars Reconnaissance Orbiter,
    using the NASA Horizons Web Interface (https://ssd.jpl.nasa.gov/horizons.cgi)."""

    obstime = Time("2020-01-10T11:30:00.0003", scale="tdb")

    v_eq_horz = CartesianDifferential(
        [-2.061955207079347, -2.671302888480574, 0.269551765393186],
        unit=u.km / u.s,
    )
    r_eq_horz = CartesianRepresentation(
        [322.259677663235, 120.3781499221643, 3676.074343158492], unit=u.km
    )
    rv_eq_horz = init_pvt(_MarsIAU2000TODEquatorial, obstime, r_eq_horz, v_eq_horz)

    v_inert_horz = CartesianDifferential(
        [-2.062805178080304, -2.670750348642094, 0.2685217398016297],
        unit=u.km / u.s,
    )
    r_inert_horz = CartesianRepresentation(
        [321.472501283011, 119.500545946684, 3676.171898278387], unit=u.km
    )
    rv_inert_horz = init_pvt(
        _MarsIAU2000J2000Equatorial, obstime, r_inert_horz, v_inert_horz
    )

    v_body_fixed_gmat = CartesianDifferential(
        [0.622659762445902, 3.329753468204519, 0.2695517652769852],
        unit=u.km / u.s,
    )
    r_body_fixed_gmat = CartesianRepresentation(
        [-233.7124768823157, -252.4295484193542, 3676.074343157078], unit=u.km
    )
    rv_body_fixed_gmat = init_pvt(
        _MarsIAU2000BodyFixed, obstime, r_body_fixed_gmat, v_body_fixed_gmat
    )

    allowable_pos_diff = 800 * u.mm
    allowable_vel_diff = 0.35 * u.mm / u.s

    # init CRS - modify frame name to point at test coords
    pvt_crs = _init_orbit_mars()
    pvt_crs.frame.body_name = "Mars_IAU_2000"

    # Convert to TOD Equatorial
    pvt_tod_eq = pvt_crs.transform_to(_MarsIAU2000TODEquatorial)
    # Convert to J2000 Equatorial
    pvt_j2000_eq = pvt_crs.transform_to(_MarsIAU2000J2000Equatorial)
    # Convert to Body Fixed
    pvt_body_fixed = pvt_crs.transform_to(_MarsIAU2000BodyFixed)

    print(pvt_crs)
    print(pvt_tod_eq)
    print(pvt_j2000_eq)
    print(pvt_body_fixed)

    print(pos_err(pvt_body_fixed, rv_body_fixed_gmat))
    print(vel_err(pvt_body_fixed, rv_body_fixed_gmat))

    # Diff between TOD Eq values
    assert pos_err(pvt_tod_eq, rv_eq_horz) < allowable_pos_diff
    assert vel_err(pvt_tod_eq, rv_eq_horz) < allowable_vel_diff

    # Diff between J2000 Eq (Inertial) values
    assert pos_err(pvt_j2000_eq, rv_inert_horz) < allowable_pos_diff
    assert vel_err(pvt_j2000_eq, rv_inert_horz) < allowable_vel_diff

    # Diff between Body Fixed values
    assert pos_err(pvt_body_fixed, rv_body_fixed_gmat) < allowable_pos_diff
    assert vel_err(pvt_body_fixed, rv_body_fixed_gmat) < allowable_vel_diff
