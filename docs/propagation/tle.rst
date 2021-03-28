Orbit Initialisation with Two-Line Elements (TLEs)
====================================================

.. _tle-intro:

Introduction to TLEs
---------------------

A two-line element set (TLE) is a data format containing a set of TEME
(True Equator, Mean Equinox) mean orbital elements of an Earth-orbiting object
for a given point in time, called the Epoch Time.

These orbital elements are solely for use with the SGP4 propagator as the two are coupled
with the underlying analytical orbit theory. [OM1]_

See the `TLE page in Wikipedia
<https://en.wikipedia.org/wiki/Two-line_element_set>`_ or `NASA definition
<https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/SSOP_Help/tle_def.html>`_
for more information.

An example TLE is given as::

    ISS (ZARYA)
    1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
    2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537

A TLE object is usually initialised from these two lines of strings:

    >>> line1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
    >>> line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"
    >>> tle = TLE.from_tle(line1, line2, "ISS (ZARYA)")

Components of a TLE
-------------------

A TLE contains the well-known mean classical (or Keplerian) orbital elements such as
Mean Anomaly, Right Ascension of the Ascending Node or
Argument of Perigee. A semimajor axis is not directly defined, but the mean motion term
given in revolutions/day can be used to derive the semimajor axis in a straightforward fashion.
For example, a mean motion of 14.5 revolutions/day is equivalent to an orbital period of
5958.62 seconds or a mean motion (:math:`n`) of 0.00106 radians/second. Then the usual
semimajor axis (:math:`a`) can be derived from :math:`a^3 n^2=\mu`, where :math:`\mu` is
equal to the Gravitational Constant (:math:`G`) times the mass of the Earth (:math:`M`).
In this example, the semimajor axis is equal to 7079.056 km.

In addition to the usual classical orbital elements, other components of the TLE are
(adapted from the
`NASA definition <https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/SSOP_Help/tle_def.html>`_
):

* Satellite catalog number (NORAD ID): The unique string
  representation of the object in space, assigned by USSTRATCOM.
  For example, the input integer of `25544` means ISS ZARYA module

* International Designator (COSPAR ID or NSSDCA ID): The unique string
  representation of the object in space. The field in the TLE `98067A` means
  1998-067A, which means the Object A belonging to the 67th launch in the year
  1998.

* Ballistic Coefficient: Also called the first derivative of mean motion, this
  is the daily rate of change in the number of revs the object completes each day,
  divided by 2. Units are revs/day. This is "catch all" drag term used in the
  Simplified General Perturbations (SGP4) USSPACECOM predictor.

* Second Derivative of Mean Motion: The second derivative of mean motion is a
  second order drag term in the SGP4 propagator used to model terminal orbit decay.
  It measures the second time derivative in daily mean motion, divided by 6.
  Units are :math:`revs/day^3`.

* Drag Term (BSTAR): Also called the radiation pressure coefficient,
  the parameter is another drag term in the SGP4 predictor. Units are
  :math:`1/earth radii`.

Initialising the TLEs
----------------------
While a TLE can be initialised using the regular :class:`.TLE` constructor by specifying the
large number of initial parameters, by far the most usual way is to use the
:meth:`.TLE.from_tle` method, with the regular two line input from an external source
(see :ref:`tle-intro` Section for an example). Some external sources to retrieve TLEs
are listed in :ref:`tle-repositories` Section.

Another way to initialise the TLE is by orbit type. For example, initialising a
geostationary satellite TLE is as easy as defining a reference time and a target longitude:

    >>> from astropy import units as u
    >>> from satmad.propagation.tle import TLE
    >>> tle_geo = TLE.init_geo(Time("2020-06-10T12:13:14.000"), 42 * u.deg)
    >>> print(tle_GEO)
    No Name
    1 99999U 12345A   20162.50918981  .00000000  00000-0  00000+0 0    15
    2 99999   0.0000 124.6202 0000000   0.0000   0.0000  1.00273791    04

Similarly, a circular sun-synchronous orbit at 800 km reference altitude (at Equator) and at a Local Time of the
Ascending Node (LTAN) at 10:30 can be initialised simply with:

    >>> from astropy import units as u
    >>> from satmad.propagation.tle import TLE
    >>> alt = 800 * u. km
    >>> ltan = 10.5
    >>> tle_sso = TLE.init_sso(Time("2020-06-10T12:13:14.000"), alt, ltan)
    >>> print(tle_sso)
    No Name
    1 99999U 12345A   20162.50918981  .00000000  00000-0  00000+0 0    15
    2 99999  98.6032  56.8119 0000000   0.0000   0.0000 14.27530325    07


Other parameters such as eccentricity, argument of perigee or mean anomaly can be optionally set to initialise with
the requested values.

.. _tle-orbit_properties:

Checking the Orbit Properties
-----------------------------
Once the TLE is initialised, it is then possible to query orbit parameters such as inclination, eccentricity
or node rotation rate (to check whether the satellite is sun-synchronous):

    >>> tle.sm_axis()
    <Quantity 6796.50623984 km>
    >>> tle.period()
    <Quantity 5576.21313766 s>

The parameters have their quantities attached, so it is possible to easily convert them to more convenient units:

    >>> tle.inclination.to(u.deg)
    <Quantity 51.6454 deg>

Please see the example `Analysis of a Repeating Sun-Synchronous Orbit <../examples/tutorials/sso_analysis.ipynb>`_
for more information. See [OM2]_ for more information on Sun-Synchronous Orbits.


.. _tle-repositories:

Common TLE Repositories
-----------------------

* `Space-Track <https://www.space-track.org/>`_

* `Celestrak <http://celestrak.com/NORAD/elements/>`_

* `N2YO <https://www.n2yo.com/>`_

Reference/API
-------------
.. automodule:: satmad.propagation.tle
    :members:
    :undoc-members:
