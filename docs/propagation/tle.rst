Orbit Initialisation with Two-Line Elements (TLEs)
====================================================

Introduction
------------
A two-line element set (TLE) is a data format containing a set of TEME
(True Equator, Mean Equinox) mean orbital elements of an Earth-orbiting object
for a given point in time, called the Epoch Time.

These orbital elements are solely for use with the SGP4 propagator as the two are coupled
with the underlying analytical orbit theory.

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
:meth:`.TLE.from_tle` method, with the regular two line input from an external source.

Reference/API
-------------
.. automodule:: satmad.propagation.tle
    :members:
    :undoc-members:
