TLE
=====

Introduction
------------
A two-line element set (TLE) is a data format encoding a list of TEME
(True Equator, Mean Equinox) mean orbital elements
of an Earth-orbiting object for a given point in time, called the Epoch Time.

These orbital elements are solely for use with the SGP4 propagator due to the
analytical orbit theory used in its derivation.

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

In addition to the usual classical orbital elements, the components of the TLE are
(adapted from the `NASA definition <https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/SSOP_Help/tle_def.html>`_
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
  second order drag term in the SGP4 predictor used to model terminal orbit decay.
  It measures the second time derivative in daily mean motion, divided by 6.
  Units are :math: `revs/day^3`.

* Drag Term (BSTAR): Also called the radiation pressure coefficient,
  the parameter is another drag term in the SGP4 predictor. Units are
  :math: `1/earth radii`.

Reference/API
-------------
.. automodule:: satmad.propagation.tle
    :members:
    :undoc-members:
