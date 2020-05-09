Additional Frames (J2000 and TIRS)
==================================

Introduction
------------
The built-in frames offered by `Astropy <https://docs.astropy.org/en/latest/coordinates/index.html>`_
do not include some frames that are used in satellite applications. To bridge this gap, this package
offers Terrestrial Intermediate Reference System (:class:`.TIRS`) and
Mean Pole and Equinox at J2000.0 Reference System (:class:`.J2000`).

The definitions and conversions for :class:`.TIRS` are from IERS Conventions 2010,
Chapter 5 [Ref2]_. Put simply, when the TIRS coordinate is multiplied by the Polar Motion Matrix, the final
coordinates are in ITRS:

.. math:: \vec{r}_{ITRS} = W \times \vec{r}_{TIRS}

:class:`.J2000` coordinate frame is similar to GCRS but rotated by a constant frame bias
[Ref1]_:

.. math:: \vec{r}_{J2000} = B \times \vec{r}_{GCRS}

This rotation is applicable only to the equinox based approach, and is only an approximation.
The difference between GCRS and J2000 is less than 1m for the Low Earth Orbit, therefore these two
can be used interchangeably with a small error.





Usage
------
The :class:`.J2000` and :class:`.TIRS` classes
are similar to (and compatible with) the `Astropy Built-in Frames
<https://docs.astropy.org/en/latest/coordinates/index.html#built-in-frame-classes>`_.

.. code-block:: python
    :linenos:

    from astropy import units as u
    from astropy.coordinates import CartesianRepresentation, CartesianDifferential
    from satmad.coordinates.j2000 import J2000

    v_j2000 = CartesianDifferential([-4.7432196000, 0.7905366000, 5.5337561900], unit=u.km / u.s)
    r_j2000 = CartesianRepresentation([5102.50960000, 6123.01152000, 6378.13630000], unit=u.km)
    rv_j2000 = J2000(r_j2000.with_differentials(v_j2000), obstime=time, representation_type="cartesian", differential_type="cartesian")


Reference / API
---------------

.. automodule:: satmad.coordinates.frames
    :members:
    :undoc-members:






