SGP4 Propagator
==================

.. _sgp4-intro:

Basics of SGP4
---------------

SGP4 is an analytical propagator that is a part of the `Simplified perturbations models
<https://en.wikipedia.org/wiki/Simplified_perturbations_models>`_. The input is TEME mean
classical orbital elements, usually known as `Two-Line-Elements or TLEs
<https://en.wikipedia.org/wiki/Two-line_element_set>`_ (represented in SatMAD with a
:class:`.TLE` object).

The force model contains the zonal harmonics J\ :sub:`2`, J\ :sub:`3` and J\ :sub:`4` with a
simple drag model. By definition, the propagation is Earth bound, therefore it cannot be used
with other celestial bodies.

The output of the raw propagation is True-Equator, Mean-Equinox (TEME), though the outputs are
converted to :class:`.GCRS`, which is the coordinate frame of the output trajectory.

Usage
------

The SGP4 propagator is initialised with a stepsize (e.g. 120 seconds for a Low-Earth Orbit satellite).
To generate the trajectory :meth:`.SGP4Propagator.gen_trajectory` method is called with
a correctly initialised TLE and an interval. A trajectory is then generated through this interval with
the required stepsize.

The following example generates a TLE, initialises a SGP4 with a stepsize of 120 seconds and
generates a trajectory, starting from the TLE epoch and with a duration of 1.5 days.

.. code-block:: python
    :linenos:

    # init TLE
    name = "VANGUARD 1"
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"
    tle = TLE.from_tle(line1, line2, name)

    # init SGP4
    sgp4 = SGP4Propagator(stepsize=120 * u.s)

    #  Generate trajectory
    interval = TimeInterval(tle.epoch, 1.5 * u.day)
    traj = sgp4.gen_trajectory(tle, interval)

    print(traj)

Reference/API
-------------
.. automodule:: satmad.propagation.sgp4_propagator
    :members:
    :undoc-members: