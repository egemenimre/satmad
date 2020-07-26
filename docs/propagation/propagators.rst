Propagators
==================

Introduction
------------

An *orbit propagator* (or *propagator* in short) takes in an initial condition
(e.g. satellite position and velocity at a given time) and outputs where the satellite is
going to be at a target time, usually with intermediate points depicting the trajectory
while proceeding towards this target time.

The basic properties of the orbit propagators are:

*   **Initial condition type:** This could be cartesian position and velocity vector
    or orbital elements. Both of them are defined at a reference time.
*   **Integration type:** This could be analytical or numerical. Analytical propagators can compute
    the coordinates at a target time directly, without computing intermediate points. Numerical
    integrators must compute intermediate points numerically until target time is reached.
*   **Force model:** Each propagator has a force model that defines which forces (or accelerations)
    are taken into account during the propagation.

The accuracy of an orbit propagator is a function of both the force model and the analytical
or numerical model. As a general rule, a numerical propagator is more accurate than
an analytical one with a similar force model. And, among the numerical propagators, a higher order
propagation scheme generally yields better accuracy than a lower order scheme. However, there
are a lot of other parameters that determine the short and long term accuracy of a propagation
scheme.

Propagators in SatMAD
---------------------

A propagator in SatMAD accepts an initial *orbital state* (cartesian coordinates or orbital elements
at a reference time) and outputs a :class:`.Trajectory` object through a `gen_trajectory()` method.
This trajectory starts at the initial time and ends around the target time, with intermediate points
at each stepsize (given by the common :meth:`.AbstractPropagator.stepsize` property).
Note that, the trajectory may exceed the target time to ensure enough output points exist for a proper
interpolation.

Currently the following propagators are implemented in SatMAD:

.. toctree::
    :maxdepth: 1

    sgp4_propagator
    numerical_propagation

Reference/API
-------------
.. automodule:: satmad.propagation.base_propagator
    :members:
    :undoc-members: