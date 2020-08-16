Numerical Propagation (with Scipy)
==================================

.. _numprop-intro:

Propagating Orbits with Scipy
---------------------------------------

Numerical propagation means integrating the equations of motion using numerical Ordinary
Differential Equation (ODE) Solvers. If this sounds a mouthful, it is because it really is - there
is a significant literature describing how to achieve high accuracy and computational performance
for different families of problems.

The good news is that,
`Scipy ODE integration algorithms <https://docs.scipy.org/doc/scipy/reference/integrate.html>`_
handle most of the hard work and all we need to do is to describe the differential equations
describing the motion.

The ODE that describes the two-body motion is:

    .. math::

        \ddot{\vec{r}} = - \dfrac{\mu}{r^3} \vec{r}

where :math:`\ddot{\vec{r}}` is the inertial acceleration, :math:`\vec{r}` is the position vector
(with :math:`r` its norm). :math:`\mu` is equal to the constant :math:`GM` (Gravitational Constant times
the Mass of the main attracting body). For more complicated force models, this differential equation becomes
significantly more complicated, but the procedure is still the same.

While the available propagation methods are those described by the
`Scipy ODE solvers <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp>`_,
we have shown with `extensive analyses <examples/tutorials/num_prop_performance_1.ipynb>`_ that DOP853
(Explicit Runge-Kutta Order 8(5,3) due to Dormand and Prince) gives the best results in terms of accuracy
and computational performance. Hence, this is the default propagator.

Usage
------

The Numerical Propagator is initialised with an output stepsize (e.g. 120 seconds for a Low-Earth Orbit satellite)
as well as relative and absolute tolerances that determine the propagation error.
How to choose these tolerances have been investigated here in
`this notebook <examples/analysis/num_prop_performance_2.ipynb>`_.

To generate the trajectory :meth:`.NumericalPropagator.gen_trajectory` method is called with
an initial state and a propagation interval. A trajectory is then generated through this interval with
the required output stepsize.

The following example generates an initial state, initialises a DOP853 numerical propagator with a stepsize
of 120 seconds and `rtol` and `atol` of 1e-13 and 1-15, respectively. Then it generates a trajectory, starting
from 1 hour after the initial state and ending 1 day later.

.. code-block:: python
    :linenos:

    from astropy import units as u
    from astropy.coordinates import (
        GCRS,
        CartesianDifferential,
        CartesianRepresentation,
        SkyCoord,
    )
    from astropy.time import Time
    from satmad.propagation.numerical_propagators import ODESolverType, NumericalPropagator
    from satmad.utils.timeinterval import TimeInterval

    time = Time("2004-04-06T07:51:28.386009", scale="utc")

    v_gcrs = CartesianDifferential(
        [-4.7432201610, 0.7905364950, 5.5337557240], unit=u.km / u.s
    )
    r_gcrs = CartesianRepresentation(
        [5102.50895290, 6123.01139910, 6378.13693380], unit=u.km
    )
    rv_init = SkyCoord(
        r_gcrs.with_differentials(v_gcrs),
        obstime=time,
        frame=GCRS,
        representation_type="cartesian",
        differential_type="cartesian",
    )

    # Set up propagation config
    stepsize = 120 * u.s
    solver_type = ODESolverType.DOP853
    rtol = 1e-13
    atol = 1e-15

    prop_start = rv_init.obstime + 1 * u.hr
    prop_duration = 1.0 * u.day

    # init propagator
    prop = NumericalPropagator(stepsize, rtol=rtol, atol=atol)

    # run propagation and get trajectory
    trajectory = prop.gen_trajectory(rv_init, TimeInterval(prop_start, prop_duration))

    print(trajectory)


Reference/API
-------------
.. automodule:: satmad.propagation.numerical_propagators
    :members:
    :undoc-members: