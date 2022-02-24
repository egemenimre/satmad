# Numerical Propagation (with Scipy)

## Propagating Orbits with Scipy

Numerical propagation means integrating the equations of motion using numerical Ordinary Differential Equation (ODE) Solvers. If this sounds a mouthful, it is because it really is - there is a significant literature describing how to achieve high accuracy and computational performance for different families of problems.

The good news is that, [Scipy ODE integration algorithms](https://docs.scipy.org/doc/scipy/reference/integrate.html) handle most of the hard work and all we need to do is to describe the differential equations describing the motion.

The ODE that describes the two-body motion is:

$$
\ddot{\vec{r}} = - \dfrac{\mu}{r^3} \vec{r}
$$

where $\ddot{\vec{r}}$ is the inertial acceleration, $\vec{r}$ is the position vector (with $r$ its norm). $\mu$ is equal to the constant $GM$ (Gravitational Constant times the Mass of the main attracting body). For more complicated force models, this differential equation becomes significantly more complicated, but the procedure is still the same.

While the available propagation methods are those described by the [Scipy ODE solvers](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp), we have shown with [extensive analyses](https://satmad-applications.readthedocs.io/en/latest/analyses/propagation/num_prop_performance_1.html) that DOP853 (Explicit Runge-Kutta Order 8(5,3) due to Dormand and Prince) gives the best results in terms of accuracy and computational performance. Hence, this is the default propagator.

## Usage

The Numerical Propagator is initialised with an output stepsize (e.g. 120 seconds for a Low-Earth Orbit satellite) as well as relative and absolute tolerances that determine the propagation error. How to choose these tolerances have been investigated in this [numerical propagation analysis](https://satmad-applications.readthedocs.io/en/latest/analyses/propagation/num_prop_performance_2.html).

To generate the trajectory {py:meth}`.NumericalPropagator.gen_trajectory` method is called with an initial state and a propagation interval. A trajectory is then generated through this interval with the required output stepsize.

The following example generates an initial state, initialises a DOP853 numerical propagator with a stepsize of 120 seconds and `rtol` and `atol` of 1e-13 and 1-15, respectively. Then it generates a trajectory, starting from 1 hour after the initial state and ending 1 day later.

```python
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
```

Propagating an orbit around some other planet or celestial body is carried out similarly. The satellite has to be initialised around the celestial body (e.g. inertial frame around Moon), and the propagator has to be set up accordingly (with a force model for the Moon).

The following code initialises a satellite around Moon (frame set to `MoonCRS` or the ICRF centred at Moon), sets up the numerical propagator and propagates the satellite around the Moon starting 0.5 days after the initial condition and for a duration of 3 days.

```python
from astropy import units as u

from astropy.coordinates import (
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord
)
from astropy.time import Time

from satmad.core.celestial_bodies_lib import MOON, MoonCRS
from satmad.propagation.numerical_propagators import NumericalPropagator, ODESolverType
from satmad.utils.timeinterval import TimeInterval

# Initialises a Moon low-orbit satellite.
time = Time("2020-01-01T11:00:00.000", scale="utc")

v_moon_crs = CartesianDifferential([1, -1, 0.6], unit=u.km / u.s)
r_moon_crs = CartesianRepresentation([1000, 1000, 2000], unit=u.km)
rv_init = SkyCoord(
    r_moon_crs.with_differentials(v_moon_crs),
    obstime=time,
    frame=MoonCRS,
    representation_type="cartesian",
    differential_type="cartesian",
)

# Set up propagation config
stepsize = 60 * u.s
solver_type = ODESolverType.DOP853
rtol = 1e-13
atol = 1e-15
init_time_offset = 0.5 * u.day
duration = 3.00 * u.day

# init propagator
prop = NumericalPropagator(
    stepsize,
    solver_type=solver_type,
    rtol=rtol,
    atol=atol,
    name="",
    central_body=MOON,
)

prop_start = rv_init.obstime + init_time_offset
prop_duration = duration

# run the propagation
trajectory = prop.gen_trajectory(rv_init, TimeInterval(prop_start, prop_duration))

# print the final element
print(trajectory.coord_list[-1])
```

## Reference/API

```{eval-rst}
.. automodule:: satmad.propagation.numerical_propagators
    :members:
    :undoc-members:
```
