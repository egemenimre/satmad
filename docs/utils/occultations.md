# Shadow and Occultations

## Computation of Penumbra, Umbra and Light

Shadow and occultation computations are very important for satellites for two main usecases:

1. Computation of eclipse intervals
2. Computation of the amount of light impinging on the satellite

The first one is useful for operations and satellite design (particularly in power and thermal subsystem sizing). The second one, to compute the amount of "light" on the satellite at any given instant, is useful during force computations in satellite propagation as well as instantaneous power generation computations.

### Computation of Eclipse Intervals

For the usecase of Umbra and Penumbra computation, we first need a satellite trajectory - the output stepsize will be the stepsize used to compute the umbra and penumbra times. Finally, the occulting body (in this example, the Earth) is specified. The {py:meth}`.occultation_intervals` method computes the umbra and penumbra intervals as {py:class}`.TimeIntervalList` classes. Note the optional specification of the illuminating body, and the ephemeris type.

    >>> # Init trajectory
    >>> pvt0 = _init_orbit_leo()
    >>>
    >>> # Set up propagation config
    >>> stepsize = 120 * u.s
    >>> prop_interval = TimeInterval(pvt0.obstime, 2.0 * u.day)
    >>>
    >>> # run propagation and get trajectory
    >>> trajectory = _init_trajectory(pvt0, stepsize, prop_interval)
    >>>
    >>> # Compute occultation intervals
    >>> occulting_body = EARTH
    >>> umbra_intervals, penumbra_intervals = occultation_intervals(trajectory, occulting_body, illum_body=SUN, ephemeris="jpl")

While the {py:meth}`.occultation_intervals` method is simple and useful, for some cases shadow intervals from multiple occulting bodies are required (e.g. the times when the GEO satellites are eclipsed by the Earth as well as the Moon). For this, the {py:meth}`.multi_body_occultation_intervals` method is used. The setup is very similar to {py:meth}`.occultation_intervals`, but the occulting bodies are defined as a list or a tuple and the results are retrieved in a dict, with the name of the occulting body being the key.

    >>> # run propagation and get trajectory
    >>> trajectory = _init_trajectory(pvt0, stepsize, prop_interval)
    >>>
    >>> # Compute occultation intervals
    >>> occult_bodies = (EARTH, MOON)
    >>> output_dict = multi_body_occultation_intervals(trajectory, occult_bodies, illum_body=SUN, ephemeris="jpl")
    >>>
    >>> # Retrieve the results
    >>> umbra_intervals_earth, penumbra_intervals_earth = output_dict[EARTH.name]
    >>> umbra_intervals_moon, penumbra_intervals_moon = output_dict[MOON.name]

As the results are given as {py:class}`.TimeIntervalList` classes, they can be intersected with other intervals (such as ground communication times) or "merged" with other intervals (such as combining Earth and Moon eclipse intervals to compute the total duration).

### Computation of the Instantaneous Occultation Geometry

For the second usecase, the instantaneous occultation geometry is computed through the {py:meth}`.compute_occultation` method. The method requires the position vectors of the satellite, occulting body (e.g. Earth), and the illuminating body (usually the Sun) as `SkyCoord` objects, as well as the {py:class}`.CelestialBody` classes of the occulting and illuminating objects. These are required to compute the sizes of the disks and flattening of the occulting body.

    >>> (illum_status, illum_ratio, penumbra_param, umbra_param, proj_distance) = compute_occultation(pos_obj, pos_occult_body, pos_illum_body, occulting_body=EARTH, illum_body=SUN)

Note that the input position vectors can be in any frame - they will be converted to the local inertial frame of the occulting body (e.g. GCRS for the Earth) to carry out the calculations. However, it seems to be much faster to convert to the correct coordinates outside this method if multiple results are requested. Furthermore, the times of the vectors are not checked to increase speed. Therefore, it is the user's responsibility to make sure that they match.

The result is a tuple containing the shadow geometry parameters as well as illumination status (light, penumbra or umbra) and illumination percentage (1 for illuminated, 0 for umbra, and between 0 and 1 for penumbra). The `penumbra_param`, `umbra_param` and `proj_distance` parameters contain critical parameters of the shadow geometry:

1. Shadow terminators may only be encountered when `proj_distance >= 0`
2. However,the object will still be illuminated if `penumbra_param > 0`
3. Object in penumbra if  `penumbra_param < 0 < umbra_param`
4. Object in umbra if `0 < umbra_param`
5. `penumbra_param = 0` is penumbra terminator
6. `umbra_param = 0` is umbra terminator

The illumination intervals are computed through these shadow geometry parameters, but for many applications they are not needed.

## Under the Hood: Occultation Model

The occultation geometry model is based on the geometry analysis presented in the paper [NASA Technical Paper 3547, "Method for the Calculation of Spacecraft Umbra and Penumbra Shadow Terminator Points", Carlos R. Ortiz Longo, Steven L. Rickman. Apr 1995](https://ntrs.nasa.gov/api/citations/19950023025/downloads/19950023025.pdf?attachment=true). A similar model is given in [Montenbruck's Satellite Orbits book [OM4]](../references.md#orbital-mechanics). The model is based on computing the "instantaneous shadow plane at satellite location" that is normal to the relative position vector between the occulting body and the illuminating body. The radii of umbra and penumbra circles are compared against the position vector of the satellite. The compared values are called `umbra_param` and `penumbra_param`, respectively.

In addition to the paper, this function takes into account the oblateness of the occulting body, assuming that the oblateness is along the local Z axis (i.e., equatorial radius is the maximum radius and polar radius is the minimum radius).

Illumination intervals are computed by collecting all umbra and penumbra parameters through the computation range and estimating the zero-crossing times through the {py:class}`.DiscreteTimeEvents` class. The accuracy of the zero-crossings are dependent on the stepsize: The more the intermediate points, the higher the accuracy but the longer the computation time. For a LEO satellite on a polar orbit, the difference between a stepsize of 180 seconds and 60 seconds is several milliseconds. Therefore, there is little point in using very small stepsizes as the gains in accuracy is likely rather limited.

Illumination percentage is also based on [Montenbruck's book (Section 3.4.2) [OM4]](../references.md#orbital-mechanics). The percentage of the area of apparent disc of the illuminating body that is still visible behind the occulting body is calculated and used as the illumination ratio.


## Reference/API

```{eval-rst}
.. automodule:: satmad.core.occultation
    :members:
    :undoc-members:
```
