# Trajectory


## Representing Trajectories in SatMAD


In SatMAD, all trajectories, from aircraft to satellites to entire planets, are represented with the {py:class}`.Trajectory` class. It not only keeps the discrete position (and velocity) points in time, it also enables the user to compute interpolated coordinates within these discrete data points.

The :class:`.Trajectory` class initialised with an Astropy {py:class}`astropy.coordinates.SkyCoord` object to keep these discrete points of the trajectory. These discrete points can be retrieved by the `coord_list` property. Once the trajectory is initialised, these points should *not* be changed, as the underlying interpolators will not be updated. Note that, the input discrete points should be representing a smooth and continuous trajectory to ensure a reliable and accurate interpolation. As such, any quasi-discontinuous motion (such as sharp turns or near-instantaneous increases in velocity) cannot be modeled perfectly and will be "smoothed out".

Any point (or set of points) in time  can be computed by a direct call to the trajectory object - as long as they are within the discrete data points:

    >>> # single point in time yielding a single interpolated position
    >>> time = Time(2458826.3, format="jd", scale="utc")
    >>> rvt = trajectory(time)
    >>>
    >>> # multiple points in time yielding a multiple element interpolated position
    >>> # in a single SkyCoord object.
    >>> times = ['2010-01-01T02:00:00.000', '2010-01-01T02:10:00']
    >>> time_list = Time(times, format="isot", scale="utc")
    >>> rvt_list = trajectory(time_list)

The result is a {py:class}`astropy.coordinates.SkyCoord` object containing all the requested coordinates. The interval of validity can be queried with the `interval` property. If the requested time is out of bounds for the interpolator, a `ValueError` will be raised.

## Interpolation Method and Accuracy


The actual interpolation is carried out by the {py:class}`.CartInterpolator3D` class for position and velocity. The interpolation algorithm is a 5<sup>th</sup> order Spline Interpolator supplied by {py:class}`scipy.interpolate.InterpolatedUnivariateSpline`. The interpolation error is usually larger within the first (and last) few data points, quickly converging to a smaller value as more data points are made available to the interpolator. As the interpolator goes through each *knot point* (in this case each discrete point initialising the trajectory), there is no error in these points.

The tests show that the errors between the first two points (or edge) can be as high as 16 mm, quickly converging to less than 0.1 mm for the remaining points in the trajectory. This ISS case represents the worst case and other tests yield maximum errors on the order of 0.04 mm for the edge and 0.0003 mm nominal for a Highly Elliptical Orbit and less than 0.0005 mm at the edge for a Geostationary Orbit.

## Output of Trajectory Data


Trajectory data can be exported to a {py:class}`astropy.timeseries.TimeSeries` object through the {py:meth}`.Trajectory.to_timeseries` method. The resulting `TimeSeries` data then can be written to a file or can be manipulated. While this sort of "time, position, velocity" information is useful for visualisation and data transfer, the metadata information such as coordinate system and information regarding how this data is generated cannot be stored in this way.

## Reference / API

```{eval-rst}
.. automodule:: satmad.coordinates.trajectory
    :members:
    :undoc-members:
```