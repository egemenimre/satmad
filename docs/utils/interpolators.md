# Interpolators


## Interpolation for 3D Vectors


Interpolators package currently has only the {py:class}`.CartInterpolator3D` class. This sets up an interpolator for a set of time-varying 3D vectors (for example a trajectory). This enables retrieval of another 3D vector for a given time, as long as it is within the initial set of values used to set up the interpolator.

The interpolation is carried out individually for each axis. The input values should be smooth and continuous. The actual interpolation algorithm is {py:class}`scipy.interpolate.InterpolatedUnivariateSpline`, with the `spline_degree` defined upon
initialisation. The number of data points must be larger than the `spline_degree`.

The {py:class}`.CartInterpolator3D` class is not intended to be called or instantiated by the user. It is used for position and velocity interpolation within the {py:class}`.Trajectory` class.


## Reference/API

```{eval-rst}
.. automodule:: satmad.utils.interpolators
    :members:
    :undoc-members:
```