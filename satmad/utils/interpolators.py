# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Interpolators used throughout the project.

"""
from scipy import interpolate as ip


class CartInterpolator3D:
    """
    One-dimensional interpolating spline for a given set of 3D data points.

    For each of the data sets, fits a spline y = spl(t) of degree
    `spline_degree` to the provided `t`, `y` data.
    Uses a `InterpolatedUnivariateSpline` inside, so the spline function
    passes through all provided points.

    Parameters
    ----------
    t : (N,) array_like
        Input dimension of data points -- must be strictly increasing
    x : (N,) array_like
        input dimension of data points
    y : (N,) array_like
        input dimension of data points
    z : (N,) array_like
        input dimension of data points
    spline_degree : int, optional
        Degree of the smoothing spline.  Must be 1 <= `k` <= 5.
    extrapolate_action : int or str, optional
        Controls the extrapolation mode for elements
        not in the interval defined by the knot sequence.

        * if ext=0 or 'extrapolate', return the extrapolated value.
        * if ext=1 or 'zeros', return 0
        * if ext=2 or 'raise', raise a ValueError
        * if ext=3 of 'const', return the boundary value.

        The default value is 0.

     See Also
     --------
     scipy.interpolate.InterpolatedUnivariateSpline : The interpolator
        used inside this class.

     Notes
     -----
     The number of data points must be larger than the `spline_degree`.
    """

    def __init__(self, t, x, y, z, spline_degree=5, extrapolate_action="raise"):
        # init interpolators
        self._r_x_interpol = ip.InterpolatedUnivariateSpline(
            t, x, k=spline_degree, ext=extrapolate_action
        )
        self._r_y_interpol = ip.InterpolatedUnivariateSpline(
            t, y, k=spline_degree, ext=extrapolate_action
        )
        self._r_z_interpol = ip.InterpolatedUnivariateSpline(
            t, z, k=spline_degree, ext=extrapolate_action
        )

        # set the interpolator class name
        self._interpolator_name = type(self._r_x_interpol).__name__

    @property
    def interpolator_name(self) -> str:
        """Returns the name of the interpolator."""
        return self._interpolator_name

    def __call__(self, t, nu=0, ext=None) -> list:
        """
         Evaluate spline (or its nu-th derivative) at positions `t`.

         Parameters
         ----------
         t : array_like
             A 1-D array of points at which to return the value of the
             smoothed spline or its derivatives. Note: t can be unordered
             but the evaluation is more efficient if t is (partially) ordered.
         nu  : int
             The order of derivative of the spline to compute.
         ext : int
             Controls the value returned for elements of ``t`` not in the
             interval defined by the knot sequence.

             * if ext=0 or 'extrapolate', return the extrapolated value.
             * if ext=1 or 'zeros', return 0
             * if ext=2 or 'raise', raise a ValueError
             * if ext=3 or 'const', return the boundary value.

             The default value is 0, passed from the initialization of
             the spline.

         Returns
         -------
         r: list
            List containing interpolated x, y and z values

        Raises
        ------
        ValueError
            If the requested `t` value is out of bounds for the interpolator
            (if `extrapolate_action` is set to `raise`)
         """

        # t limit check is carried out within the interpolators themselves

        # generate the interpolated values
        r = [
            self._r_x_interpol(t, nu, ext),
            self._r_y_interpol(t, nu, ext),
            self._r_z_interpol(t, nu, ext),
        ]

        return r

    def __str__(self):
        return (
            f"3D Cartesian Interpolator ({self.interpolator_name}) "
            f"with {len(self._r_x_interpol.get_knots())} "
            f"knots for each axis."
        )
