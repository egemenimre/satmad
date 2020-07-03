# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2020 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Trajectory of an object with interpolators.

"""
import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, SkyCoord
from astropy.timeseries import TimeSeries

from satmad.utils.interpolators import CartInterpolator3D
from satmad.utils.timeinterval import TimeInterval


class Trajectory:
    """
    Class that keeps the coordinates of a trajectory and interpolates
    coordinates in between.

    This class keeps the coordinate points in time and lazy inits an
    interpolated a set of new
    in between coordinates when requested. This enables sub-sampling or
    supersampling the coordinates or finding special points (e.g. Equator
    crossings) with ease.

    The class uses an `InterpolatedUnivariateSpline` inside, so the spline
    function is guaranteed to pass through all provided points. The default
    spline degree is set to 5.

    """

    # Do not init interpolators by default - lazy init when requested
    _interpolators_initialised = False
    """Flag to indicate whether interpolators are already initialised."""

    _has_velocity = False
    """Flag to indicate whether there is velocity information in the supplied
    trajectory coordinates. Controls whether velocity interpolation is to
    take place or not."""

    # Default settings for the interpolator
    _spline_degree = 5
    _extrapolate_action = "raise"

    _EPS_TIME = 10 * u.us
    """Allowable time threshold, this much 'out of bounds' is allowed when
    handling requested
    interpolation times. This helps with floating point artifacts such as
    round-off errors."""

    # Shorthand to save rendering of the derived unit
    _u_km_per_s = u.km / u.s

    def __init__(self, coords_list: SkyCoord, replicate=False):
        """
        Parameters
        ----------
        coords_list : SkyCoord
            `SkyCoord` object containing the trajectory
        replicate : bool, optional
            If `True` replicates the `coords_list` into the object, otherwise
            creates a shallow copy
            of the `coords_list`. Default is False.
        """

        if replicate:
            # replicate internal coordinate data
            self._coord_list = coords_list.replicate()
        else:
            # shallow copy the internal coordinate data
            self._coord_list = coords_list

        # save the frame name for easy access when creating new instances
        self._frame_name = self._coord_list.frame[0].name

        # save the begin and end times for the internal trajectory
        # and the interpolator
        self._interval = TimeInterval(
            self._coord_list[0].obstime, self._coord_list[-1].obstime
        )

        # Check whether we have velocity data
        if self._coord_list.data.differentials:
            self._has_velocity = True

    def __call__(self, t):
        """
        Computes the interpolated coordinates at the given time(s).

        Parameters
        ----------
        t : Time
            Time or list of times where an interpolated coordinate is required

        Returns
        -------
        coords : SkyCoord
           A `SkyCoord` object with the requested coordinate(s), in the
           original frame

        Raises
        ------
        ValueError
            If the requested `t` value is out of bounds for the interpolator

        """
        if not self._interpolators_initialised:
            # Lazy init interpolators only when required
            self._init_interpolators()

        # **** get coords at requested time ****

        # compute raw r and v vectors and
        # fill Astropy cartesian vectors (with velocities if available)
        if t.isscalar:
            # there is a single time instance
            coords = CartesianRepresentation(self._compute_pos(t), copy=False)
            if self._has_velocity:
                v = CartesianDifferential(self._compute_vel(t), copy=False)
                coords = coords.with_differentials(v)
        else:
            # there is more than one time instance
            r = [self._compute_pos(time) for time in t]
            coords = CartesianRepresentation(
                np.asarray(r), unit=r[0].unit, copy=False, xyz_axis=1
            )

            if self._has_velocity:
                v = [self._compute_vel(time) for time in t]
                v = CartesianDifferential(
                    np.asarray(v), unit=v[0].unit, copy=False, xyz_axis=1
                )
                coords = coords.with_differentials(v)

        return SkyCoord(
            coords,
            obstime=t,
            frame=self._frame_name,
            representation_type="cartesian",
            differential_type="cartesian",
        )

    def _compute_pos(self, t):
        """
        Computes the position vector at the interpolated time, also checks for
        the small float artifacts at beginning or end.

        Parameters
        ----------
        t : Time
            Time or list of times where an interpolated coordinate is required

        Returns
        -------
        Position vector at the requested time.

        Raises
        ------
        ValueError
            If the requested `t` value is out of bounds for the interpolator
        """
        # check whether the requested times are _almost_ equal to begin and
        # end times.
        # this saves against floating point rounding errors
        # is the requested time very close to first time
        if abs((t - self._interval.start).to(u.us)) < self._EPS_TIME:
            # set position to first position component
            r = self._coord_list[0].cartesian.xyz

        # is the requested time very close to last time
        elif abs((t - self._interval.end).to(u.us)) < self._EPS_TIME:
            # set position to last position component
            r = self._coord_list[-1].cartesian.xyz

        # target time is either completely out of bounds or comfortably
        # inside the interpolation range
        else:
            # Convert time to "days since epoch"
            t_req = (t - self._interval.start).jd
            # interpolate to get position component at target time
            r = self._r_interpol(t_req) * u.km

        return r

    def _compute_vel(self, t):
        """
        Computes the velocity vector at the interpolated time, also checks for
        the small float artifacts at beginning or end.

        Parameters
        ----------
        t : Time
            Time or list of times where an interpolated coordinate is required

        Returns
        -------
        Velocity vector at the requested time.

        Raises
        ------
        ValueError
            If the requested `t` value is out of bounds for the interpolator
        """
        # check whether the requested times are _almost_ equal to begin and
        # end times
        # this saves against floating point rounding errors
        # is the requested time very close to first time
        if abs((t - self._interval.start).to(u.us)) < self._EPS_TIME:
            # set velocity to first velocity component
            v = self._coord_list[0].velocity.d_xyz

        # is the requested time very close to last time
        elif abs((t - self._interval.end).to(u.us)) < self._EPS_TIME:
            # set velocity to last velocity component
            v = self._coord_list[-1].velocity.d_xyz

        # target time is either completely out of bounds or comfortably
        # inside the interpolation range
        else:
            # Convert time to "days since epoch"
            t_req = (t - self._interval.start).jd

            # interpolate to get velocity  component at target time
            v = self._v_interpol(t_req) * self._u_km_per_s

        return v

    def to_timeseries(self):
        """
        Converts the internal `coord_list` of the trajectory to a table for easy
        manipulation and file output.

        Returns
        -------
        TimeSeries
            Coordinate list as `TimeSeries`
        """
        r = self.coord_list.cartesian

        # generate timeseries

        # add velocity info if possible
        if self._has_velocity:
            v = self.coord_list.velocity
            ts = TimeSeries(
                time=self.coord_list.obstime.isot,
                data={
                    "r_X": r.x,
                    "r_Y": r.y,
                    "r_Z": r.z,
                    "v_X": v.d_x,
                    "v_Y": v.d_y,
                    "v_Z": v.d_z,
                },
            )
        else:
            ts = TimeSeries(
                time=self.coord_list.obstime.isot,
                data={"r_X": r.x, "r_Y": r.y, "r_Z": r.z},
            )

        return ts

    @property
    def coord_list(self):
        """Returns the list of underlying coordinates forming the
        trajectory."""
        return self._coord_list

    @property
    def interval(self):
        """Returns the interval of validity for the
        trajectory."""
        return self._interval

    @property
    def interpolator_name(self) -> str:
        """Returns the name of the interpolator."""
        return self._r_interpol.interpolator_name

    @classmethod
    def reqd_min_elements(cls):
        """Returns the required minimum number of elements in the trajectory."""
        return cls._spline_degree

    def _init_interpolators(self):
        """
        Initialises position and velocity interpolators.
        """
        # Init time list in days
        t_list = (self._coord_list.obstime - self._interval.start).jd
        # Init pos interpolators
        self._init_pos_interpolators(t_list)

        # init velocity interpolators if velocity info present
        if self._has_velocity:
            self._init_vel_interpolators(t_list)

        # Finally, set interpolators flag
        self._interpolators_initialised = True

    def _init_pos_interpolators(self, t_list):
        """
        Initialises position interpolators (all values in km).

        Parameters
        ----------
        t_list : (N,) array_like
            Input dimension of data points -- must be strictly increasing

        """
        r_x_list = self._coord_list.cartesian.xyz[0, :].to(u.km).value
        r_y_list = self._coord_list.cartesian.xyz[1, :].to(u.km).value
        r_z_list = self._coord_list.cartesian.xyz[2, :].to(u.km).value

        self._r_interpol = CartInterpolator3D(
            t_list,
            r_x_list,
            r_y_list,
            r_z_list,
            spline_degree=self._spline_degree,
            extrapolate_action=self._extrapolate_action,
        )

    def _init_vel_interpolators(self, t_list):
        """
        Initialises velocity interpolators (all values in km/s).

        Parameters
        ----------
        t_list : (N,) array_like
            Input dimension of data points -- must be strictly increasing

        """
        v_x_list = self._coord_list.velocity.d_xyz[0].to(self._u_km_per_s).value
        v_y_list = self._coord_list.velocity.d_xyz[1].to(self._u_km_per_s).value
        v_z_list = self._coord_list.velocity.d_xyz[2].to(self._u_km_per_s).value

        self._v_interpol = CartInterpolator3D(
            t_list,
            v_x_list,
            v_y_list,
            v_z_list,
            spline_degree=self._spline_degree,
            extrapolate_action=self._extrapolate_action,
        )

    def __str__(self):
        """String representation of the object."""
        return (
            f"Trajectory from {self._interval.start} to {self._interval.end} i"
            f"n frame {self._frame_name}. "
            f"(Interpolators initialised: "
            f"{self._interpolators_initialised})"
        )
