User Guide
==========

Everything Starts with a :py:class:`~satmad.coordinates.trajectory.Trajectory`
------------------------------------------------------------------------------

At the heart of the SatMAD data structures is the :py:class:`~satmad.coordinates.trajectory.Trajectory` class.
It contains an Astropy `SkyCoord` object to keep the discrete points of the trajectory. However, it can output
the coordinates for any requested time within the time bounds of the original `SkyCoord` trajectory object.
This enables many other functionalities not otherwise possible.

Other Topics
------------

.. toctree::
   :maxdepth: 2

   coordinates/trajectory
   coordinates/frames
   utils/index
