User Guide
==========

Everything Starts with a :class:`.Trajectory`
------------------------------------------------------------------------------

At the heart of the SatMAD data structures is the :class:`.Trajectory` class.
It contains an Astropy :class:`astropy.coordinates.SkyCoord` object to keep the discrete points of the
trajectory. However, it can output the coordinates for any requested time within the time bounds of the original
:class:`astropy.coordinates.SkyCoord` trajectory object. This enables many other functionalities
not otherwise possible.


Other Topics
------------

.. toctree::
   :maxdepth: 2

   coordinates/trajectory
   propagation/tle
   coordinates/frames
   utils/timeinterval
   utils/index
