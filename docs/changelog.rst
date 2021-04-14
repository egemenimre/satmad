Changelog
=========

Versions here correspond to those in `PyPI`_ and `conda-forge`_.

.. _`PyPI`: https://pypi.org/project/satmad/
.. _`conda-forge`: https://anaconda.org/conda-forge/satmad

Development Version
-------------------

The major functionalities under development are:

- Add Classical Orbital Elements
- Add resample functionality to :class:`.TimeIntervalList`
- Pass and Access finding with a simple elevation mask (e.g. Groundstation communication times)



.. _changelog-latest:

Latest Version
-----------------

Version 0.1 (05 Apr 2021)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added :class:`.TleStorage` and :class:`.TleTimeSeries` classes to load and filter multiple TLEs
  (see :ref:`tle_storage-intro`)
- Improved units handling of TLEs
- Moved analyses to a dedicated project called SatMAD Applications,
  `available at Github <https://github.com/egemenimre/satmad_applications>`_ (for Jupyter notebooks)
  and `in plain document format <https://satmad-applications.readthedocs.io/>`_.




Previous Versions
-----------------

Version 0.0.7 (24 Mar 2021)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added :class:`.GroundLocation` to model the Ground Locations (e.g. groundstations) on any planet.
- Initialised non-Earth-bound propagation
- Modified the dependencies to point to the new `pyERFA` library

Version 0.0.6 (10 Nov 2020)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Deleted the TEME Coordinate System in favour of the native Astropy version.
- Added user defined inertial coordinate systems (:class:`.CelestialBodyCRS`)

Version 0.0.5 (27 Jul 2020)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Introduced basic numerical propagation (see :ref:`numprop-intro` Section)
- Added operators (union, intersect etc.) functionality to :class:`.TimeIntervalList`

Version 0.0.4 (11 Jul 2020)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- First experimental release.

Prior to 0.0.4
^^^^^^^^^^^^^^^^^^^^^^^^^^^
All versions before 0.04 have been wildily experimental.

