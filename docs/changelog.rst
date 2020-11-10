Changelog
=========

Versions here correspond to those in `PyPI`_ and `conda-forge`_.

.. _`PyPI`: https://pypi.org/project/satmad/
.. _`conda-forge`: https://anaconda.org/conda-forge/satmad

Development Version
-------------------

The major functionalities under development are:

- Add resample functionality to :class:`.TimeIntervalList`
- Pass and Access finding with a simple elevation mask (e.g. Groundstation communication times)


.. _changelog-latest:

Latest Version
-----------------

Version 0.0.6 (10 Nov 2020)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Deleted the TEME Coordinate System in favour of the native Astropy version.
- Added user defined inertial coordinate systems (:class:`.CelestialBodyCRS`)

Version 0.0.5 (27 Jul 2020)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Introduced basic numerical propagation (see :ref:`numprop-intro` Section)
- Added operators (union, intersect etc.) functionality to :class:`.TimeIntervalList`


Previous Versions
-----------------
Version 0.0.4 (11 Jul 2020)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- First experimental release.

Prior to 0.0.4
^^^^^^^^^^^^^^^^^^^^^^^^^^^
All versions before 0.04 have been wildily experimental.

