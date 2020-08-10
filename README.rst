SatMAD: Satellite Mission Analysis and Design
---------------------------------------------
|CircleCI Status| |Codecov Status| |Documentation Status| |Astropy Badge|

SatMAD is an open source Python package, aiming at providing the base functionality to solve
satellite mission analysis and design as well as orbital mechanics problems with enough precision and performance
to be used in the design and operation of real satellites. The target audience is academics and amateur satellite
community, including Cubesats (and anyone else who might be interested).

The focus is on good documentation and good test coverage to produce a reliable
flight dynamics library.

Current functionality is:

-   Initialising an orbit with a TLE
-   Propagating the orbit analytically via the SGP4 propagator
-   Propagating the orbit numerically via the Scipy ODE Solvers (two-body force model)
-   Good infrastructure for time interval management and event finding
-   Harnessing the extensive `Astropy <http://www.astropy.org>`_ functionalities
    (e.g. finding azimuth and elevation over a ground location and coordinate frame
    transformations)

Documentation and Examples
--------------------------

The documentation for SatMAD is here: https://satmad.readthedocs.io/

You can find some hands-on Jupyter examples in the
`examples directory <https://github.com/egemenimre/satmad/tree/master/docs/examples>`_ (or
in the `documentation <https://satmad.readthedocs.io/en/latest/examples.html>`_ for a
text version).

Installing SatMAD
-----------------

The SatMAD package is on `PyPI`_ and you can install it simply by running::

    pip install satmad

You can also install it via `conda-forge`_::

    conda install -c conda-forge satmad

Do not install `satmad` using `sudo`.

You can find the source code on GitHub: https://github.com/egemenimre/satmad

.. _`PyPI`: https://pypi.org/project/satmad/
.. _`conda-forge`: https://github.com/conda-forge/satmad-feedstock

Requirements
------------

-   NumPy and SciPy are used for the underlying mathematical algorithms
-   Matplotlib is used for plots
-   Astropy handles all time and reference frame computations
-   `Portion <https://github.com/AlexandreDecan/portion>`_ handles the
    time interval mechanics
-   `Sgp4 <https://pypi.org/project/sgp4>`_ provides the TLE manipulation
    and SGP4 propagation engine
-   Pytest provides the testing framework


License
-------

This project is Copyright (c) Egemen Imre and licensed under
the terms of the GNU GPL v3+ license.

.. |Documentation Status| image:: https://readthedocs.org/projects/satmad/badge/?version=latest&token=645e1945f952813df0bb16427c4cf410850811214e4c7b6269e869291d7d8cc4
    :target: https://satmad.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |Astropy Badge| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

.. |CircleCI Status| image::  https://img.shields.io/circleci/build/github/egemenimre/satmad/master?logo=circleci&label=CircleCI
    :target: https://circleci.com/gh/satmad/satmad
    :alt: SatMAD CircleCI Status

.. |Codecov Status| image::  https://codecov.io/gh/egemenimre/satmad/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/egemenimre/satmad
    :alt: SatMAD Codecov Status
