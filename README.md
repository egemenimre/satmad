# SatMAD: Satellite Mission Analysis and Design

[![CircleCI Status](https://img.shields.io/circleci/build/github/egemenimre/satmad/master?logo=circleci&label=CircleCI)](https://circleci.com/gh/egemenimre/satmad)
[![Codecov Status](https://codecov.io/gh/egemenimre/satmad/branch/master/graph/badge.svg)](https://codecov.io/gh/egemenimre/satmad)
[![Documentation Status](https://readthedocs.org/projects/satmad/badge/?version=latest)](https://satmad.readthedocs.io/en/latest/?badge=latest)
[![Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat-square)](http://www.astropy.org/)

SatMAD is an open source Python package, aiming at providing the base functionality to solve satellite mission analysis and design as well as orbital mechanics problems with enough precision and performance to be used in the design and operation of real satellites. The target audience is academics and amateur satellite community, including Cubesats (and anyone else who might be interested).

The focus is on good documentation and good test coverage to produce a reliable flight dynamics library.

Current functionality is:

-   Loading a list of TLEs from file
-   Initialising an orbit with a TLE or Keplerian Orbital Elements
-   Propagating the orbit analytically via the SGP4 propagator
-   Propagating an orbit numerically via the Scipy ODE Solvers (two-body force model) around the Earth or another planet
-   Good infrastructure for time interval management and event finding
-   Harnessing the extensive [Astropy](http://www.astropy.org) functionalities
    (e.g. finding azimuth and elevation over a ground location and coordinate frame
    transformations)

## Documentation and Examples

The documentation for SatMAD is here: <https://satmad.readthedocs.io/>

You can find some hands-on Jupyter examples in the [examples directory](https://github.com/egemenimre/satmad/tree/master/docs/examples) (or in the [documentation](https://satmad.readthedocs.io/en/latest/examples.html) for a text version).

In addition, there is a repository of tutorials, how-to guides and analyses in a dedicated project called SatMAD Applications,
[available at Github](https://github.com/egemenimre/satmad_applications) (for Jupyter notebooks) and [in plain document format](https://satmad-applications.readthedocs.io/).


## Installing SatMAD

The SatMAD package is on [PyPI](https://pypi.org/project/satmad/) and you can install it simply by running:

    pip install satmad

You can also install it via [conda-forge](https://github.com/conda-forge/satmad-feedstock):

    conda install -c conda-forge satmad

Do not install `satmad` using `sudo`.

You can find the source code on GitHub: <https://github.com/egemenimre/satmad>


## Requirements

-   NumPy and SciPy are used for the underlying mathematical algorithms
-   Matplotlib is used for plots
-   Astropy handles all time and reference frame computations
-   [PyERFA](https://github.com/liberfa/pyerfa) handles lower level coordinate and time transformations
-   [Portion](https://github.com/AlexandreDecan/portion) handles the time interval mechanics
-   [Sgp4](https://pypi.org/project/sgp4) provides the TLE manipulation and SGP4 propagation engine
-   Pytest provides the testing framework


## License

This project is Copyright (c) Egemen Imre and licensed under the terms of the GNU GPL v3+ license.



