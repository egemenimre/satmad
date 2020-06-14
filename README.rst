SatMAD: Satellite Mission Analysis and Design
---------------------------------------------
|CircleCI Status| |Documentation Status| |Astropy Badge|

SatMAD is an open source Python package, aiming at providing the base functionality to solve
satellite mission analysis and design as well as orbital mechanics problems with enough precision and performance
to be used in the design and operation of real satellites. The target audience is academics and amateur satellite
community, including Cubesats (and anyone else who might be interested).

The focus is on good documentation and good test coverage to produce a reliable
flight dynamics library.

While current functionality is limited, the aim looks like this:

#. Operations support
    a) Groundstation communication times
    b) Target imaging times
#. Satellite propagation
    a) Numerical Propagation with full force model: Geopotentials, Solar Radiation Pressure, 3rd body interactions, atmospheric drag
    b) Analytical Propagation: SGP4 and Keplerian
#. Satellite orbit design and analysis
    a) GEO, Sun-synch and repeating orbits
    b) Analysis of deviations from the ideal orbits
#. Orbit change
    a) Manoeuvres and Delta-V calculations
#. Attitude kinematics modelling and basic attitude laws
    a) Sun pointing, yaw compensation, spin
#. Satellite design
    a) Power generation with solar arrays
    b) Power consumption and battery sizing
    c) Propellant budget


Documentation
-------------

The documentation for SatMAD is here:
https://satmad.readthedocs.io/



Requirements
------------

- NumPy and SciPy are used for the underlying mathematical algorithms
- Astropy handles all time and reference frame computations
- Portion handles the time interval mechanics
- python-sgp4 provides the SGP4 propagation engine
- Pytest provides the testing framework


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
