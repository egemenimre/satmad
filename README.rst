SatMAD: Satellite Mission Analysis and Design
---------------------------------------------
.. image:: https://readthedocs.org/projects/satmad/badge/?version=latest
    :target: https://satmad.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

SatMAD is an open source Python package, aiming at providing the base functionality to solve
satellite mission analysis and design as well as orbital mechanics problems with enough precision and performance
to be used in the design and operation of real satellites. The target audience is academics and amateur satellite
community, including Cubesats.

While current functionality is limited, the aim looks like this:

#. Operations support
    a) Groundstation communication times
    b) Target imaging times
#. Satellite propagation
    a) Numerical Propagation with full force model: Geopotentials, Solar Radiation Pressure, 3rd body interactions, atmospheric drag
    b) Analytical Propagation: SGP4 and Keplerian
#. Satellite orbit design and analysis
    GEO, Sun-synch and repeating orbits as well as deviations from the ideal
#. Orbit change
    Manoeuvres and Delta-V calculations
#. Attitude kinematics modelling and basic attitude laws
    Sun pointing, yaw compensation, spin
#. Satellite design
    a) Power generation with solar arrays
    b) Power consumption and battery sizing

TODO insert simple run example here

.. code-block:: python

    from satmad.examples import XXX

    XXX.runX()

Documentation
-------------

The documentation lives here:
https://satmad.readthedocs.io/


Examples
--------

You can find Jupyter notebooks with sample applications of Satmad under the `examples` directory.

Requirements
------------

- NumPy and SciPy are used for the underlying mathematical algorithms
- Astropy handles all time and reference frame computations


License
-------

This project is Copyright (c) Egemen Imre and licensed under
the terms of the GNU GPL v3+ license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause license. See the licenses folder for
more information.