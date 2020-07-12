Project Details
===============

What is SatMAD?
---------------

SatMAD is an open source Python package, aiming at providing the base functionality to solve
satellite mission analysis and design as well as orbital mechanics problems with enough precision and performance
to be used in the design and operation of real satellites. The target audience is academics and amateur satellite
community, including (but not limited to) Cubesats.


Current Status
---------------

Current functionality is:

-   Initialising an orbit with a TLE (see :ref:`tle-intro` Section)
-   Propagating the orbit via SGP4 propagator (see :ref:`sgp4-intro` Section)
-   Good infrastructure for time interval management and event finding
-   Harnessing the extensive `Astropy <http://www.astropy.org>`_ functionalities
    (e.g. finding azimuth and elevation over a ground location and coordinate frame
    transformations)

What's New?
-----------

Check the :ref:`changelog-latest` for the changelog and recently added functionalities.

Future Functionality
--------------------

1. Operations support
    a) Groundstation communication times
    b) Target imaging times
2. Satellite propagation
    a) Numerical Propagation with full force model: Geopotentials, solar radiation pressure, 3rd body interactions, atmospheric drag
    b) Analytical Propagation: Two-Body
3. Satellite orbit design and analysis
    a) GEO, Sun-synch and repeating orbits
    b) Analysis of deviations from the ideal orbits
4. Orbit change
    a) Manoeuvres and Delta-V calculations
5. Attitude kinematics modelling and basic attitude laws
    a) Sun pointing, yaw compensation, spin
6. Satellite design
    a) Power generation with solar arrays
    b) Power consumption and battery sizing
    c) Propellant budget


License
-------

This project is Copyright (c) Egemen Imre and licensed under
the terms of the GNU GPL v3+ license.

About the Author
----------------
My name is `Egemen Imre <https://twitter.com/uyducusirin>`_, the author of SatMAD.
While this project has really been an excuse to learn Python for me,
for 20 years I have been developing Orbital Mechanics software professionally for topics ranging from
satellite mission analysis and design to actual satellite operations.

Acknowledgements
----------------
While setting up the code environment and continuous integration, I have taken a lot of inspiration from
`Astropy <https://www.astropy.org/>`_ and `Poliastro <https://github.com/poliastro/poliastro>`_.
Astropy is an amazing Astronomy Library for Python (which is doing a lot of the tedious infrastructure job for SatMAD)
and Poliastro is probably the most actively developed Astrodynamics and Orbital Mechanics library in Python,
with a focus on interplanetary flight.

