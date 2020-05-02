SatMAD: Satellite Mission Analysis and Design
---------------------------------------------------

SatMAD is an open source Python package, aiming at providing the base functionality to solve
satellite mission analysis and design as well as orbital mechanics problems with enough precision and performance
to be used in the design and operation of real satellites. The target audience is academics and amateur satellite
community, including Cubesats.

Check the source code here: https://github.com/egemenimre/satmad

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


.. toctree::
  :maxdepth: 2

  frames/index.rst

