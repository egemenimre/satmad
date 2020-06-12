User Guide
==========

Orbit Initialisation
--------------------

The easiest way to start is to initialise an orbit from an existing repository, for example see the TLE of
`ISS from Celestrak <https://celestrak.com/satcat/tle.php?CATNR=25544>`_:

    >>> line1 = "1 25544U 98067A   20164.72025505  .00000382  00000-0  14906-4 0  9996"
    >>> line2 = "2 25544  51.6454  10.5753 0002538  44.8752  74.7852 15.49438622231324"
    >>> tle = TLE.from_tle(line1, line2, "ISS (ZARYA)")

It is then possible to query orbital parameters of the satellite:

    >>> tle.sm_axis()
    <Quantity 6796.50623984 km>
    >>> tle.period()
    <Quantity 5576.21313766 s>

The parameters have their quantities attached, so it is possible to easily convert them to more convenient units:

    >>> from astropy import units as u
    >>> tle.inclination.to(u.deg)
    <Quantity 51.6454 deg>

