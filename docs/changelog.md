# Changelog

Versions here correspond to those in [PyPI](https://pypi.org/project/satmad/) and [conda-forge](https://anaconda.org/conda-forge/satmad).


## Development Version

The major functionalities under development are:

- Introduction of Body Fixed Frame
- Development of a more generic AltAz class.

## Latest Version

### Version 0.1.2 (29 May 2021)

- Added [occultations and shadows](utils/occultations.md)
- Introduced [True-of-Date and J2000 Equatorial frames](coordinates/frames.md)
- Initialised the [DiscreteTimeData interpolator](utils/interpolators.md) 
- Migrate Documentation to [MyST Markdown](https://myst-parser.readthedocs.io)

## Previous Versions

### Version 0.1.1 (14 Apr 2021)

- Added [Keplerian Orbital Elements](propagation/classical_orb_elems.md)
- Added new tutorials on [Numerical Orbit Propagation](tutorials/numerical_prop_1.ipynb)


### Version 0.1 (05 Apr 2021)

- Added {py:class}`.TleStorage` and {py:class}`.TleTimeSeries` classes to load and filter [multiple TLEs](propagation/tle_storage.md)
- Improved units handling of TLEs
- Moved analyses to a dedicated project called SatMAD Applications, [available at Github](https://github.com/egemenimre/satmad_applications) (for Jupyter notebooks) and [in plain document format](https://satmad-applications.readthedocs.io/).


### Version 0.0.7 (24 Mar 2021)

- Added {py:class}`.GroundLocation` to model the Ground Locations (e.g. groundstations) on any planet.
- Initialised non-Earth-bound propagation
- Modified the dependencies to point to the new `pyERFA` library

###  Version 0.0.6 (10 Nov 2020)

- Deleted the TEME Coordinate System in favour of the native Astropy version.
- Added user defined inertial coordinate systems ({py:class}`.CelestialBodyCRS`)

### Version 0.0.5 (27 Jul 2020)

- Introduced basic numerical propagation (see [Numerical Propagators](propagation/numerical_propagation.md))
- Added operators (union, intersect etc.) functionality to {py:class}`.TimeIntervalList`

### Version 0.0.4 (11 Jul 2020)

- First experimental release.

### Prior to 0.0.4

All versions before 0.04 have been wildly experimental.
