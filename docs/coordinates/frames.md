# Coordinate Systems and Frames

In SatMAD, frames and coordinate systems as well as conversions between them are handled through [Astropy](https://docs.astropy.org/en/latest/coordinates/index.html). This section introduces the additional frames defined by SatMAD.

## Local Frames of Celestial Bodies

In addition to Earth, it is possible to define a Celestial Body in space through the {py:class}`.CelestialBody` class. For each such Celestial Body, it is then possible to realise local reference frames.

The first of them is a local inertial or "Celestial Reference System" (CRS) ({py:class}`.CelestialBodyCRS`). This class simply translates the ICRS coordinate system to the centre of the celestial body. As such, the XY plane does not correspond to the equator plane.  This enables the user to define coordinates in this local inertial coordinate system and then (as an example) run an orbit propagation around it or run an analysis. For example, the following would define a "Sun CRS" (equivalent to Heliocentric Celestial Reference System), and a "Moon CRS" by simply subclassing {py:class}`.CelestialBodyCRS`. 

```python
from satmad.coordinates.frames import CelestialBodyCRS

class SunCRS(CelestialBodyCRS):
    body = "Sun"


class MoonCRS(CelestialBodyCRS):
    body = "Moon"
    ephemeris_type = "jpl"
```

In the latter, the optional `ephemeris_type` parameter determines which ephemeris to use when computing the location of the Celestial Bodies. Note that, while some frames like {py:class}`.MoonCRS` and {py:class}`.MarsCRS` are predefined, and other Planets can be easily created as shown.

Similarly, three other frames can be defined:

- {py:class}`.CelestialBodyFixed`: The planet Body Fixed Equatorial frame, similar to GCRS for the Earth.
- {py:class}`.CelestialBodyJ2000Equatorial`: An inertial equatorial frame with the alignment fixed at J2000 Epoch. This is similar to J2000 frame for the Earth.
- {py:class}`.CelestialBodyTODEquatorial`: An equatorial frame with its orientation computed using the instantaneous orientation due to planetary nutation and precession:

$$
\vec{r}_{CRS} = R_x(90+ \alpha) R_z(90- \delta)\times \vec{r}_{TOD}
$$

Inspecting the preset Mars frames, new frames for other planets can be defined:

```python
from satmad.coordinates.frames import CelestialBodyJ2000Equatorial, CelestialBodyTODEquatorial
from satmad.core.celestial_bodies_lib import MarsCRS


class MarsTODEquatorial(CelestialBodyTODEquatorial):
    body_name = "Mars"
    cb_crs = MarsCRS

class MarsJ2000Equatorial(CelestialBodyJ2000Equatorial):
    body_name = "Mars"
    cb_crs = MarsCRS

class MarsBodyFixed(CelestialBodyTODEquatorial):

    body_name = "Mars"
    cb_crs = MarsCRS

```
Similar to the {py:class}`.CelestialBody` class, `body_name` defines the name of the body, such that its coordinates can be computed within the ICRS using builtin or JPL DE4XX ephemerides. The `cb_crs` parameter links these frames to CRS inertial frame of this body, as SatMAD (or Astropy) would not know how to chain the coordinate conversions from this frame to (for example) HCRS.

The definition of these coordinate systems require the modelling of how the planet is oriented with respect to ICRS. This is done via the Celestial Body Rotation parameters: right ascension and declination of the celestial body North Pole, and prime meridian angle. As can be imagined, these are different for each planet - the values used in this model are taken from [Report on IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015 [TCF2]](../references.md#time-and-coordinate-frames), except for the Moon, which is taken from the 2009 version of the same report. As such, these values are incompatible with, yet much more up-to-date than GMAT and [NASA HORIZONS Web Interface](https://ssd.jpl.nasa.gov/horizons.cgi), which use the much simpler 2000 version of these models. The TOD and J2000 Equatorial frame definitions are given in [[OM1], [OM3] and [OM4]](../references.md#time-and-coordinate-frames).


The {py:class}`.CelestialBodyCRS` class corresponds to setting the frame as "ICRF" around a central body in GMAT. The {py:class}`.CelestialBodyJ2000Equatorial` is equivalent to "Body Inertial" and {py:class}`.CelestialBodyTODEquatorial` to "Body Equatorial" in GMAT. 

## Earth Based Additional Frames (J2000)

The built-in frames offered by [Astropy](https://docs.astropy.org/en/latest/coordinates/index.html) do not include some frames that are used in satellite applications. To bridge this gap, this package offers Mean Pole and Equinox at J2000.0 Reference System ({py:class}`.J2000`).

{py:class}`.J2000` coordinate frame is similar to GCRS but rotated by a constant frame bias
[[TCF1]](../references.md#time-and-coordinate-frames):

$$
\vec{r}_{J2000} = B \times \vec{r}_{GCRS}
$$

This rotation is applicable only to the equinox based approach, and is only an approximation. The difference between GCRS and J2000 is less than 1 metre for the Low Earth Orbit, therefore these two can be used interchangeably with a small error.

The {py:class}`.J2000` class is similar to (and compatible with) the [Astropy Built-in Frames](https://docs.astropy.org/en/latest/coordinates/index.html#built-in-frame-classes).


```python
from astropy import units as u
from astropy.coordinates import CartesianRepresentation, CartesianDifferential
from satmad.coordinates.frames import J2000
from astropy.time import Time

time = Time("2020-01-11T11:00:00.000", scale="utc")
v_j2000 = CartesianDifferential([-4.7432196000, 0.7905366000, 5.5337561900], unit=u.km / u.s)
r_j2000 = CartesianRepresentation([5102.50960000, 6123.01152000, 6378.13630000], unit=u.km)
rv_j2000 = J2000(r_j2000.with_differentials(v_j2000), obstime=time, representation_type="cartesian", differential_type="cartesian")
```

## Reference / API

```{eval-rst}
.. automodule:: satmad.coordinates.frames
    :members:
    :undoc-members:
    
.. automodule:: satmad.coordinates.planetary_rot_elems
    :members:
    :undoc-members:
```
