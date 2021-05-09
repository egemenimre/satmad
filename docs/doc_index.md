# User Documentation

A word about documentation approach: There is no dedicated API documentation as such. SatMAD employs the [Astropy](https://docs.astropy.org/en/latest/coordinates/index.html) approach, where each package is introduced with its own API documentation at the end.

## Orbits and Trajectories

```{toctree} 
---
maxdepth: 2
---
coordinates/trajectory
propagation/classical_orb_elems
propagation/tle
propagation/tle_storage
propagation/propagators
propagation/force_models
```


## Time and Coordinates

```{toctree} 
---
maxdepth: 2
---
utils/timeinterval
coordinates/frames
```


## Utilities

```{toctree} 
---
maxdepth: 2
---
utils/discrete_time_events
utils/interpolators
core/ground_location
core/celestial_bodies
```

