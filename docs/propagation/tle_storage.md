# Working with Multiple TLEs


## Lists of Mixed TLEs: {py:class}`.TleStorage`

Many applications require working with TLEs from a number of satellites (for example satellites from the same launch or belonging to the same constellation). SatMAD enables loading and flexible filtering of such TLE lists using the {py:class}`.TleStorage` class.

For this usecase, the most common application is to load a file (or a string) containing TLE data in text format and filtering satellites with a certain parameter:

```python
from pathlib import Path
from astropy.time import Time
from satmad.propagation.tle_storage import TleFilterParams, TleStorage

# use case 1, load from TLE file, filter for a certain satellite number
file_path = Path("my_data_dir", "my_tle_file.txt")

tle_storage_1 = TleStorage.from_path(file_path)
filtered_list_1 = tle_storage_1.filter_by_value(TleFilterParams.SAT_NUMBER, 46495)

# use case 2, load from TLE string, filter for a certain epoch
with open(file_path, "r") as f:
    tle_source_str = f.read()

tle_storage_2 = TleStorage.from_string(tle_source_str)

threshold_time = Time("2021-04-01T00:00:00.000")
filtered_list_2 = tle_storage_2.filter_by_range(TleFilterParams.EPOCH, min_value=threshold_time)
```

The example above shows a number of important functionalities. A {py:class}`.TleStorage` object can be initialised using {py:meth}`.TleStorage.from_path` method (from a TLE file) or {py:meth}`.TleStorage.from_string` method (from a string containing a list of TLEs). The latter can be useful to download TLE data from a remote location without saving it to a file. The example shows filtering functionality, which is detailed in the next section.

## Extracting Specific Data from the Lists (Filtering)

Once initialised, the TLE list can be filtered using the enumerator {py:class}`.TleFilterParams` and a filtering value. To filter with an element like an identifier (e.g. name or catalogue number) where an exact match can be found, `filter_by_value` method should be used. In the example above, a satellite with the catalogue number "46495" is extracted from the list.

Conversely, if a range rather than an exact match is sought (e.g. semimajor axis, epoch or eccentricity), then `filter_by_range` method should be used. In the example above, a threshold epoch is given and all TLE values after this threshold are extracted. This method can accept a minimum, a maximum or both, such that:

    max_value > param > min_value

Continuing the example above, we can illustrate a few other filtering examples:

```
# filtering with a max and min epoch dates
min_threshold_time = Time("2021-04-01T00:00:00.000")
max_threshold_time = Time("2021-06-01T00:00:00.000")
filtered_list_3 = tle_storage_2.filter_by_range(TleFilterParams.EPOCH, min_value=min_threshold_time, max_value=max_threshold_time)

from astropy import units as u

# filtering with a min inclination
min_inclination = 90 * u.deg
filtered_list_4 = tle_storage_2.filter_by_range(TleFilterParams.INCLINATION, min_value=min_inclination)

# filtering with a max semimajor axis
max_sma = 7000 * u.km
filtered_list_5 = tle_storage_2.filter_by_range(TleFilterParams.SM_AXIS, max_value=max_sma)
```

One important thing to note is that, the filtering value should be given as a {py:class}`astropy.units.quantity.Quantity` i.e., with a unit as the values are kept as Quantities under the hood. The code simply does not know how to interpret "7000" without units (kilometres or metres?). This saves from the all too usual interfacing and unit specification errors.

The second thing to note is that the resulting filtered list is another `TleStorage` object, with its internal `tle_list` backed by the source object. In other words, the filtering does not create new TLE objects. Clearly, if the backing list is somehow changed, then all the other filtered lists may change as well.

Finally, if the filter results in an empty list, then a new `TleStorage` object is still returned, with an empty internal `tle_list`.

In addition to the `filter_by_value` and `filter_by_range` methods, there is a third and very powerful method to filter the TLEs through user defined functions. In `filter_by_func`, a user-defined function takes one parameter related to the TLE and returns `True` or `False` through some test. It is even possible for this filter function to accept the entire TLE (as opposed to a single parameter, for more complicated computations). The semimajor axis filter above can be written with this method, sending the semimajor axis
or the entire TLE.

```python
from pathlib import Path
from astropy import units as u
from satmad.propagation.tle_storage import TleFilterParams, TleStorage

# load TLEs from file
file_path = Path("my_data_dir", "my_tle_file.txt")
tle_storage = TleStorage.from_path(file_path)

# define the filter function and filter the list
def sma_filter(a):
    """Semimajor axis filter min/max."""
    return True if 7100 * u.km > a > 7000 * u.km else False

filtered_list_sma_1 = tle_storage.filter_by_func(TleFilterParams.SM_AXIS, sma_filter)

# define the filter function and filter the list
def tle_filter(tle):
    """Semimajor axis filter min/max."""
    return True if 7100 * u.km > tle.sm_axis > 7000 * u.km else False

filtered_list_sma_2 = tle_storage.filter_by_func(TleFilterParams.TLE, tle_filter)
```

Note that, `TleFilterParams.TLE` is available for `filter_by_func` method only, as the behaviour with an entire TLE is not well defined with `filter_by_value` and `filter_by_range` methods.

## Lists of TLEs of the Same Satellite: {py:class}`.TleTimeSeries`

A specific class :class:`.TleTimeSeries` exists to store the multiple TLEs from a single satellite with time ordering (epoch values). This is particularly useful to plot the orbit or to feed it to a propagator to propagate the satellite orbit through multiple TLEs. The ideal way to initialise it is to initialise a `TleStorage` from a file or another source and then filter for a single orbit.

The code below initialises a `TleStorage` from a file and filters for the satellites with the catalogue number `37791`.

```python
from pathlib import Path
from satmad.propagation.tle_storage import TleStorage

file_path = Path("my_data_dir", "my_tle_file.txt")
tle_timeseries = TleStorage.from_path(file_path).to_tle_timeseries(37791)
````
Filtering is possible using the same methods as `TleStorage` given in the [TLE Filtering](#extracting-specific-data-from-the-lists-filtering) section.

## Reference/API

```{eval-rst}
.. automodule:: satmad.propagation.tle_storage
    :members:
    :undoc-members:
```