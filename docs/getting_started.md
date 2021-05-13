# Getting Started

## Getting the Package and the Code

The SatMAD package is on [PyPI](https://pypi.org/project/satmad/) and you can install it simply by running:

    pip install satmad

You can also install it via [conda-forge](https://anaconda.org/conda-forge/satmad):

    conda install -c conda-forge satmad

Do not install `satmad` using `sudo`.

You can find the source code on GitHub: <https://github.com/egemenimre/satmad>


## Dependencies

-   NumPy and SciPy are used for the underlying mathematical algorithms
-   Matplotlib is used for plots
-   Astropy handles all time and reference frame computations
-   [PyERFA](https://github.com/liberfa/pyerfa) handles lower level coordinate and time transformations
-   [Portion](https://github.com/AlexandreDecan/portion) handles the time interval mechanics
-   [Sgp4](https://pypi.org/project/sgp4) provides the TLE manipulation and SGP4 propagation engine
-   [Jplephem](https://github.com/brandon-rhodes/python-jplephem/) for JPL Ephemeris calculations