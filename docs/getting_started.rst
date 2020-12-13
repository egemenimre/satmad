Getting Started
===============

Getting the Package and the Code
--------------------------------

The SatMAD package is on `PyPI`_ and you can install it simply by running::

    pip install satmad

You can also install it via `conda-forge`_::

    conda install -c conda-forge satmad

Do not install `satmad` using `sudo`.

You can find the source code on GitHub: https://github.com/egemenimre/satmad

.. _`PyPI`: https://pypi.org/project/satmad/
.. _`conda-forge`: https://anaconda.org/conda-forge/satmad

Dependencies
------------
- NumPy and SciPy are used for the underlying mathematical algorithms
- Matplotlib is used for plots
- Astropy handles all time and reference frame computations
- PyERFA handles lower level coordinate and time transformations
- Portion handles the time interval mechanics
- python-sgp4 provides the TLE manipulation and SGP4 propagation engine
- Pytest provides the testing framework
