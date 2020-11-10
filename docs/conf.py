# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------

project = "satmad"
copyright = "2020, Egemen Imre"
author = "Egemen Imre"

# Version Info
# ------------
# The short X.Y version.
version = "0.0.6"
# The full version, including alpha/beta/rc tags.
release = "0.0.6"

# -- General configuration ---------------------------------------------------
# By default, highlight as Python 3.
highlight_language = "python3"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",  # auto api generation
    "sphinx.ext.napoleon",  # numpy support
    "sphinx.ext.mathjax",  # LaTex style math
    "sphinx.ext.graphviz",  # Dependency diagrams
    "sphinx.ext.intersphinx",  # Link mapping to external projects
    # "sphinx.ext.doctest",  # Doctest
    "nbsphinx",  # Jupyter notebook support
]
numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The suffix of source filenames.
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# Intersphinx configuration
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    "matplotlib": ("https://matplotlib.org", None),
}

if os.environ.get("READTHEDOCS") == "True":
    nbsphinx_execute = "never"
else:
    nbsphinx_execute = "always"

# Controls when a cell will time out (defaults to 30; use -1 for no timeout):
nbsphinx_timeout = 1000

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
