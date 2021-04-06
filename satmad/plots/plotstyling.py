# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Useful plot styles.

"""
from matplotlib.ticker import AutoLocator, AutoMinorLocator


def format_axis_with_grid(ax, x_label=None, y_label=None, x_rotation=False):
    """
    Adds a major and minor grid to axes, adds labels.

    Parameters
    ----------
    ax : Axes
        Axes object containing the figure axis (modified)
    x_label : str
        X axis label
    y_label : str
        Y axis label
    x_rotation : bool
        True rotates the X axis labels by 90 deg (useful for ISOT dates)


    Returns
    -------
    Axes
        Axes object containing the figure axis (modified)
    """

    # Change major ticks
    ax.xaxis.set_major_locator(AutoLocator())
    ax.yaxis.set_major_locator(AutoLocator())

    # Change minor ticks to show with divisor of 4
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.grid(which="major", color="#CCCCCC", linestyle="--")
    ax.grid(which="minor", color="#CCCCCC", linestyle=":")

    # Set axis labels
    if x_label:
        ax.set_xlabel(x_label)

    if y_label:
        ax.set_ylabel(y_label)

    # Rotate x axis labels
    if x_rotation:
        ax.tick_params(axis="x", labelrotation=90)

    return ax
