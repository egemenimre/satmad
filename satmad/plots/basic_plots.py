# SatMAD: Satellite Mission Analysis and Design for Python
#
# Copyright (C) 2021 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.rst for more info.
"""
Basic/Common plots for quick visualisation of data.

"""
from astropy.visualization import quantity_support, time_support
from matplotlib import pyplot as plt

from satmad.plots.plotstyling import format_axis_with_grid


def plot_time_param(
    time_list, param_list, x_label=None, y_label=None, x_rotation=False
):
    """
    Plots a parameter against time.

    `Quantity` type is supported.

    Parameters
    ----------
    time_list: Time
        list of times
    param_list : ndarray or list
        list of parameter values on the y axis
    x_label : str
        X axis label
    y_label : str
        Y axis label
    x_rotation : bool
        True rotates the X axis labels by 90 deg (useful for ISOT dates)

    """

    quantity_support()
    time_support(format="isot")

    fig1, ax1 = plt.subplots(figsize=(12, 8), dpi=100)

    # default with grid
    format_axis_with_grid(ax1, x_label, y_label, x_rotation)

    # plot the data and show the plot
    ax1.plot(time_list, param_list)

    plt.show()
