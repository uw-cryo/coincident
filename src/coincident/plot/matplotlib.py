"""
Selection of visualizations using matplotlib.
"""

from __future__ import annotations

import warnings

import numpy as np
import pystac
import xarray as xr

try:
    import matplotlib.colors
    import matplotlib.pyplot as plt
    from matplotlib import cm
except ImportError:
    warnings.warn(
        "'matplotlib' package not found. Install for plotting functions: https://matplotlib.org/stable/install/index.html",
        stacklevel=2,
    )

from coincident._utils import depends_on_optional
from coincident.datasets.planetary_computer import WorldCover
from coincident.io.xarray import open_maxar_browse


@depends_on_optional("matplotlib")
def plot_esa_worldcover(ds: xr.Dataset) -> plt.Axes:
    """
    Plots ESA WorldCover data using Matplotlib.

    This function visualizes the ESA WorldCover dataset by mapping the class values to their recommended colors
    and descriptions.

    Parameters
    ----------
    ds : xr.Dataset
        An xarray Dataset containing the ESA WorldCover data.

    Returns
    -------
    plt.Axes
        The Matplotlib Axes object with the plot.

    References
    ----------
    https://planetarycomputer.microsoft.com/dataset/esa-worldcover#Example-Notebook
    """
    classmap = WorldCover().classmap

    colors = ["#000000" for r in range(256)]
    for key, value in classmap.items():
        colors[int(key)] = value["hex"]
    cmap = matplotlib.colors.ListedColormap(colors)

    # sequences needed for an informative colorbar
    values = list(classmap)
    boundaries = [(values[i + 1] + values[i]) / 2 for i in range(len(values) - 1)]
    boundaries = [0, *boundaries, 255]
    ticks = [
        (boundaries[i + 1] + boundaries[i]) / 2 for i in range(len(boundaries) - 1)
    ]
    tick_labels = [value["description"] for value in classmap.values()]

    fig, ax = plt.subplots()
    normalizer = matplotlib.colors.Normalize(vmin=0, vmax=255)

    da = ds.to_dataarray().squeeze()
    da.plot(ax=ax, cmap=cmap, norm=normalizer)

    colorbar = fig.colorbar(
        cm.ScalarMappable(norm=normalizer, cmap=cmap),
        boundaries=boundaries,
        values=values,
        cax=fig.axes[1].axes,
    )
    colorbar.set_ticks(ticks, labels=tick_labels)
    return ax


@depends_on_optional("matplotlib")
def plot_maxar_browse(item: pystac.Item, overview_level: int = 0) -> plt.Axes:
    """
    Plots a browse image from a STAC item using Matplotlib.

    Parameters
    ----------
    item : pystac.Item
        The STAC item containing the browse image to be plotted.
    overview_level : int, optional
        The overview level of the browse image to be opened, by default 0.

    Returns
    -------
    plt.Axes
        The Matplotlib Axes object with the plot.
    """

    da = open_maxar_browse(item, overview_level=overview_level)
    mid_lat = da.y[int(da.y.size / 2)].to_numpy()  # PD011

    fig, ax = plt.subplots(figsize=(8, 11))  # pylint: disable=unused-variable
    da.plot.imshow(rgb="band", add_labels=False, ax=ax)
    ax.set_aspect(aspect=1 / np.cos(np.deg2rad(mid_lat)))
    ax.set_title(item.id)
