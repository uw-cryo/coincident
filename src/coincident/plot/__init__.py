"""
selection of useful plots with matplotlib
"""

from __future__ import annotations

from coincident.plot.matplotlib import (
    boxplot_aspect,
    boxplot_elevation,
    boxplot_slope,
    boxplot_terrain_diff,
    compare_dems,
    hist_esa,
    plot_altimeter_points,
    plot_dem,
    plot_diff_hist,
    plot_esa_worldcover,
    plot_maxar_browse,
)
from coincident.plot.utils import (
    get_elev_diff,
    hillshade,
)

__all__ = [
    "plot_esa_worldcover",
    "plot_maxar_browse",
    "hillshade",
    "plot_dem",
    "plot_altimeter_points",
    "get_elev_diff",
    "plot_diff_hist",
    "compare_dems",
    "boxplot_terrain_diff",
    "boxplot_slope",
    "boxplot_elevation",
    "boxplot_aspect",
    "hist_esa",
]
