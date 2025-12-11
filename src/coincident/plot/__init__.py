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
    get_aspect,
    get_haversine_distance,
    get_scale,
    get_tiles,
    sample_dem_at_points,
)

__all__ = [
    "boxplot_aspect",
    "boxplot_elevation",
    "boxplot_slope",
    "boxplot_terrain_diff",
    "compare_dems",
    "get_aspect",
    "get_haversine_distance",
    "get_scale",
    "get_tiles",
    "hist_esa",
    "plot_altimeter_points",
    "plot_dem",
    "plot_diff_hist",
    "plot_esa_worldcover",
    "plot_maxar_browse",
    "sample_dem_at_points",
]
