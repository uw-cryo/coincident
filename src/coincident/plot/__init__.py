"""
selection of useful plots with matplotlib
"""

from __future__ import annotations

from coincident.plot.matplotlib import (
    boxplot_aspect,
    boxplot_elevation,
    boxplot_slope,
    boxplot_terrain_diff,
    clear_labels,
    compare_dems,
    create_stats_legend,
    get_elev_diff,
    hillshade,
    hist_esa,
    plot_altimeter_points,
    plot_dem,
    plot_diff_hist,
    plot_esa_worldcover,
    plot_maxar_browse,
    plot_worldcover_custom_ax,
)

__all__ = [
    "plot_esa_worldcover",
    "plot_maxar_browse",
    "hillshade",
    "clear_labels",
    "plot_dem",
    "plot_altimeter_points",
    "get_elev_diff",
    "plot_diff_hist",
    "plot_worldcover_custom_ax",
    "compare_dems",
    "boxplot_terrain_diff",
    "boxplot_slope",
    "boxplot_elevation",
    "boxplot_aspect",
    "create_stats_legend",
    "hist_esa",
]
