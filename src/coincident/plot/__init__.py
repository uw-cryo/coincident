"""
selection of useful plots with matplotlib
"""

from __future__ import annotations

from coincident.plot.matplotlib import (
    clear_labels,
    compare_dems,
    get_elev_diff,
    hillshade,
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
]
