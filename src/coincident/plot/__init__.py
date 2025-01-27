"""
selection of useful plots with matplotlib
"""

from __future__ import annotations

from coincident.plot.matplotlib import (
    clear_labels,
    compare_dems,
    create_hillshade,
    get_elev_diff,
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
    "create_hillshade",
    "clear_labels",
    "plot_dem",
    "plot_altimeter_points",
    "get_elev_diff",
    "plot_diff_hist",
    "plot_worldcover_custom_ax",
    "compare_dems",
]
