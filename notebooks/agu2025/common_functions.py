"""
Common plots for each PCD site
"""

from __future__ import annotations

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import contextily as ctx
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from folium import Element
from folium.plugins import MiniMap
from matplotlib.patches import Patch
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry import box

from coincident.overlaps import geographic_area
from coincident.plot.utils import get_scale

# always use the same colors across plots
sensor2color = {
    "ALS": "blue",
    "Maxar": "gold",
    "ICESat-2": "green",
    "GEDI": "red",
    "Overlap Area": "black",
}


def load_geodataframes(site_meta):
    return (
        site_meta["als"],
        site_meta["maxar"],
        site_meta["is2"],
        site_meta["gedi"],
        site_meta["overlap"],
    )


def interactive_site_map(site_meta, title=None):
    gf_als, gf_maxar, gf_is2, gf_gedi, gf_overlap = load_geodataframes(site_meta)

    # first plot sets zoom level and tiles
    style_args = {"fillOpacity": 0.15, "weight": 2.5}

    m = gf_als.explore(
        name="ALS", color="gray", style_kwds=style_args, tiles="Esri.WorldImagery"
    )

    gf = gpd.pd.concat(
        [gf_als, gf_maxar, gf_is2, gf_gedi], join="inner", ignore_index=True
    )
    m = gf.explore(m=m, column="collection", title=title)
    m = gf_overlap.explore(
        m=m, name="Overlap Area", color="black", style_kwds=style_args
    )

    # Overview of location
    minimap = MiniMap(
        position="topright",  # Default position
        zoom_level_offset=-7,
        # tiles="OpenStreetMap",
        # tiles="Esri.WorldImagery",
        toggle_display=True,  # Makes it collapsible
    )
    minimap.add_to(m)

    if title:
        title_html = f"""
<div style="position: fixed; top: 10px; left: 50%; transform: translateX(-50%); z-index:1000; font-size:16px; font-weight:bold; background-color: white; padding: 5px; border-radius: 5px;">
    {title}
</div>
"""
        # title_html = f'<h3 align="center" style="font-size:20px"><b>{title}</b></h3>'
        m.get_root().html.add_child(Element(title_html))

    return m


def static_site_map(site_meta, title=None):
    _, ax = plt.subplots(figsize=(8, 8))

    gf_als, gf_maxar, gf_is2, gf_gedi, gf_overlap = load_geodataframes(site_meta)

    area_sqkm = geographic_area(gf_overlap) * 1e-6

    gf_is2_i = gf_als.overlay(gf_is2, how="intersection")
    gf_gedi_i = gf_als.overlay(gf_gedi, how="intersection")

    alpha = 0.3
    lw = 1
    legend_handles = [
        Patch(
            facecolor=sensor2color["ALS"],
            edgecolor=sensor2color["ALS"],
            alpha=alpha,
            linewidth=lw,
        ),
        Patch(
            facecolor=sensor2color["Maxar"],
            edgecolor=sensor2color["Maxar"],
            alpha=alpha,
            linewidth=lw,
        ),
        Patch(
            facecolor=sensor2color["ICESat-2"],
            edgecolor=sensor2color["ICESat-2"],
            alpha=alpha,
            linewidth=lw,
        ),
        Patch(
            facecolor=sensor2color["GEDI"],
            edgecolor=sensor2color["GEDI"],
            alpha=alpha,
            linewidth=lw,
        ),
        Patch(facecolor="none", edgecolor=sensor2color["Overlap Area"], linewidth=lw),
    ]
    legend_labels = [
        "ALS",
        "Maxar Stereo",
        "ICESat-2",
        "GEDI",
        f"Overlap Area ({int(area_sqkm.to_numpy()[0])} km²)",
    ]

    gf_als.plot(
        ax=ax,
        facecolor=sensor2color["ALS"],
        edgecolor=sensor2color["ALS"],
        alpha=alpha,
        linewidth=lw,
    )
    gf_maxar.plot(
        ax=ax,
        facecolor=sensor2color["Maxar"],
        edgecolor=sensor2color["Maxar"],
        alpha=alpha,
        linewidth=lw,
    )
    gf_is2_i.plot(
        ax=ax,
        facecolor=sensor2color["ICESat-2"],
        edgecolor=sensor2color["ICESat-2"],
        alpha=alpha,
        linewidth=lw,
    )
    gf_gedi_i.plot(
        ax=ax,
        facecolor=sensor2color["GEDI"],
        edgecolor=sensor2color["GEDI"],
        alpha=alpha,
        linewidth=lw,
    )
    gf_overlap.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=lw)

    # suppress unequal aspect ratio warning
    scalebar = ScaleBar(
        get_scale(ax), units="m", location="lower right", rotation="horizontal-only"
    )
    ax.add_artist(scalebar)

    ctx.add_basemap(ax, crs=gf_als.crs, source=ctx.providers.Esri.WorldImagery)

    # Add inset map showing site location within United States
    # Create inset axes in upper right corner
    axins = ax.inset_axes([0.01, 0.77, 0.3, 0.3], projection=ccrs.PlateCarree())
    # axins = ax.inset_axes([0.65, 0.65, 0.3, 0.3], projection=ccrs.PlateCarree())

    # Add US map features
    axins.add_feature(cfeature.LAND, facecolor="lightgray")
    axins.add_feature(cfeature.COASTLINE, linewidth=0.5)
    axins.add_feature(cfeature.BORDERS, linewidth=0.5)
    axins.add_feature(cfeature.STATES, linewidth=0.3, edgecolor="gray")

    # Set US extent
    axins.set_extent([-125, -66, 24, 50], crs=ccrs.PlateCarree())

    # Add rectangle showing the site bounds
    main_bounds = ax.get_xlim() + ax.get_ylim()  # (xmin, xmax, ymin, ymax)
    # Plot bounding box
    site_bbox = box(main_bounds[0], main_bounds[2], main_bounds[1], main_bounds[3])
    axins.add_geometries(
        [site_bbox],
        ccrs.PlateCarree(),
        facecolor="none",
        edgecolor="red",
        linewidth=1,
        zorder=10,
    )

    # Get centroid of the site in lat/lon
    # centroid = gf_overlap.to_crs(epsg=4326).geometry.centroid.iloc[0]
    # # Plot site location as a star
    # axins.plot(centroid.x, centroid.y, marker='*', color='red',
    #             markersize=3, transform=ccrs.PlateCarree(), zorder=10)

    ax.legend(legend_handles, legend_labels)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)
    plt.tight_layout()
    plt.show()
    return ax


def plot_timeline(site_meta, searchdates, title=None, neon=False):
    fig, ax = plt.subplots(figsize=(8.5, 3))
    gf_als, gf_maxar, gf_is2, gf_gedi, _ = load_geodataframes(site_meta)

    # Fix for NEON
    if neon:
        gf_als["collection"] = "NEON"
        gf_als["start_datetime"] = pd.to_datetime(gf_als["start_datetime"])
        gf_als["end_datetime"] = pd.to_datetime(gf_als["end_datetime"])

    als_start = gf_als["start_datetime"].iloc[0]
    als_end = gf_als["end_datetime"].iloc[0]
    # ax.axvspan(als_start, als_end, color="black", alpha=0.5, label="ALS Window")

    # Add date annotations for ALS start and end
    annotation_fontsize = 9
    annotation_fontweight = "bold"
    rotation = 90
    ax.text(
        als_start,
        ax.get_ylim()[1],
        als_start.strftime("%Y-%m-%d"),
        rotation=rotation,
        verticalalignment="top",
        horizontalalignment="right",
        fontsize=annotation_fontsize,
        fontweight=annotation_fontweight,
    )
    ax.text(
        als_end,
        ax.get_ylim()[1],
        als_end.strftime("%Y-%m-%d"),
        rotation=rotation,
        verticalalignment="top",
        horizontalalignment="right",
        fontsize=annotation_fontsize,
        fontweight=annotation_fontweight,
    )
    # Add date annotations for Maxar points
    for date in gf_maxar["datetime"]:
        ax.text(
            date,
            "Maxar",
            f"  {date.strftime('%Y-%m-%d')}",
            rotation=rotation,
            verticalalignment="top",
            fontsize=annotation_fontsize,
            fontweight=annotation_fontweight,
        )

    # Add date annotations for ICESat-2 points
    for date in gf_is2["datetime"]:
        ax.text(
            date,
            "ICESat-2",
            f"  {date.strftime('%Y-%m-%d')}",
            rotation=rotation,
            verticalalignment="top",
            fontsize=annotation_fontsize,
            fontweight=annotation_fontweight,
        )

    # Add date annotations for GEDI points
    for date in gf_gedi["datetime"]:
        ax.text(
            date,
            "GEDI",
            f"  {date.strftime('%Y-%m-%d')}",
            rotation=rotation,
            verticalalignment="top",
            fontsize=annotation_fontsize,
            fontweight=annotation_fontweight,
        )

    ax.axvline(
        als_start,
        color=sensor2color["ALS"],
        linestyle="--",
        label="ALS Acquisition Span",
    )
    ax.axvline(als_end, color=sensor2color["ALS"], linestyle="--")

    # Also include range of search for coincident data
    ax.axvline(
        pd.to_datetime(searchdates[0]),
        color="gray",
        linestyle="-",
        label="Search Window",
    )
    ax.axvline(pd.to_datetime(searchdates[1]), color="gray", linestyle="-")

    ax.scatter(
        gf_maxar["datetime"],
        ["Maxar"] * len(gf_maxar),
        marker="s",
        s=80,
        label="Maxar Stereo",
        color=sensor2color["Maxar"],
    )
    ax.scatter(
        gf_is2["datetime"],
        ["ICESat-2"] * len(gf_is2),
        marker="D",
        s=80,
        label="ICESat-2",
        color=sensor2color["ICESat-2"],
    )
    ax.scatter(
        gf_gedi["datetime"],
        ["GEDI"] * len(gf_gedi),
        marker="D",
        s=80,
        label="GEDI",
        color=sensor2color["GEDI"],
    )

    ax.set_title(title)
    # ax.set_xlabel("Date")
    # ax.set_ylabel("Data Collection")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    # ax.legend(loc="best")
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    fig.autofmt_xdate()
    plt.tight_layout()
