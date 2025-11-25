"""
Selection of visualizations using matplotlib.
"""

from __future__ import annotations

from typing import Any

import contextily as ctx
import geopandas as gpd
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pystac
import xarray as xr
from matplotlib import cm
from matplotlib.ticker import FuncFormatter
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry import box

from coincident.datasets.planetary_computer import WorldCover
from coincident.io.xarray import open_maxar_browse, to_dataset
from coincident.plot import utils
from coincident.search import search


def plot_esa_worldcover(
    ds: xr.Dataset,
    ax: plt.Axes | None = None,
    cax: plt.Axes | None = None,
    add_colorbar: bool = True,
    add_labels: bool = True,
) -> plt.Axes:
    """
    Map view of ESA WorldCover data.

    This function visualizes the ESA WorldCover dataset by mapping the class values to their recommended colors
    and descriptions.

    Parameters
    ----------
    ds : xr.Dataset
        An xarray Dataset containing the ESA WorldCover data.
    ax : plt.Axes
        The Matplotlib Axes on which the plot will be drawn.
    cax : plt.Axes
        The Matplotlib Axes used for the colorbar.
    add_colorbar : bool, optional
        Whether to include a colorbar in the plot, by default True.
    add_labels : bool, optional
        Whether to add labels to the plot, by default True.

    Returns
    -------
    plt.Axes
        The Matplotlib Axes object with the plot.

    Notes
    -----
    https://planetarycomputer.microsoft.com/dataset/esa-worldcover#Example-Notebook
    """
    # Custom categorical colormap for ESA WorldCover
    classmap = WorldCover().classmap

    colors = ["#000000" for r in range(256)]
    for key, value in classmap.items():
        colors[int(key)] = value["hex"]
    cmap = matplotlib.colors.ListedColormap(colors)
    values = list(classmap)
    boundaries = [(values[i + 1] + values[i]) / 2 for i in range(len(values) - 1)]
    boundaries = [0, *boundaries, 255]
    ticks = [
        (boundaries[i + 1] + boundaries[i]) / 2 for i in range(len(boundaries) - 1)
    ]
    tick_labels = [value["description"] for value in classmap.values()]
    normalizer = matplotlib.colors.Normalize(vmin=0, vmax=255)

    da = ds.to_dataarray().squeeze()

    if ax is None:
        fig, ax = plt.subplots()
        # set aspect ratio according to mid scene latitude
        if da.rio.crs.to_epsg() == 4326:
            mid_lat = da.latitude[int(da.latitude.size / 2)].to_numpy()  # PD011
            ax.set_aspect(aspect=utils.get_lonlat_aspect(mid_lat))
    else:
        fig = ax.get_figure()

    da.plot(
        ax=ax, cmap=cmap, norm=normalizer, add_labels=add_labels, add_colorbar=False
    )

    if add_colorbar:
        # Specific to panel_plots, so rather than cax, add_to_panel=True boolean?
        if cax is not None:
            # Override standard xarray colorbar with custom one
            colorbar = fig.colorbar(
                cm.ScalarMappable(norm=normalizer, cmap=cmap),
                boundaries=boundaries,
                values=values,
                ax=cax,
                orientation="horizontal",
                pad=-1.1,  # hack to get colorbar in subplot mosaic in the right place
            )
            colorbar.set_ticks(ticks, labels=tick_labels, rotation=90)
        else:
            colorbar = fig.colorbar(
                cm.ScalarMappable(norm=normalizer, cmap=cmap),
                boundaries=boundaries,
                values=values,
                ax=ax,
            )
            colorbar.set_ticks(ticks, labels=tick_labels)

    return ax


def plot_maxar_browse(
    item: pystac.Item, ax: plt.Axes | None = None, overview_level: int = 0
) -> plt.Axes:
    """
    Map view of Maxar browse image from a STAC item using Matplotlib.

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

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 11))
        mid_lat = da.latitude[int(da.latitude.size / 2)].to_numpy()  # PD011
        ax.set_aspect(aspect=utils.get_lonlat_aspect(mid_lat))
        ax.set_title(item.id)

    nbands = da.band.size
    if nbands == 1:
        da_plot = da.squeeze("band")
        da_plot.plot.imshow(add_labels=False, ax=ax, cmap="gray", add_colorbar=False)
    elif nbands == 3:
        da.plot.imshow(rgb="band", add_labels=False, ax=ax, add_colorbar=False)
    else:
        error_message = (
            f"Maxar browse image must have 1 or 3 bands. Found: {nbands} bands."
        )
        raise ValueError(error_message)

    return ax


def _clear_labels(ax: plt.Axes) -> None:
    """
    Remove x and y axis labels.

    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes
        The matplotlib axis object to remove labels from

    Returns
    -------
    None
        This function modifies the input axis object directly
    """
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.set_xlabel("")
    ax.set_ylabel("")


def plot_dem(
    da: xr.DataArray,
    ax: plt.Axes | None = None,
    cmap: str = "inferno",
    title: str = "",
    da_hillshade: xr.DataArray | None = None,
    **kwargs: Any,
) -> plt.Axes:
    """
    Map view of DEM with an option to plot a hillshade underneath

    Parameters
    ----------
    da: xarray.core.dataarray.DataArray
        A 2d xarray DataArray containing elevation data
    ax: matplotlib.axes._axes.Axes
        Axis on which to plot the DEM
    cmap: str
        Matplotlib colormap to use for the DEM
    title: str
        Title of the plot
    da_hillshade: xarray.core.dataarray.DataArray | None
        A 2d xarray DataArray containing hillshade data
    cmin: float | None
        Minimum value for the elevation colorbar
    cmax: float | None
        Maximum value for the elevation colorbar
    **kwargs: Any
        Artist properties passed to matplotlib.axes.Axes.imshow (e.g., vmin, vmax, etc.)

    Returns
    -------
    plt.Axes
        The Axes object with the DEM plot

    Notes
    -----
    alpha=0.5 set by default if you pass in a hillshade.
    """
    if ax is None:
        _, ax = plt.subplots()

    # hillshade with alpha=1.0
    if isinstance(da_hillshade, xr.DataArray):
        da_hillshade.plot.imshow(
            ax=ax, cmap="gray", alpha=1.0, add_colorbar=False, add_labels=False
        )
        if "alpha" not in kwargs:
            kwargs["alpha"] = 0.5
    # plot the DEM with alpha=0.5
    da.squeeze().plot.imshow(ax=ax, cmap=cmap, **kwargs)
    ax.set_title(title)

    return ax


def plot_altimeter_points(
    gf: gpd.GeoDataFrame,
    column: str,
    ax: plt.Axes | None = None,
    cmap: str = "inferno",
    title: str = "",
    da_hillshade: xr.DataArray | None = None,
    **kwargs: Any,
) -> plt.Axes:
    """
    Map view of laser altimeter point data with an optional hillshade background.

    Parameters
    ----------
    gf: geopandas.GeoDataFrame
        GeoDataFrame containing altimeter point data
    ax: matplotlib.axes._axes.Axes
        Axis on which to plot the points
    column: str
        Name of the GeoDataFrame elevation column to plot (e.g. 'h_li' for ICESat-2 ATL06)
    cmap: str
        Matplotlib colormap to use for the points
    alpha_p: float
        the alpha of the points
    title: str
        Title of the plot
    da_hillshade: xarray.core.dataarray.DataArray | None
        A 2d xarray DataArray containing hillshade data
    p_alpha: float | None
        Alpha value for the points, default is 0.5
    **kwargs: Any
        Artist properties passed to geopandas.GeoDataFrame.plot (e.g., markersize, legend, etc.)

    Returns
    -------
    plt.Axes
        The Axes object with the altimeter points plot

    Notes
    -----
    alpha=0.5 set by default if you pass in a hillshade.
    """
    if ax is None:
        _, ax = plt.subplots()

    # Plot the hillshade if available
    if isinstance(da_hillshade, xr.DataArray):
        da_hillshade.plot.imshow(
            ax=ax,
            cmap="gray",
            alpha=1.0,
            add_colorbar=False,
            add_labels=False,
        )
        if "alpha" not in kwargs:
            kwargs["alpha"] = 0.5
    # NOTE: geopandas automatically scales aspect ratio
    gf.plot(ax=ax, column=column, cmap=cmap, **kwargs)
    ax.set_title(title)
    return ax


def plot_diff_hist(
    s: pd.Series,
    ax: plt.Axes | None = None,
    range: tuple[float, float] | None = None,
    add_stats: bool = True,
    **kwargs: Any,
) -> plt.Axes:
    """
    Histogram of elevation differences with median and mean lines.

    Parameters
    ----------
    s: pd.Series
        Pandas Series containing elevation differences to plot
    ax: plt.Axes
        Axis on which to plot the histogram
    range: (float, float) | None
        Minimum and maximum values for histogram range
    add_stats: bool
        Whether to add annotations and legend with statistics to the plot, default True
    **kwargs: Any
        Artist properties passed to ax.hist()

    Returns
    -------
    plt.Axes
        The Axes object with the histogram plot
    """

    if ax is None:
        _, ax = plt.subplots()
        ax.set_xlabel("Elevation difference (m)")

    s.hist(ax=ax, bins=100, range=range, color="gray", **kwargs)
    ax.axvline(0, color="k", lw=1)

    if add_stats:
        stats_dict = utils.calc_stats(s)

        median = stats_dict["Median"]
        mean = stats_dict["Mean"]

        # NOTE: hardcoded unit label of meters
        ax.axvline(mean, label=f"{'Mean':<7}{mean:>5.2f}", color="magenta", lw=0.5)
        ax.axvline(median, label=f"{'Median':<7}{median:>5.2f}", color="cyan", lw=0.5)

        _create_stats_legend(ax, stats_dict)

    ax.grid(False)

    def _simple_label(x: float, pos: int) -> str:
        return f"{x:.1e}".replace("e+0", "e") if pos > 0 else f"{x:.0f}"

    yticks = ax.get_yticks()
    if yticks[1] > 1000:
        ax.yaxis.set_major_formatter(FuncFormatter(_simple_label))

    return ax


def compare_dems(
    dem_dict: dict[str, xr.Dataset],
    gdf_dict: dict[str, tuple[gpd.GeoDataFrame, str]] | None = None,
    add_hillshade: bool = False,
    altimetry_basemap: str | None = None,
    elevation_cmap: str = "plasma",
    elevation_clim: tuple[float, float] | None = None,
    diff_cmap: str = "RdBu",
    diff_clim: tuple[float, float] = (-10, 10),
    figsize: tuple[float, float] | None = None,
) -> dict[str, plt.Axes]:
    """
    Create a panel figure comparing DEM and altimetry elevations.

    The first row shows elevation maps and altimeter points over first DEM hillshade (if basemap=='hillshade')
    The second row shows elevation differences against the first DEM, plus ESA Worldcover for context
    The third row shows the histograms of the elevation differences

    Parameters
    ----------
    dem_dict: dict
        Dictionary of xr.Datasets. Each Dataset should have .elevation and .hillshade variables.
        The first DEM in this list will be the 'source' DEM and will served as a reference for differencing
        e.g.  {'dem_1':ds1,'dem_2':ds2,'dem_3':d3} will result in diff_1 = dem_2 - dem_1 and diff_2 = dem_3 - dem_1
    gdf_dict: Optional Dict[str, Tuple[gpd.GeoDataFrame, str]]
        Dictionary where keys are subplot names and values are (GeoDataFrame, column_name) pairs.
        column_name denotes elevation values to plot (e.g. h_li for ICESat-2 ATL06)
    add_hillshade: Optional bool
        Overlay elevation data on hillshade. If True, your DEM Dataset must share a 'hillshade' variable.
    altimetery_basemap: Optional str
        If 'hillshade', your reference DEM needs a 'hillshade' variable. Or pass an xyzservices provider name like 'Esri.WorldImagery'
    elevation_cmap: Optional str
        Colormap for elevation. Default: 'inferno'
    elevation_clim: Optional Tuple[float, float]
        Tuple for elevation color limits, otherwise scaled to 2nd and 98th percentiles
    diff_cmap: Optional str
        Colormap for elevation differences. Default: 'RdBu'
    diff_clim: Optional Tuple[float, float]
        Tuple for difference color limits. Default (-10,10)
    figsize: Optional Tuple[int, int]
        Figure size (width, height) tuple.

    Returns
    ----------
        Dictionary containing the axes objects for each subplot

    Notes
    -----
    All inputs assumed to have the same CRS and to be aligned
    All DEMs must have variable 'elevation'
    A minimum requirement of 2 DEMs
    A maximum of 5 total datasets (dem_list + gdf_dict) to be passed
    """
    # FIGURE SETUP
    # ---------------------------
    # TODO: check for consistent CRS?
    # Calculate number of columns for the plot (DEMs + GeoDataFrames if provided)
    n_dems = len(dem_dict)
    # for convenience
    dem_list = list(dem_dict.values())
    dem_names = list(dem_dict.keys())
    first_dem = dem_list[0]
    gdf_keys = list(gdf_dict.keys()) if gdf_dict else []
    n_columns = max(2, n_dems) + len(gdf_keys)
    n_columns = min(n_columns, 5)

    # Scale figure size with number of columns
    if figsize is None:
        figsize = (3 * n_columns, 11)

    mosaic = [
        [f"dem_{i}" for i in range(n_dems)] + [f"dem_{g}" for g in gdf_keys],
        ["worldcover"]
        + [f"diff_{i}" for i in range(1, n_dems)]
        + [f"diff_{g}" for g in gdf_keys],
        ["wc_legend"]
        + [f"hist_{i}" for i in range(1, n_dems)]
        + [f"hist_{g}" for g in gdf_keys],
    ]
    # All columns get same width, but first two rows have 2x height
    gs_kw = {"width_ratios": [1] * len(mosaic[0]), "height_ratios": [2, 2, 1]}
    fig, axd = plt.subplot_mosaic(
        mosaic, gridspec_kw=gs_kw, figsize=figsize, layout="constrained"
    )

    # Determine aspect ratio
    if first_dem.rio.crs.is_geographic:
        # coordinate name: y or latitude?
        mid_lat = first_dem.y[int(first_dem.y.size / 2)].to_numpy()  # PD011
        aspect = utils.get_lonlat_aspect(mid_lat)
    elif first_dem.rio.crs.is_projected:
        aspect = 1.0
    else:
        aspect = None  # Not sure how matplotlib does default aspect...

    reference_ax = axd["dem_0"]

    for key in axd:
        if key.startswith(("dem_", "diff_", "worldcover")):
            axd[key].sharex(reference_ax)
            axd[key].sharey(reference_ax)

    def _style_ax(
        ax: plt.Axes,
        facecolor: str | None = None,
        aspect: float | None = None,
        altimetry_basemap: str | None = None,
    ) -> None:
        """Helper function to style axes (remove ticks and set face color)."""
        ax.set_xticks([])
        ax.set_yticks([])
        if facecolor:
            ax.set_facecolor(facecolor)
        if aspect:
            ax.set_aspect(aspect)
        if altimetry_basemap and altimetry_basemap != "hillshade":
            ctx.add_basemap(
                ax,
                crs=first_dem.rio.crs,
                attribution=False,
                source=utils.get_tiles(altimetry_basemap),
            )

    # Get Landcover for sceene
    bounds = first_dem.rio.bounds()
    gf_search = gpd.GeoDataFrame(
        geometry=[box(*bounds)],
        crs=first_dem.rio.crs,
    )
    gf_wc = search(
        dataset="worldcover",
        intersects=gf_search,
        datetime=["2021"],
        # NOTE: 2020 throwing an error...
    )
    # TODO: set requested resolution based on AOI / DEM resolution?
    # Or always keep at native 10m resolution ?
    ds_wc = to_dataset(
        gf_wc,
        bands=["map"],
        aoi=gf_search,
    ).compute()

    # Calculate colorbar limits
    if elevation_clim is None:
        vmin = np.nanpercentile(first_dem.elevation.to_numpy(), 2)
        vmax = np.nanpercentile(first_dem.elevation.to_numpy(), 98)
    else:
        vmin, vmax = elevation_clim

    # TOP ROW: Elevation maps
    # ---------------------------
    # raster elevs
    for i, (title, dem) in enumerate(dem_dict.items()):
        axi = axd[f"dem_{i}"]
        plot_dem(
            da=dem.elevation,
            ax=axi,
            cmap=elevation_cmap,
            da_hillshade=dem.hillshade if add_hillshade else None,
            add_colorbar=False,
            vmin=vmin,
            vmax=vmax,
            title=title,
        )
        _clear_labels(axi)
        _style_ax(axi, aspect=aspect)
        # Add a scalebar to first plot
        if i == 0:
            rotation = "horizontal-only" if aspect != 0 else None
            scalebar = ScaleBar(
                utils.get_scale(axi),
                units="m",
                location="lower right",
                rotation=rotation,
            )
            axi.add_artist(scalebar)
            # Keep track of mappable for row colorbar
            elevation_mappable = axi.images[1] if add_hillshade else axi.images[0]

    # point elevs
    # need this for MyPy to handle the chance the gdf_dict is None
    # despite the logic in the loop "if gdf_dict else []"
    if gdf_dict:
        # Plot GeoDataFrames (e.g., IS2 and GEDI)
        for key, (gdf, column) in gdf_dict.items() if gdf_dict else []:
            axi = axd[f"dem_{key}"]
            plot_altimeter_points(
                gdf,
                column,
                ax=axi,
                title=key,
                da_hillshade=first_dem.hillshade
                if altimetry_basemap == "hillshade"
                else None,
                cmap=elevation_cmap,
                vmin=vmin,
                vmax=vmax,
                s=(figsize[0] * figsize[1])
                / n_columns**2,  # Dynamic point size based on figure dimensions
            )
            _style_ax(
                axi,
                aspect=aspect,
                facecolor="black",
                altimetry_basemap=altimetry_basemap,
            )

    # Add colorbar to last plot in top row
    # https://matplotlib.org/stable/users/explain/axes/colorbar_placement.html#colorbars-attached-to-fixed-aspect-ratio-axes
    cax = axi.inset_axes([1.04, 0.1, 0.1, 0.8])  # [x0, y0, width, height]
    fig.colorbar(elevation_mappable, cax=cax, label="Elevation (m)")

    # MIDDLE ROW: Differences
    # ---------------------------
    plot_esa_worldcover(
        ds_wc.isel(time=0),
        ax=axd["worldcover"],
        cax=axd["wc_legend"],
        add_colorbar=True,
        add_labels=False,
    )
    _style_ax(axd["worldcover"], aspect=aspect)
    axd["worldcover"].set_title("ESA WorldCover 2021")

    # raster diffs
    for i, dem in enumerate(dem_list[1:], start=1):
        diff = utils.get_elev_diff(dem, first_dem)
        axi = axd[f"diff_{i}"]
        diff.elev_diff.squeeze().plot.imshow(
            ax=axi,
            cmap=diff_cmap,
            add_colorbar=False,
            add_labels=False,
            vmin=diff_clim[0],
            vmax=diff_clim[1],
        )
        _style_ax(axi, aspect=aspect, facecolor="black")
        axi.set_title(f"{dem_names[i]} minus {dem_names[0]}")
        if i == 1:
            # keep track of mappable for colorbar
            diff_mappable = axi.images[0]

    if gdf_dict:
        for _, (key, (gdf, column)) in enumerate(gdf_dict.items() if gdf_dict else []):
            axi = axd[f"diff_{key}"]
            gf_diff = utils.get_elev_diff(gdf, first_dem, source_col=column)
            plot_altimeter_points(
                gf_diff,
                "elev_diff",
                axi,
                da_hillshade=first_dem.hillshade
                if altimetry_basemap == "hillshade"
                else None,
                cmap=diff_cmap,
                vmin=diff_clim[0],
                vmax=diff_clim[1],
                title=f"{column} minus {dem_names[0]}",
                s=(figsize[0] * figsize[1])
                / n_columns**2,  # Dynamic sizing based on figure dimensions
            )
            _style_ax(
                axi,
                aspect=aspect,
                facecolor="black",
                altimetry_basemap=altimetry_basemap,
            )

    # Add colorbar to last plot in middle row
    cax = axi.inset_axes([1.04, 0.1, 0.1, 0.8])  # [x0, y0, width, height]
    fig.colorbar(
        diff_mappable,
        cax=cax,
        label="Elevation Difference (m)",
    )

    # BOTTOM ROW: Histograms
    # ---------------------------
    axd["wc_legend"].axis("off")

    # raster hists
    for i, dem in enumerate(dem_list[1:], start=1):
        axi = axd[f"hist_{i}"]
        diff = utils.get_elev_diff(dem, first_dem).elev_diff.to_series()
        plot_diff_hist(
            diff,
            ax=axi,
            range=(diff_clim[0], diff_clim[1]),
            add_stats=False,
        )
        axi.set_xlabel("Elevation Difference (m)")
        stats = utils.calc_stats(diff)
        stats_title = rf"$\mu$={stats['Mean']:.2f}, $\sigma$={stats['Std']:.2f}, $n$={stats['Count']:.1e}".replace(
            "e+0", "e"
        )
        axi.set_title(stats_title, fontsize=8)

    # gdf hists
    if gdf_dict:
        for _, (key, (gdf, column)) in enumerate(gdf_dict.items() if gdf_dict else []):
            axi = axd[f"hist_{key}"]
            diff = utils.get_elev_diff(gdf, first_dem, source_col=column).elev_diff
            plot_diff_hist(
                diff,
                ax=axi,
                range=(diff_clim[0], diff_clim[1]),
                add_stats=False,
            )
            stats = utils.calc_stats(diff)
            stats_title = rf"$\mu$={stats['Mean']:.2f}, $\sigma$={stats['Std']:.2f}, $n$={stats['Count']:.1e}".replace(
                "e+0", "e"
            )
            axi.set_title(stats_title, fontsize=8)

    return axd  # type: ignore[no-any-return]


def boxplot_terrain_diff(
    dem_list: list[xr.Dataset | gpd.GeoDataFrame],
    ax: plt.Axes | None = None,
    terrain_v: str = "slope",
    terrain_groups: np.ndarray | None = None,
    ylims: list[float] | None = None,
    title: str = "Elevation Differences (m) by Terrain",
    ylabel: str = "Elevation Difference (m)",
    elev_col: str | None = None,
) -> plt.Axes:
    """
    Box plots of elevation differences by terrain group for two input elevation datasets.
    (e.g. elevation differences over varying slopes for two DEMs)
    This also shows the counts of the grouped elevation differences

    Parameters
    ----------
    dem_list : list[xr.Dataset | gpd.GeoDataFrame]
        List containing exactly two elevation datasets to compare. The first
        dataset must be a xr.Dataset with an 'elevation' variable and a variable
        corresponding to the terrain_v variable (e.g. 'slope').
    terrain_v : str
        Terrain variable of interest (e.g. 'slope') that exists in first DEM. This
        is what your elevation differences will be grouped by
    terrain_groups : np.ndarray
        Array defining the edges of terrain groups. e.g. if set to np.arange(0, 46, 1),
        the groups will be [0, 1, 2, ..., 45]
    ylims : list
        The y-axis limits for the plot. If None, then the y-axis limits will be
        automatically adjusted to fit the whiskers of each boxplot
    title : str
        The title of the plot. Default is 'Elevation Differences (m) by Terrain'
    ylabel : str
        The yabel of the plot. Default is 'Elevation Difference (m)'
    elev_col : str
        The name of the column containing elevation values to difference in
        the GeoDataFrame if a GeoDataFrame is provided

    Returns
    -------
    plt.Axes
        The matplotlib axes object containing the plot.

    Notes
    -----
    This function assumes that the datasets in dem_list are in the same CRS and aligned
    This function can also work with point geometry GeoDataFrames (e.g. ICESat-2 points)
    If using a GeoDataFrame, you must also provide the elev_col parameter which is the column name containing elevation values you wish to compare
    This function requires there to be EXACTLY 2 datasets in the dem_list.
    The first dataset in dem_list MUST be a xr.Dataset with an 'elevation' variable and the corresponding terrain_v variable (e.g. 'slope')
    """
    if len(dem_list) != 2:
        msg_len = "dem_list must contain exactly two datasets"
        raise ValueError(msg_len)
    # ruff B008 if default is set to np.arange(0, 46, 1)
    if terrain_groups is None:
        msg_groups = "terrain_groups is undefined"
        raise ValueError(msg_groups)
    # difference (second - first) based on type (xr.dataset vs gdf)
    if isinstance(dem_list[1], gpd.GeoDataFrame):
        if elev_col is None:
            msg_col_arg = "elev_col must be provided when using point data"
            raise ValueError(msg_col_arg)
        da_points = dem_list[1].get_coordinates().to_xarray()
        samples = (
            dem_list[0]
            .interp(da_points)
            .drop_vars(["band", "spatial_ref"])
            .to_dataframe()
        )
        dem_list[1]["elev_diff"] = (
            dem_list[1][elev_col].to_numpy() - samples["elevation"].to_numpy()
        )
        dem_list[1][terrain_v] = samples[terrain_v].to_numpy()
        diff_data = dem_list[1]["elev_diff"]
    else:
        diff_data = dem_list[1]["elevation"] - dem_list[0]["elevation"]

    if ax is None:
        _, ax = plt.subplots(figsize=(14, 6))
    box_data = []  # difference data per group
    box_labels = []  # xtick labels for box groups
    counts = []  # group observation counts
    # loop over terrain groups and extract differences by group
    for i in range(len(terrain_groups) - 1):
        b1, b2 = terrain_groups[i], terrain_groups[i + 1]

        if isinstance(dem_list[1], gpd.GeoDataFrame):
            mask = (dem_list[1][terrain_v] > b1) & (dem_list[1][terrain_v] <= b2)
            data = dem_list[1].loc[mask, "elev_diff"].dropna()
        else:
            mask = (dem_list[0][terrain_v] > b1) & (dem_list[0][terrain_v] <= b2)
            data = diff_data.where(mask).to_numpy()
            data = data[~np.isnan(data)]

        # minimum 30 observations per group
        if len(data) >= 30:
            box_data.append(data)
            box_labels.append(f"{b1}-{b2}")
            counts.append(len(data))

    if len(box_data) == 0:
        msg_counts = "No groups satisfied the minimum 30 observation count threshold."
        raise ValueError(msg_counts)

    ax.boxplot(
        box_data,
        orientation="vertical",
        patch_artist=True,
        showfliers=True,
        boxprops={"facecolor": "lightgray", "color": "black"},
        flierprops={"marker": "x", "color": "black", "markersize": 5, "alpha": 0.6},
        medianprops={"color": "black"},
        positions=np.arange(len(box_data)),
    )

    # second axis for observation counts
    ax2 = ax.twinx()
    ax2.plot(np.arange(len(counts)), counts, "o", color="orange", alpha=0.6)
    # dynamic y-ticks for observation counts
    magnitude = 10 ** np.floor(np.log10(max(counts)))
    min_count = np.floor(min(counts) / magnitude) * magnitude
    max_count = np.ceil(max(counts) / magnitude) * magnitude
    ticks = np.linspace(min_count, max_count, 11)
    ax2.set_yticks(ticks)
    ax2.spines["right"].set_color("orange")
    ax2.tick_params(axis="y", colors="orange")
    ax2.set_ylabel("Count", color="orange", fontsize=12)

    # original axis elements
    ax.axhline(0, color="black", linestyle="dashed", linewidth=0.7)
    global_median = np.nanmedian(diff_data.to_numpy())
    ax.axhline(
        global_median,
        color="magenta",
        linestyle="dashed",
        linewidth=0.7,
        label=f"Global Med: {global_median:.2f}m",
        alpha=0.8,
    )
    # Dynamic ylims for the elevation differences
    # Zooms into whiskers extent
    if ylims is None:
        # get whiskers to fit the ylims of the plot
        data_concat = np.concatenate(box_data)
        q1, q3 = np.percentile(data_concat, [25, 75])
        iqr = q3 - q1
        ymin = np.floor(q1 - 1.5 * iqr)
        ymax = np.ceil(q3 + 1.5 * iqr)
        # force include 0 in y axis
        ymin = min(ymin, 0)
        ymax = max(ymax, 0)
        n_ticks = 11
        spacing = max(abs(ymin), abs(ymax)) * 2 / (n_ticks - 1)
        spacing = np.ceil(spacing)
        ymin = np.floor(ymin / spacing) * spacing
        ymax = np.ceil(ymax / spacing) * spacing
        yticks = np.arange(ymin, ymax + spacing, spacing)
        ylims = [ymin, ymax]
        ax.set_yticks(yticks)
    ax.set_ylim(ylims)
    ax.set_xticks(np.arange(len(box_labels)))
    ax.set_xticklabels(box_labels, rotation=45, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel(terrain_v, fontsize=12)
    ax.set_title(title, fontsize=14)

    lines1, labels1 = ax.get_legend_handles_labels()
    ax.legend(lines1, labels1, loc="best", fontsize=10)

    return ax


# slope wrapper for boxplot_terrain_diff()
def boxplot_slope(
    dem_list: list[xr.Dataset | gpd.GeoDataFrame],
    ax: plt.Axes | None = None,
    slope_bins: np.ndarray | None = None,
    ylims: list[float] | None = None,
    title: str = "Elevation Differences (m) by Slope",
    ylabel: str = "Elevation Difference (m)",
    elev_col: str | None = None,
) -> plt.Axes:
    """
    Boxplots of elevation differences grouped by elevation values.
    Grouped on the first dem in dem_list, which must have an 'slope' variable

    This is a wrapper around boxplot_terrain_diff() specifically for slope analysis
    with groups from 0-45 degrees in 1 degree increments.

    Parameters
    ----------
    dem_list : list[xr.Dataset | gpd.GeoDataFrame]
        List containing exactly two elevation datasets to compare
    slope_bins : np.ndarray
        Array defining the edges of terrain groups. e.g. if set to np.arange(0, 46, 1),
        the groups will be [0, 1, 2, ..., 45]
        Default is np.arange(0, 46, 1)
    ylims : list[float] | None
        The y-axis limits for the plot
    title : str
        The title of the plot
    ylabel : str
        The ylabel of the plot
    elev_col : str | None
        Column name containing elevation values if using GeoDataFrame
    show : bool
        Whether to display the plot

    Returns
    -------
    plt.Axes
        The matplotlib axes object containing the plot
    """
    if slope_bins is None:
        slope_bins = np.arange(0, 46, 1)
    return boxplot_terrain_diff(
        dem_list=dem_list,
        ax=ax,
        terrain_v="slope",
        terrain_groups=np.arange(0, 46, 1),
        ylims=ylims,
        title=title,
        ylabel=ylabel,
        elev_col=elev_col,
    )


def boxplot_elevation(
    dem_list: list[xr.Dataset | gpd.GeoDataFrame],
    ax: plt.Axes | None = None,
    elevation_bins: np.ndarray | None = None,
    ylims: list[float] | None = None,
    title: str = "Elevation Differences (m) by Source Elevation",
    ylabel: str = "Elevation Difference (m)",
    elev_col: str | None = None,
) -> plt.Axes:
    """
    Boxplots of elevation differences grouped by elevation values.
    Grouped on the first dem in dem_list, which must have an 'elevation' variable

    This is a wrapper around boxplot_terrain_diff() specifically for elevation analysis
    with groups defined by elevation_bins. If elevation_bins is None, bins are created
    from the min/max elevations of the first DEM in steps of 300m.

    Parameters
    ----------
    dem_list : list[xr.Dataset | gpd.GeoDataFrame]
        List containing exactly two elevation datasets to compare
    elevation_bins : np.ndarray | None
        Array defining the edges of elevation groups
        If None, bins are created from min/max of source DEM in 300m steps
    ylims : list[float] | None
        The y-axis limits for the plot
    title : str
        The title of the plot
    ylabel : str
        The ylabel of the plot
    elev_col : str | None
        Column name containing elevation values if using GeoDataFrame

    Returns
    -------
    plt.Axes
        The matplotlib axes object containing the plot
    """
    if elevation_bins is None:
        elev_min = np.floor(np.nanmin(dem_list[0].elevation) / 100) * 100
        elev_max = np.ceil(np.nanmax(dem_list[0].elevation) / 100) * 100
        elevation_bins = np.arange(elev_min, elev_max + 300, 300)

    return boxplot_terrain_diff(
        dem_list=dem_list,
        ax=ax,
        terrain_v="elevation",
        terrain_groups=elevation_bins,
        ylims=ylims,
        title=title,
        ylabel=ylabel,
        elev_col=elev_col,
    )


def boxplot_aspect(
    dem_list: list[xr.Dataset | gpd.GeoDataFrame],
    ax: plt.Axes | None = None,
    aspect_bins: np.ndarray | None = None,
    ylims: list[float] | None = None,
    title: str = "Elevation Differences (m) by Source Aspect",
    ylabel: str = "Elevation Difference (m)",
    elev_col: str | None = None,
) -> plt.Axes:
    """
    Boxplots of elevation differences grouped by aspect values.
    Grouped on the first dem in dem_list, which must have an 'aspect' variable

    This is a wrapper around boxplot_terrain_diff() specifically for aspect analysis
    with groups defined by aspect_bins. If aspect_bins is None, bins are created
    from 0-360 degrees in steps of 10 degrees.

    Parameters
    ----------
    dem_list : list[xr.Dataset | gpd.GeoDataFrame]
        List containing exactly two elevation datasets to compare
    aspect_bins : np.ndarray | None
        Array defining the edges of aspect groups
        If None, bins are created from 0-360 in steps of 10 degrees
    ylims : list[float] | None
        The y-axis limits for the plot
    title : str
        The title of the plot
    ylabel : str
        The ylabel of the plot
    elev_col : str | None
        Column name containing elevation values if using GeoDataFrame

    Returns
    -------
    plt.Axes
        The matplotlib axes object containing the plot
    """
    if aspect_bins is None:
        aspect_bins = np.arange(0, 370, 10)

    return boxplot_terrain_diff(
        dem_list=dem_list,
        ax=ax,
        terrain_v="aspect",
        terrain_groups=aspect_bins,
        ylims=ylims,
        title=title,
        ylabel=ylabel,
        elev_col=elev_col,
    )


def _create_stats_legend(
    ax: plt.Axes,
    stats_dict: dict[str, float],
) -> None:
    """
    Helper function for hist_esa. Create a unified legend with statistics for elevation difference plots.

    Parameters
    ----------
    ax : plt.Axes
        The matplotlib axes object to add the legend to
    stats_dict : dict[str, float]
        Dictionary containing statistics to display

    Returns
    -------
    None
        Modifies the input axes object directly
    """
    # Assume mean and median already plotted?
    handles, _ = ax.get_legend_handles_labels()

    # Append Std, NMAD, and Count to legend
    handles.append(
        plt.Line2D([0], [0], color="none", label=f"{'Std':<7}{stats_dict['Std']:>5.2f}")
    )
    handles.append(
        plt.Line2D(
            [0], [0], color="none", label=f"{'NMAD':<7}{stats_dict['NMAD']:>5.2f}"
        )
    )
    handles.append(
        plt.Line2D([0], [0], color="none", label=f"{'Min':<7}{stats_dict['Min']:>5.1f}")
    )
    handles.append(
        plt.Line2D([0], [0], color="none", label=f"{'Max':<7}{stats_dict['Max']:>5.1f}")
    )
    handles.append(
        plt.Line2D(
            [0],
            [0],
            color="none",
            label=f"{'Count':<7}{stats_dict['Count']:>5.1e}".replace("e+0", "e"),
        )
    )

    # NOTE: monospace font is key for legend alignment
    ax.legend(
        handles=handles, loc="upper right", fontsize=8, prop={"family": "monospace"}
    )


def hist_esa(
    dem_list: list[xr.Dataset | gpd.GeoDataFrame],
    elev_col: str | None = None,
    min_count: int = 30,
) -> plt.Figure:
    """
    Histogram of elevation differences between DEMs or point data, grouped by ESA World Cover 2021 land cover class.

    Parameters
    ----------
    dem_list : list[xr.Dataset | gpd.GeoDataFrame]
        List containing two elevation datasets to compare
    elev_col : str | None
        Column name containing elevation values if using point data
    min_count : int
        Minimum number of points required for a land cover class to be included

    Returns
    -------
    np.array
        An array containing the matplotlib axes objects
    """
    if len(dem_list) != 2:
        msg_len = "dem_list must contain exactly two datasets"
        raise ValueError(msg_len)

    # calculate elevation differences if second dataset is a GDF
    if isinstance(dem_list[1], gpd.GeoDataFrame):
        if elev_col is None:
            col_msg = "elev_col must be provided when using point data"
            raise ValueError(col_msg)

        da_points = dem_list[1].get_coordinates().to_xarray()
        samples = (
            dem_list[0]
            .interp(da_points)
            .drop_vars(["band", "spatial_ref"])
            .to_dataframe()
        )
        dem_list[1]["elev_diff"] = (
            dem_list[1][elev_col].to_numpy() - samples["elevation"].to_numpy()
        )
        diff_data = dem_list[1]
    else:  # Both are rasters
        diff_data = dem_list[1]["elevation"] - dem_list[0]["elevation"]

    classmap = WorldCover().classmap

    # Get WorldCover dataset
    bounds = dem_list[0].rio.bounds()
    gf_search = gpd.GeoDataFrame(
        geometry=[box(*bounds)], crs=dem_list[0].rio.crs
    ).to_crs(epsg=4326)
    gf_wc = search(dataset="worldcover", intersects=gf_search, datetime=["2021"])
    # For tiny datasets (e.g. tests), mask=True can fail..
    # TODO: revisit masking here
    ds_wc = to_dataset(gf_wc, bands=["map"], aoi=gf_search)

    ds_wc = ds_wc.rio.reproject(dem_list[0].rio.crs).rio.reproject_match(dem_list[0])

    # plotting logic for raster data (xarray)
    if isinstance(diff_data, xr.DataArray):
        bin_width = 0.25
        bins = np.arange(-30, 30 + bin_width, bin_width)
        filtered_classes = []

        for class_value, class_info in classmap.items():
            class_mask = ds_wc.isel(time=0)["map"].to_numpy() == class_value
            class_mask_da = xr.DataArray(
                class_mask, coords=[diff_data.y, diff_data.x], dims=["y", "x"]
            )
            class_data = diff_data.where(class_mask_da, drop=True)
            class_data_flat = class_data.to_numpy().flatten()
            class_data_flat = class_data_flat[~np.isnan(class_data_flat)]

            if len(class_data_flat) >= min_count:
                filtered_classes.append((class_value, class_info, class_data_flat))

        num_classes = len(filtered_classes)
        ncols, nrows = 2, (num_classes + 1) // 2
        _, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(12, 4 * nrows),
            constrained_layout=True,
            sharex=True,
            sharey=True,
        )

        for idx, (ax, (_class_value, class_info, class_data_flat)) in enumerate(
            zip(axes.ravel(), filtered_classes, strict=False)
        ):
            ax.hist(
                class_data_flat, bins=bins, alpha=0.5, color=class_info["hex"], log=True
            )
            ax.set_title(class_info["description"], fontsize=12)
            ax.axvline(0, color="k", linestyle="dashed", linewidth=0.5)
            ax.set_xlim(-20, 20)
            # only set ylabel if axis is in first column
            if idx % ncols == 0:
                ax.set_ylabel("Log(Count)")
            else:
                ax.set_ylabel("")

            # calculate statistics
            stats_dict = utils.calc_stats(pd.Series(class_data_flat))
            median = stats_dict["Median"]
            mean = stats_dict["Mean"]
            ax.axvline(mean, label=f"{'Mean':<7}{mean:>5.2f}", color="magenta", lw=0.5)
            ax.axvline(
                median, label=f"{'Median':<7}{median:>5.2f}", color="cyan", lw=0.5
            )
            _create_stats_legend(ax, stats_dict)

        plt.suptitle(
            "Elevation Differences by Land Cover (ESA World Cover 2021)", fontsize=14
        )

    # plotting logic if second dataset is a GDF
    elif isinstance(diff_data, gpd.GeoDataFrame):
        # sample worldcover values at points
        land_cover = ds_wc.isel(time=0)["map"].interp(
            x=("points", diff_data.geometry.x.to_numpy()),
            y=("points", diff_data.geometry.y.to_numpy()),
            method="nearest",
        )
        diff_data["land_cover_class"] = land_cover.to_numpy()

        groups = diff_data.groupby("land_cover_class")
        # filter out groups with less than min_count points
        filtered_groups = [
            (label, group) for label, group in groups if len(group) >= min_count
        ]

        num_groups = len(filtered_groups)
        ncols, nrows = 2, (num_groups + 1) // 2
        _, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(12, 4 * nrows),
            constrained_layout=True,
            sharex=True,
            sharey=True,
        )

        for idx, (ax, (label, group)) in enumerate(
            zip(axes.ravel(), filtered_groups, strict=False)
        ):
            color = classmap.get(label, {}).get("hex", "#cccccc")
            ax.hist(
                group["elev_diff"].dropna(), bins=256, alpha=0.5, color=color, log=True
            )
            ax.axvline(0, color="k", linestyle="dashed", linewidth=1)
            ax.set_title(classmap.get(label, {}).get("description", ""), fontsize=12)
            ax.set_xlim(-20, 20)
            # only set ylabel if axis is in first column
            if idx % ncols == 0:
                ax.set_ylabel("Log(Count)")
            else:
                ax.set_ylabel("")

            group_no_nan = group["elev_diff"].dropna()

            stats_dict = utils.calc_stats(group_no_nan)
            median = stats_dict["Median"]
            mean = stats_dict["Mean"]
            ax.axvline(mean, label=f"{'Mean':<7}{mean:>5.2f}", color="magenta", lw=0.5)
            ax.axvline(
                median, label=f"{'Median':<7}{median:>5.2f}", color="cyan", lw=0.5
            )
            _create_stats_legend(ax, stats_dict)

        plt.suptitle(
            "Elevation Differences by Land Cover (ESA World Cover 2021)", fontsize=14
        )

    return axes
