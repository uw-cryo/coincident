"""
Selection of visualizations using matplotlib.
"""

from __future__ import annotations

import warnings
from typing import Any

import geopandas as gpd
import matplotlib.figure
import numpy as np
import pandas as pd
import pystac
import xarray as xr

# from matplotlib_scalebar.scalebar import ScaleBar
from scipy import stats
from shapely.geometry import box

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
from coincident.io.xarray import open_maxar_browse, to_dataset
from coincident.plot.utils import get_elev_diff
from coincident.search import search


@depends_on_optional("matplotlib")
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
            ax.set_aspect(aspect=1 / np.cos(np.deg2rad(mid_lat)))
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


@depends_on_optional("matplotlib")
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
    mid_lat = da.y[int(da.y.size / 2)].to_numpy()

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 11))

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

    ax.set_aspect(aspect=1 / np.cos(np.deg2rad(mid_lat)))
    ax.set_title(item.id)
    return ax


@depends_on_optional("matplotlib")
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


@depends_on_optional("matplotlib")
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


@depends_on_optional("matplotlib")
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
    cmin: float | None
        Minimum value for the elevation colorbar
    cmax: float | None
        Maximum value for the elevation colorbar
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
            ax=ax, cmap="gray", alpha=1.0, add_colorbar=False, add_labels=False
        )
        if "alpha" not in kwargs:
            kwargs["alpha"] = 0.5
    gf.plot(ax=ax, column=column, cmap=cmap, **kwargs)
    ax.set_title(title)
    return ax


def plot_diff_hist(
    da: xr.DataArray | pd.Series,
    ax: plt.Axes | None = None,
    dmin: float | None = -15,
    dmax: float | None = 15,
    **kwargs: Any,
) -> plt.Axes:
    """
    Histogram of elevation differences with median and mean lines.

    Parameters
    ----------
    da: xarray.DataArray | pd.Series
        DataArray or series containing elevation differences to plot
    ax: plt.Axes
        Axis on which to plot the histogram
    dmin: float
        Minimum value for histogram x-axis range, default -15
    dmax: float
        Maximum value for histogram x-axis range, default 15
    **kwargs: Any
        Artist properties passed to ax.hist()

    Returns
    -------
    plt.Axes
        The Axes object with the histogram plot
    """
    if ax is None:
        _, ax = plt.subplots()

    if isinstance(da, xr.DataArray):
        _ = da.plot.hist(ax=ax, bins=100, range=(dmin, dmax), color="gray", **kwargs)
    else:
        da.hist(ax=ax, bins=100, range=(dmin, dmax), color="gray", **kwargs)
    med = np.nanmedian(da)
    mean = np.nanmean(da)
    ax.axvline(0, color="k", lw=1)
    # NOTE: hardcoded unit label of meters
    ax.axvline(med, label=f"med: {med:.2f}m", color="cyan", lw=0.5)
    ax.axvline(mean, label=f"mean: {mean:.2f}m", color="magenta", lw=0.5)
    ax.set_xlabel("Elevation Difference (m)")
    ax.set_ylabel("Frequency")
    ax.legend(loc="best", fontsize="small")
    ax.grid(False)
    return ax


# TODO: Add add_hillshade=False?
@depends_on_optional("matplotlib")
def compare_dems(
    dem_list: list[xr.Dataset],
    gdf_dict: dict[str, tuple[gpd.GeoDataFrame, str]] | None = None,
    ds_wc: xr.Dataset | None = None,
    elevation_cmap: str = "inferno",
    add_hillshade: bool = True,
    diff_cmap: str = "RdBu",
    elevation_clim: tuple[float, float] | None = None,
    diff_clim: tuple[float, float] = (-10, 10),
    figsize: tuple[float, float] = (17, 11),
) -> plt.Axes:
    """
    Create a panel figure comparing DEM and altimetry elevations.

    Where the first row shows elevation maps over hillshades (if add_hillshade is True) altimetry points will be plotted over the hillshade of the first DEM (if add_hillshade is True)
    The second row shows elevation differences, with the first plot being the Landcover of the scene First DEM provided will be the 'source' DEM where all other elevations will be compared to
    The third row shows the histograms of these differences

    Parameters
    ----------
    dem_list: list
        List of xr.Datasets. Each Dataset should have .elevation and .hillshade variables.
        The first DEM in this list will be the 'source' DEM and will served as a reference for differencing
        e.g.  [dem_1,dem_2,dem_3] will result in diff_1 = dem_2 - dem_1 and diff_2 = dem_3 - dem_1
    gdf_dict: Optional Dict[str, Tuple[gpd.GeoDataFrame, str]]
        Dictionary where keys are subplot names and values are (GeoDataFrame, column_name) pairs.
        column_name denotes elevation values to plot (e.g. h_li for ICESat-2 ATL06)
    ds_wc: Optional xr.Dataset
        WorldCover dataset to plot. If this is not provided, ESA WorldCover is used.
    elevation_cmap: Optional str
        Colormap for elevation. Default: 'inferno'
    add_hillshade: Optional bool
        Whether to plot hillshades under your DEMs. Default True.
        If True, your DEMs should have a .hillshade variable
    diff_cmap: Optional str
        Colormap for elevation differences. Default: 'RdBu'
    elevation_clim: Optional Tuple[float, float]
        Tuple for elevation color limits
    diff_clim: Optional Tuple[float, float]
        Tuple for difference color limits. Default (-10,10)
    figsize: Optional Tuple[int, int]
        Figure size tuple. Will be dynamically adjusted based on number of DEMs and GeoDataFrames
        so you typically don't want to adjust this

    Returns
    ----------
        Dictionary containing the axes objects for each subplot

    Notes
    -----
    All inputs assumed to have the same CRS and to be aligned
    All DEMs must have variable 'elevation'
    A minimum requirement of 2 DEMs in dem_list
    A maximum of 5 total datasets (dem_list + gdf_dict) to be passed
    """
    # TODO: better DEM titles in plots
    # Calculate number of columns for the plot (DEMs + GeoDataFrames if provided)
    n_dems = len(dem_list)
    gdf_keys = list(gdf_dict.keys()) if gdf_dict else []
    n_columns = max(2, n_dems) + len(gdf_keys)
    n_columns = min(n_columns, 5)

    # Dynamically adjust figure size
    figsize = (3 * n_columns, 11)
    mosaic = [
        [f"dem_{i}" for i in range(n_dems)] + gdf_keys,
        ["worldcover"]
        + [f"diff_{i}" for i in range(1, n_dems)]
        + [f"diff_{g}" for g in gdf_keys],
        ["wc_legend"]
        + [f"hist_{i}" for i in range(1, n_dems)]
        + [f"hist_{g}" for g in gdf_keys],
    ]
    gs_kw = {"width_ratios": [1] * len(mosaic[0]), "height_ratios": [2, 2, 1]}
    _, axd = plt.subplot_mosaic(
        mosaic, gridspec_kw=gs_kw, figsize=figsize, layout="constrained"
    )
    reference_ax = axd["dem_0"]
    for key in axd:
        if not key.startswith("hist") and not key.startswith("wc"):
            axd[key].sharex(reference_ax)
            axd[key].sharey(reference_ax)

    # TODO: move this to a top-level helper function, and set aspect based on mid scene latitude
    # or equal if UTM, etc. ?
    def _style_ax(
        ax: plt.Axes,
        facecolor: str | None = None,
        aspect: str | None = None,
    ) -> None:
        """Helper function to style axes (remove ticks and set face color)."""
        ax.set_xticks([])
        ax.set_yticks([])
        if facecolor:
            ax.set_facecolor(facecolor)
        if aspect:
            ax.set_aspect(aspect)

    # Generate WorldCover dataset if not provided
    if ds_wc is None:
        bounds = dem_list[0].rio.bounds()
        gf_search = gpd.GeoDataFrame(
            geometry=[box(*bounds)],
            crs=dem_list[0].rio.crs,
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

    # TOP ROW: Elevation maps
    # ---------------------------
    # raster elevs
    for i, dem in enumerate(dem_list):
        dem_key = f"dem_{i}"
        plot_dem(
            dem.elevation,
            axd[dem_key],
            elevation_cmap,
            da_hillshade=dem.hillshade if add_hillshade else None,
            add_colorbar=False,
            vmin=(elevation_clim[0] if elevation_clim else None),
            vmax=(elevation_clim[1] if elevation_clim else None),
            title=dem_key,
        )
        _clear_labels(axd[dem_key])
        _style_ax(axd[dem_key])
        if i == n_columns - 1:  # Add colorbar to the last GeoDataFrame plot
            # images[1] since [0] would be the hillshade
            plt.colorbar(axd[dem_key].images[1], ax=axd[dem_key], label="Elevation (m)")
        # TODO: add scalebar later, allow for lon/lat or meters
        # elif i == 0:
        #     dx = dem.rio.resolution()[0]
        #     scalebar = ScaleBar(
        #         dx,  # meters per pixel from resolution
        #         "m",
        #         length_fraction=0.25,
        #         location="lower right",
        #         box_alpha=0.5,
        #         box_color="white",
        #         color="black",
        #         pad=0.5,
        #     )
        #     axd[dem_key].add_artist(scalebar)

    # point elevs
    # need this for MyPy to handle the chance the gdf_dict is None
    # despite the logic in the loop "if gdf_dict else []"
    if gdf_dict:
        # Plot GeoDataFrames (e.g., IS2 and GEDI)
        for key, (gdf, column) in gdf_dict.items() if gdf_dict else []:
            plot_altimeter_points(
                gdf,
                column,
                ax=axd[key],
                title=key,
                da_hillshade=dem_list[0].hillshade if add_hillshade else None,
                vmin=(elevation_clim[0] if elevation_clim else None),
                vmax=(elevation_clim[1] if elevation_clim else None),
                s=(figsize[0] * figsize[1])
                / n_columns**2,  # Dynamic point size based on figure dimensions
            )
            _style_ax(axd[key])
            if key == list(gdf_dict.keys())[-1]:  # Add colorbar to the last gdf plot
                plt.colorbar(
                    axd[key].collections[0], ax=axd[key], label="Elevation (m)"
                )

    # MIDDLE ROW: Differences
    # ---------------------------
    plot_esa_worldcover(
        ds_wc.isel(time=0),
        ax=axd["worldcover"],
        cax=axd["wc_legend"],
        add_colorbar=True,
        add_labels=False,
    )
    _style_ax(axd["worldcover"])
    axd["worldcover"].set_title("ESA WorldCover 2021")

    # raster diffs
    for i, dem in enumerate(dem_list[1:], start=1):
        diff = get_elev_diff(dem, dem_list[0])
        diff_key = f"diff_{i}"
        if add_hillshade:
            dem_list[0].hillshade.plot.imshow(
                ax=axd[diff_key],
                cmap="gray",
                alpha=1.0,
                add_colorbar=False,
                add_labels=False,
            )
        diff.elev_diff.squeeze().plot.imshow(
            ax=axd[diff_key],
            cmap=diff_cmap,
            add_colorbar=False,
            add_labels=False,
            vmin=diff_clim[0],
            vmax=diff_clim[1],
            alpha=0.5 if add_hillshade else 1.0,
        )
        _style_ax(axd[diff_key], facecolor="black")
        axd[diff_key].set_title(f"dem_{i} minus dem_0")
        if i == n_columns - 1:  # Add colorbar to the last gdf plot
            plt.colorbar(
                axd[diff_key].images[1] if add_hillshade else axd[diff_key].images[0],
                ax=axd[diff_key],
                label="Elevation Difference (m)",
            )
    # gdf diffs
    # need this for MyPy to handle the chance the gdf_dict is None
    # despite the logic in the loop "if gdf_dict else []"
    if gdf_dict:
        for _, (key, (gdf, column)) in enumerate(gdf_dict.items() if gdf_dict else []):
            diff_key = f"diff_{key}"
            if add_hillshade:
                dem_list[0].hillshade.plot.imshow(
                    ax=axd[diff_key],
                    cmap="gray",
                    alpha=1.0,
                    add_colorbar=False,
                    add_labels=False,
                )
            gf_diff = get_elev_diff(gdf, dem_list[0], source_col=column)
            plot_altimeter_points(
                gf_diff,
                "elev_diff",
                axd[diff_key],
                cmap=diff_cmap,
                vmin=diff_clim[0],
                vmax=diff_clim[1],
                title=f"{column} minus dem_0",
                s=(figsize[0] * figsize[1])
                / n_columns**2,  # Dynamic sizing based on figure dimensions
            )
            _style_ax(axd[diff_key], facecolor="black")
            if key == list(gdf_dict.keys())[-1]:  # Add colorbar to the last gdf plot
                plt.colorbar(
                    axd[diff_key].collections[0],
                    ax=axd[diff_key],
                    label="Elevation Difference (m)",
                )

    # BOTTOM ROW: Histograms
    # ---------------------------
    axd["wc_legend"].axis("off")

    # raster hists
    for i, dem in enumerate(dem_list[1:], start=1):
        hist_key = f"hist_{i}"
        plot_diff_hist(
            get_elev_diff(dem, dem_list[0]).elev_diff,
            ax=axd[hist_key],
            dmin=diff_clim[0],
            dmax=diff_clim[1],
        )
        axd[hist_key].set_title("")

    # gdf hists
    for _, (key, (gdf, column)) in enumerate(gdf_dict.items() if gdf_dict else []):
        hist_key = f"hist_{key}"
        diff = get_elev_diff(gdf, dem_list[0], source_col=column)
        plot_diff_hist(
            diff.elev_diff, ax=axd[hist_key], dmin=diff_clim[0], dmax=diff_clim[1]
        )
        axd[hist_key].set_title("")

    return axd


@depends_on_optional("matplotlib")
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
@depends_on_optional("matplotlib")
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


@depends_on_optional("matplotlib")
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


@depends_on_optional("matplotlib")
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


@depends_on_optional("matplotlib")
def _create_stats_legend(
    ax: plt.Axes, stats_dict: dict[str, float], color_dict: dict[str, str]
) -> None:
    """
    Helper function for hist_esa. Create a unified legend with statistics for elevation difference plots.

    Parameters
    ----------
    ax : plt.Axes
        The matplotlib axes object to add the legend to
    stats_dict : dict[str, float]
        Dictionary containing statistics to display
    color_dict : dict[str, str]
        Dictionary mapping statistic names to colors

    Returns
    -------
    None
        Modifies the input axes object directly
    """
    legend_elements = []

    for label in ["Mean", "Med"]:
        color = color_dict.get(label, "black")
        legend_elements.append(
            plt.Line2D([0], [0], color=color, label=f"{label}: {stats_dict[label]:.2f}")
        )

    for label in ["Std", "NMAD"]:
        legend_elements.append(
            plt.Line2D(
                [0], [0], color="none", label=f"{label}: {stats_dict[label]:.2f}"
            )
        )

    legend_elements.append(
        plt.Line2D([0], [0], color="none", label=f"Count: {int(stats_dict['Count'])}")
    )

    legend = ax.legend(
        handles=legend_elements,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        fontsize=10,
        facecolor="lightgray",
        edgecolor="none",
    )
    legend.get_frame().set_alpha(0.5)


@depends_on_optional("matplotlib")
def hist_esa(
    dem_list: list[xr.Dataset | gpd.GeoDataFrame],
    elev_col: str | None = None,
    min_count: int = 30,
    color_dict: dict[str, str] | None = None,
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
    color_dict : dict[str, str]
        Dictionary mapping statistics to colors for plotting

    Returns
    -------
    np.array
        An array containing the matplotlib axes objects
    """
    if len(dem_list) != 2:
        msg_len = "dem_list must contain exactly two datasets"
        raise ValueError(msg_len)
    if color_dict is None:
        color_dict = {
            "Mean": "magenta",
            "Med": "deepskyblue",
            "Std": "black",
            "NMAD": "black",
            "Count": "black",
        }
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
            ax.axvline(0, color="k", linestyle="dashed", linewidth=1)
            ax.set_xlim(-20, 20)
            # only set ylabel if axis is in first column
            if idx % ncols == 0:
                ax.set_ylabel("Log(Count)")
            else:
                ax.set_ylabel("")

            # calculate statistics
            stats_dict = {
                "Mean": np.mean(class_data_flat),
                "Med": np.median(class_data_flat),
                "Std": np.std(class_data_flat),
                "NMAD": stats.median_abs_deviation(class_data_flat, scale="normal"),
                "Count": len(class_data_flat),
            }

            _create_stats_legend(ax, stats_dict, color_dict)

            # vertical lines for med and mean
            ax.axvline(stats_dict["Med"], color="deepskyblue", lw=0.5)
            ax.axvline(stats_dict["Mean"], color="magenta", lw=0.5)

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
            stats_dict = {
                "Mean": np.mean(group_no_nan),
                "Med": np.median(group_no_nan),
                "Std": np.std(group_no_nan),
                "NMAD": stats.median_abs_deviation(group_no_nan, scale="normal"),
                "Count": len(group_no_nan),
            }

            _create_stats_legend(ax, stats_dict, color_dict)

            # vertical lines for med and mean
            ax.axvline(stats_dict["Med"], color="deepskyblue", lw=0.5)
            ax.axvline(stats_dict["Mean"], color="magenta", lw=0.5)

        plt.suptitle(
            "Elevation Differences by Land Cover (ESA World Cover 2021)", fontsize=14
        )

    return axes
