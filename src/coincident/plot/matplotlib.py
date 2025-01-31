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
from coincident.search import search


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

    # set aspect ratio according to mid scene latitude
    if ds.rio.crs.to_epsg() == 4326:
        mid_lat = da.latitude[int(da.latitude.size / 2)].to_numpy()  # PD011
        ax.set_aspect(aspect=1 / np.cos(np.deg2rad(mid_lat)))

    colorbar = fig.colorbar(
        cm.ScalarMappable(norm=normalizer, cmap=cmap),
        boundaries=boundaries,
        values=values,
        cax=fig.axes[1].axes,  # hard to get right w/ subplot_mosaic panel...
    )
    colorbar.set_ticks(ticks, labels=tick_labels)
    return ax


# NOTE: below is a temporary function
# this function is only used in the compare_dem function to make the
# horizontal colorbar and such easier
# hoping to merge this with plot_esa_worldcover to make a more-dynamic function
@depends_on_optional("matplotlib")
def plot_worldcover_custom_ax(
    ds: xr.Dataset,
    ax: plt.Axes,
    cax: plt.Axes,
    add_colorbar: bool = True,
    add_labels: bool = True,
) -> plt.Axes:
    """
    Plots ESA WorldCover data on a custom Matplotlib Axes.
    THIS FUNCTION IS PURELY TO MAKE compare_dems() EASIER
    DUE TO DIFFICULTIES OF MAKING plot_esa_worldcover() INCLUDE
    AN OPTION TO PLOT A HORIZONTAL CBAR BELOW THE AXIS OF INTEREST

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

    References
    ----------
    https://planetarycomputer.microsoft.com/dataset/esa-worldcover#Example-Notebook
    """
    # Map class values to colors and descriptions
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

    normalizer = matplotlib.colors.Normalize(vmin=0, vmax=255)

    da = ds.to_dataarray().squeeze()
    da.plot(
        ax=ax, cmap=cmap, norm=normalizer, add_labels=add_labels, add_colorbar=False
    )
    ax.set_title("")

    # set aspect ratio according to mid scene latitude
    if ds.rio.crs.to_epsg() == 4326:
        mid_lat = da.latitude[int(da.latitude.size / 2)].to_numpy()  # PD011
        ax.set_aspect(aspect=1 / np.cos(np.deg2rad(mid_lat)))

    if add_colorbar:
        # Override standard xarray colorbar with custom one
        fig = ax.get_figure()
        # if cax==None:
        #    cax=ax.inset_axes([0, -0.35, 1, 0.1])
        colorbar = fig.colorbar(
            matplotlib.cm.ScalarMappable(norm=normalizer, cmap=cmap),
            boundaries=boundaries,
            values=values,
            # cax=fig.axes[1].axes, # hard to get right w/ subplot_mosaic panel...
            ax=cax,
            # cax=cax,
            orientation="horizontal",
            # location='top'
            pad=-1.3,  # hack to get colorbar in subplot mosaic in the right place
        )
        colorbar.set_ticks(ticks, labels=tick_labels, rotation=90)

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


# NOTE: xdem.DEM.hillshade method expects the same x and y resolution
@depends_on_optional("matplotlib")
def hillshade(
    da: xr.DataArray, azi: int = 45, alt: int = 45
) -> tuple[tuple[str, ...], np.ndarray]:
    """
    Create hillshade array from DEM.
    Pass with dem['hillshade'] = coincident.plot.create_hillshade(dem.elevation)

    Parameters
    ----------
    da: xarray.DataArray
        Dataarray containing elevation values
    azi: int
        The shading azimuth in degrees (0-360°) going clockwise, starting from north.
    alt: int
        The shading altitude in degrees (0-90°). 90° is straight from above.

    Returns
    -------
    tuple
        (dims, data) tuple where dims are dimension names and data is uint8 array
    """
    # Ensure 2D array by squeezing extra dimensions
    da = da.squeeze()

    x, y = np.gradient(da)
    slope = np.pi / 2.0 - np.arctan(np.sqrt(x * x + y * y))
    aspect = np.arctan2(-x, y)
    azimuth = 360.0 - azi
    azimuthrad = azimuth * np.pi / 180.0
    altituderad = alt * np.pi / 180.0

    shaded = np.sin(altituderad) * np.sin(slope) + np.cos(altituderad) * np.cos(
        slope
    ) * np.cos(azimuthrad - aspect)

    data = 255 * (shaded + 1) / 2
    data[data < 0] = 0  # set hillshade values to min of 0.

    return da.dims, data.astype(np.uint8)


@depends_on_optional("matplotlib")
def clear_labels(ax: matplotlib.axes.Axes) -> None:
    """
    Removes x and y axis labels from a matplotlib plot.

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
    ax: matplotlib.axes.Axes,
    cmap: str = "inferno",
    title: str = "",
    da_hillshade: xr.DataArray | None = None,
    **kwargs: Any,
) -> plt.Axes:
    """
    Plots a DEM with an option to plot a hillshade underneath
    NOTE: set kwarg alpha=0.5 if you pass in a hillshade.
        hillshade alpha is set to 1.0, and nee to thus set overlaying DEM alpha to 0.5

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
    """
    # hillshade with alpha=1.0
    if isinstance(da_hillshade, xr.DataArray):
        da_hillshade.plot.imshow(
            ax=ax, cmap="gray", alpha=1.0, add_colorbar=False, add_labels=False
        )
        ax.set_aspect("equal")
    # plot the DEM with alpha=0.5
    da.plot.imshow(ax=ax, cmap=cmap, **kwargs)
    ax.set_title(title)
    ax.set_aspect("equal")
    clear_labels(ax)
    return ax


# TODO: fix kwargs implementation
# I was having issues controlling the alpha and vmin/vmax with kwargs
# so I made them their own separate args
@depends_on_optional("matplotlib")
def plot_altimeter_points(
    gf: gpd.GeoDataFrame,
    ax: matplotlib.axes.Axes,
    column: str,
    cmap: str = "inferno",
    title: str = "",
    da_hillshade: xr.DataArray | None = None,
    **kwargs: Any,
) -> plt.Axes:
    """
    Plots laser altimeter point data with an optional hillshade background.
    NOTE: set kwarg alpha=0.5 if you pass in a hillshade.
        hillshade alpha is set to 1.0, and nee to thus set overlaying DEM alpha to 0.5

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
    """
    # Plot the hillshade if available
    if isinstance(da_hillshade, xr.DataArray):
        da_hillshade.plot.imshow(
            ax=ax, cmap="gray", alpha=1.0, add_colorbar=False, add_labels=False
        )
    # Plot the altimeter data with overlay
    gf.plot(ax=ax, column=column, cmap=cmap, **kwargs)
    ax.set_title(title)
    ax.set_aspect("equal")
    clear_labels(ax)
    return ax


def get_elev_diff(
    source: xr.Dataset | gpd.GeoDataFrame,
    reference: xr.DataArray,
    source_col: str | None = None,
    diff_col_name: str | None = "elev_diff",
) -> xr.Dataset | gpd.GeoDataFrame:
    """
    Calculate elevation differences between source and reference datasets.
    This assumes that the source and reference dataset are in the same CRS.
    If the source is an xr.Dataset, this assumes the source dataset will
    have elevation values stored in variable name 'elevation'.

    Parameters
    ----------
    source: xr.Dataset | gpd.GeoDataFrame
        Source elevation data, either as raster or points
    reference: xr.DataArray
        Reference elevation data as raster
    source_col: str | None
        Column name containing elevation values to difference if source is GeoDataFrame (e.g. h_li for ICESat-2 ATL06)
    diff_col_name: str
        Resulting variable name or column name for elevation differences depending if output
        is xarray.Dataset or GeoDataFrame. Default to 'elev_diff'
    Returns
    -------
    xr.Dataset | gpd.GeoDataFrame
        Source dataset with added elevation differences
    """
    if isinstance(source, xr.Dataset):
        # raster source
        source_matched = source.rio.reproject_match(reference)
        source_elev_diff = source_matched.elevation - reference.elevation
        source_matched[diff_col_name] = source_elev_diff
        return source_matched

    if isinstance(source, gpd.GeoDataFrame):
        # point source
        if source_col is None:
            msg_val_error = "source_col must be specified when source is GeoDataFrame"
            raise ValueError(msg_val_error)
        # sample reference at point locations
        da_points = source.get_coordinates().to_xarray()
        samples = reference.interp(da_points).to_dataframe()
        source["elev_diff"] = (
            source[source_col].to_numpy() - samples["elevation"].to_numpy()
        )
        return source

    msg_type_error = "source must be xarray Dataset or GeoDataFrame"
    raise TypeError(msg_type_error)


def plot_diff_hist(
    da: xr.DataArray | pd.Series,
    ax: matplotlib.axes.Axes,
    dmin: float | None = -15,
    dmax: float | None = 15,
    **kwargs: Any,
) -> plt.Axes:
    """
    Plots a histogram of elevation differences with median and mean lines.
    x-axis is difference in elevation and y-axis is frequency.

    Parameters
    ----------
    da: xarray.DataArray | pd.Series
        DataArray or series containing elevation differences to plot
    ax: matplotlib.axes._axes.Axes
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
    if isinstance(da, xr.DataArray):
        n, bins, _ = da.plot.hist(
            ax=ax, bins=100, range=(dmin, dmax), color="gray", **kwargs
        )
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


@depends_on_optional("matplotlib")
def compare_dems(
    dem_list: list[xr.Dataset],
    gdf_dict: dict[str, tuple[gpd.GeoDataFrame, str]] | None = None,
    ds_wc: xr.Dataset | None = None,
    elevation_cmap: str = "inferno",
    hillshade: bool = True,
    diff_cmap: str = "RdBu",
    elevation_clim: tuple[float, float] | None = None,
    diff_clim: tuple[float, float] = (-10, 10),
    figsize: tuple[float, float] = (17, 11),
    show: bool = True,
) -> plt.Axes:
    """
    Generalized function to create a standard panel comparing DEMs and altimetry points.
    Where the first row shows elevation maps over hillshades (if hillshade is True)
        altimetry points will be plotted over the hillshade of the first DEM (if hillshade is True)
    The second row shows elevation differences, with the first plot being the Landcover of the scene
        First DEM provided will be the 'source' DEM where all other elevations will be compared to
    The third row shows the histograms of these differences

    NOTE: This figure changes dynamically based on the number of DEMs and GeoDataFrames provided.
    This function requires:
        * All inputs to have the same CRS and to be aligned
        * All DEMs to have variable 'elevation'
        * There is a minimum requirement of 2 DEMs to be passed
        * There is a maximum of 5 total datasets (DEMs + GeoDataFrames) to be passed

    Parameters
    ----------
    dem_list: list
        List of DEM objects. Each DEM should have .elevation and .hillshade attributes.
        The first DEM in this list will be the 'source' DEM and will served as a reference for differencing
        e.g.  [dem_1,dem_2,dem_3] will result in diff_1 = dem_2 - dem_1 and diff_2 = dem_3 - dem_1
    gdf_dict: Optional Dict[str, Tuple[gpd.GeoDataFrame, str]]
        Dictionary where keys are subplot names and values are (GeoDataFrame, column_name) pairs.
        Where column_name denotes the column name that contains elevation values of interest
            (e.g. h_li for ICESat-2 ATL06)
    ds_wc: Optional xr.Dataset
        WorldCover dataset to plot. If this is not provided, plot_worldcover_custom_ax
        will be called to create a ESA WorldCover raster. If provided, this should be
        the same CRS and aligned with your other inputs
    elevation_cmap: Optional str
        Colormap for elevation. Default inferno
    hillshade: Optional bool
        Whether to plot hillshades under your DEMs. Default True.
        If True, your DEMs should have a .hillshade attribute
    diff_cmap: Optional str
        Colormap for elevation differences. Default RdBu
    elevation_clim: Optional Tuple[float, float]
        Tuple for elevation color limits
    diff_clim: Optional Tuple[float, float]
        Tuple for difference color limits. Default (-10,10)
    figsize: Optional Tuple[int, int]
        Figure size tuple. Will be dynamically adjusted based on number of DEMs and GeoDataFrames
        so you typically don't want to adjust this
    show: bool
        Whether to display the plot. Implemented for testing purposes
        Set False for tests so plots don't display

    Returns
    ----------
        Dictionary containing the axes objects for each subplot
    """
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
    fig, axd = plt.subplot_mosaic(
        mosaic, gridspec_kw=gs_kw, figsize=figsize, layout="constrained"
    )

    def style_ax(
        ax: matplotlib.axes.Axes,
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
            geometry=[box(bounds[0], bounds[1], bounds[2], bounds[3])],
            crs=dem_list[0].rio.crs,
        )
        gf_search = gf_search.to_crs(epsg=4326)
        gf_wc = search(
            dataset="worldcover",
            intersects=gf_search,
            datetime=["2020"],
        )
        try:
            ds_wc = to_dataset(
                gf_wc,
                bands=["map"],
                aoi=gf_search,
                mask=True,
                resolution=0.00081,  # ~90m
            )
        except Exception:
            # mask=True will fail on very small aois (e.g.20mx20m)
            ds_wc = to_dataset(
                gf_wc,
                bands=["map"],
                aoi=gf_search,
                resolution=0.00081,  # ~90m
            )
        ds_wc = ds_wc.rio.reproject(dem_list[0].rio.crs)

    # TOP ROW: Elevation maps
    # ---------------------------
    # raster elevs
    for i, dem in enumerate(dem_list):
        dem_key = f"dem_{i}"
        plot_dem(
            dem.elevation.squeeze(),
            axd[dem_key],
            elevation_cmap,
            da_hillshade=dem.hillshade if hillshade else None,
            add_colorbar=False,
            vmin=(elevation_clim[0] if elevation_clim else None),
            vmax=(elevation_clim[1] if elevation_clim else None),
            title=dem_key,
            alpha=0.5 if hillshade else 1.0,
        )
        style_ax(axd[dem_key])
        if i == n_columns - 1:  # Add colorbar to the last GeoDataFrame plot
            # images[1] since [0] would be the hillshade
            plt.colorbar(axd[dem_key].images[1], ax=axd[dem_key], label="Elevation (m)")

    # point elevs
    # need this for MyPy to handle the chance the gdf_dict is None
    # despite the logic in the loop "if gdf_dict else []"
    if gdf_dict:
        # Plot GeoDataFrames (e.g., IS2 and GEDI)
        for key, (gdf, column) in gdf_dict.items() if gdf_dict else []:
            plot_altimeter_points(
                gdf,
                axd[key],
                column=column,
                title=key,
                da_hillshade=dem_list[0].hillshade if hillshade else None,
                vmin=(elevation_clim[0] if elevation_clim else None),
                vmax=(elevation_clim[1] if elevation_clim else None),
                alpha=0.5 if hillshade else 1.0,
                s=(figsize[0] * figsize[1])
                / n_columns**2,  # Dynamic sizing based on figure dimensions
            )
            style_ax(axd[key])
            if key == list(gdf_dict.keys())[-1]:  # Add colorbar to the last gdf plot
                plt.colorbar(
                    axd[key].collections[0], ax=axd[key], label="Elevation (m)"
                )

    # MIDDLE ROW: Differences
    # ---------------------------
    plot_worldcover_custom_ax(
        ds_wc.isel(time=0),
        ax=axd["worldcover"],
        cax=axd["wc_legend"],
        add_colorbar=True,
        add_labels=False,
    )
    style_ax(axd["worldcover"], aspect="equal")
    axd["worldcover"].set_title("ESA WorldCover 2020")

    # raster diffs
    for i, dem in enumerate(dem_list[1:], start=1):
        diff = get_elev_diff(dem, dem_list[0])
        diff_key = f"diff_{i}"
        diff.elev_diff.squeeze().plot.imshow(
            ax=axd[diff_key],
            cmap=diff_cmap,
            add_colorbar=False,
            add_labels=False,
            vmin=diff_clim[0],
            vmax=diff_clim[1],
        )
        style_ax(axd[diff_key], facecolor="black", aspect="equal")
        axd[diff_key].set_title(f"dem_{i} minus dem_0")
        if i == n_columns - 1:  # Add colorbar to the last gdf plot
            plt.colorbar(
                axd[diff_key].images[0],
                ax=axd[diff_key],
                label="Elevation Difference (m)",
            )
    # gdf diffs
    # need this for MyPy to handle the chance the gdf_dict is None
    # despite the logic in the loop "if gdf_dict else []"
    if gdf_dict:
        for _i, (key, (gdf, column)) in enumerate(gdf_dict.items() if gdf_dict else []):
            diff_key = f"diff_{key}"
            gf_diff = get_elev_diff(gdf, dem_list[0], source_col=column)
            plot_altimeter_points(
                gf_diff,
                axd[diff_key],
                "elev_diff",
                cmap=diff_cmap,
                vmin=diff_clim[0],
                vmax=diff_clim[1],
                alpha=1.0,
                title=f"{column} minus dem_0",
                s=(figsize[0] * figsize[1])
                / n_columns**2,  # Dynamic sizing based on figure dimensions
            )
            style_ax(axd[diff_key], facecolor="black")
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
        plot_diff_hist(get_elev_diff(dem, dem_list[0]).elev_diff, ax=axd[hist_key])
        axd[hist_key].set_title("")

    # gdf hists
    for _i, (key, (gdf, column)) in enumerate(gdf_dict.items() if gdf_dict else []):
        hist_key = f"hist_{key}"
        diff = get_elev_diff(gdf, dem_list[0], source_col=column)
        plot_diff_hist(diff.elev_diff, ax=axd[hist_key])
        axd[hist_key].set_title("")

    if show:
        plt.show()
    return axd


@depends_on_optional("matplotlib")
def boxplot_terrain_diff(
    dem_list: list[xr.Dataset | gpd.GeoDataFrame],
    terrain_v: str = "slope",
    terrain_groups: np.ndarray | None = None,
    ylims: list[float] | None = None,
    title: str = "Elevation Differences (m) by Terrain",
    xlabel: str = "Terrain Groups",
    elev_col: str | None = None,
    show: bool = True,
) -> plt.Axes:
    """
    Plots boxplots of elevation differences by terrain group for two input elevation datasets.
    (e.g. boxplots for elevation differences over varying slopes for two DEMs)
    This also plots the counts of the grouped elevation differences

    NOTE:
        - This function assumes that the datasets in dem_list are in the same CRS and aligned
        - This function can also work with point geometry GeoDataFrames (e.g. ICESat-2 points)
            - If using a GeoDataFrame, you must also provide the elev_col parameter
            which is the column name containing elevation values you wish to compare
        - This function requires there to be EXACTLY 2 datasets in the dem_list.
        - The first dataset in dem_list MUST be a xr.Dataset with an 'elevation' variable
        and the corresponding terrain_v variable (e.g. 'slope')

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
        Array defining the edges of terrain groups.
    ylims : list
        The y-axis limits for the plot. If None, then the y-axis limits will be
        automatically adjusted to fit the whiskers of each boxplot
    title : str
        The title of the plot.
    xlabel : str
        The label for the x-axis.
    elev_col : str
        The name of the column containing elevation values to difference in
        the GeoDataFrame if a GeoDataFrame is provided
    show: bool
        Whether to display the plot. Implemented for testing purposes
        Set False for tests so plots don't display
    Returns
    -------
    plt.Axes
        The matplotlib axes object containing the plot.
    """
    if len(dem_list) != 2:
        msg_len = "dem_list must contain exactly two datasets"
        raise ValueError(msg_len)
    # ruff B008 if default is set to np.arange(0, 46, 1)
    if terrain_groups is None:
        terrain_groups = np.arange(0, 46, 1)
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
    fig, ax = plt.subplots(figsize=(14, 6))
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
    ax.set_ylabel("Elevation Difference (m)", fontsize=12)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_title(title, fontsize=14)

    lines1, labels1 = ax.get_legend_handles_labels()
    ax.legend(lines1, labels1, loc="best", fontsize=10)

    if show:
        plt.tight_layout()
        plt.show()

    return ax
