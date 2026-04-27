"""
Common functions for matplotlib visualizations
"""

from __future__ import annotations

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import xyzservices
from scipy import stats
from rasterstats import zonal_stats
from affine import Affine


def get_tiles(name: str) -> xyzservices.TileProvider:
    """Return xyzservices TileProvider object for contextily basemaps"""
    return xyzservices.providers.flatten()[name]


def get_aspect(data: xr.DataArray | gpd.GeoDataFrame) -> float | None:
    """Determine best plot aspect ratio based on latitude and CRS for a xr.DataArray or gpd.GeoDataFrame"""
    if isinstance(data, gpd.GeoDataFrame):
        bounds = data.total_bounds
        crs = data.crs
    elif isinstance(data, xr.DataArray):
        bounds = data.rio.bounds(recalc=True)  # in case of xarray ops like coarsen
        crs = data.rio.crs
    else:
        message = (
            "aspect calculation requires an xarray.DataArray or geopandas.GeoDataFrame"
        )
        raise TypeError(message)

    if crs.is_geographic:
        mid_lat = (bounds[1] + bounds[3]) / 2
        aspect = float(1 / np.cos(np.deg2rad(mid_lat)))
    elif crs.is_projected:
        aspect = 1.0
    else:
        aspect = None  # Not sure how matplotlib does default aspect...

    return aspect


def get_haversine_distance(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """Approximate great-circle distance for scale bars"""
    R = 6371000

    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    a = (
        np.sin(dlat / 2) ** 2
        + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    )
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return float(R * c)


def get_scale(ax: plt.Axes) -> float:
    """distance in meters of two points for a given latitude (Y coordinate) which are one full degree of longitude (X) apart
    ax: matplotlib axis object with lon, lat coordinates
    """
    xmin, _ = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    midlat = (ymin + ymax) / 2
    return get_haversine_distance(xmin, midlat, xmin + 1, midlat)


def da_to_rasterstats_inputs(da):
    """
    Extract numpy array, affine transform, and nodata
    from an xarray DataArray for use with rasterstats.

    Parameters
    ----------
    da : xarray.DataArray
        2D array with dimensions (y, x).

    Returns
    -------
    data : np.ndarray
    affine : Affine
    nodata : float or None
    """
    data = da.values

    x = da.x.values
    y = da.y.values
    res_x = float(x[1] - x[0])
    res_y = float(y[1] - y[0])  # typically negative for top-down rasters

    # Affine: upper-left corner of the upper-left pixel
    affine = Affine(res_x, 0, x[0] - res_x / 2,
                    0, res_y, y[0] - res_y / 2)

    nodata = None
    if hasattr(da, "rio") and da.rio.nodata is not None:
        nodata = da.rio.nodata
    elif "_FillValue" in da.attrs:
        nodata = da.attrs["_FillValue"]
    elif "_FillValue" in da.encoding:
        nodata = da.encoding["_FillValue"]

    return data, affine, nodata



def sample_windowed_stats(da, gdf, window_radius=12.5,
                          stats=["mean", "median", "std", "min", "max", "count"],
                          percentiles=None):
    """
    ...
    percentiles : list of int, optional
        Percentile values to compute (e.g., [2, 98]).
    """
    data, affine, nodata = da_to_rasterstats_inputs(da)

    # Build percentile stat keys for rasterstats
    pct_keys = []
    if percentiles:
        pct_keys = [f"percentile_{p}" for p in percentiles]

    all_stats = list(stats) + pct_keys

    gdf_buf = gdf.copy()
    gdf_buf["geometry"] = gdf_buf.geometry.buffer(window_radius)

    zs = zonal_stats(
        gdf_buf,
        data,
        affine=affine,
        nodata=np.nan if nodata is None else nodata,
        stats=all_stats,
    )

    return {s: np.array([z[s] for z in zs], dtype=float) for s in all_stats}

def sample_dem_at_points(
    da_dem: xr.DataArray,
    gdf_points: gpd.GeoDataFrame,
    sampling_method: str = 'bilinear',
    diff_col: str | None = None,
    percentile: int | None = None,
) -> pd.DataFrame:
    """
    ...
    percentile : int, optional
        Percentile to compute within footprint (e.g., 98). Only used when
        sampling_method='percentile'.
    """
    assert isinstance(da_dem, xr.DataArray), "da_dem must be an xarray.DataArray"
    da_points = gdf_points.get_coordinates().to_xarray()
    # Ensure we don't return a multi-index
    da_dem = da_dem.squeeze()
    if sampling_method == 'bilinear':
        samples = (
            da_dem.drop_vars(["band", "spatial_ref"])
            .interp(da_points)
            .to_dataframe(name="dem_elevation")
        )
    elif sampling_method == 'percentile':
        if percentile is None:
            raise ValueError("Must specify `percentile` when sampling_method='percentile'")
        pct_key = f"percentile_{percentile}"
        results = sample_windowed_stats(
            da_dem, gdf_points, window_radius=12.5,
            stats=['count'], percentiles=[percentile],
        )
        count = results['count']
        output = results[pct_key]
        output[count < 30] = np.nan

        samples = gdf_points.get_coordinates()
        samples['dem_elevation'] = output
    else:
        results = sample_windowed_stats(
            da_dem, gdf_points, window_radius=12.5,
            stats=['count', sampling_method],
        )
        count = results['count']
        output = results[sampling_method]
        output[count < 30] = np.nan

        samples = gdf_points.get_coordinates()
        samples['dem_elevation'] = output
    if diff_col is not None:
        samples[diff_col] = gdf_points[diff_col].to_numpy()
        samples["elev_diff"] = (
            gdf_points[diff_col].to_numpy() - samples["dem_elevation"]
        )

    return samples

def calc_stats(s: gpd.pd.Series) -> dict[str, float]:
    """Calculate basic statistics for a GeoPandas Series"""
    s = s.dropna()
    return {
        "Mean": np.mean(s),
        "Median": np.median(s),
        "Std": np.std(s),
        "NMAD": stats.median_abs_deviation(s, scale="normal"),
        "Count": s.size,
        "Min": np.min(s) if not s.empty else np.nan,
        "Max": np.max(s) if not s.empty else np.nan,
    }
