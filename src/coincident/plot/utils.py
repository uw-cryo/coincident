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


def sample_dem_at_points(
    da_dem: xr.DataArray,
    gdf_points: gpd.GeoDataFrame,
    diff_col: str | None = None,
) -> pd.DataFrame:
    """
    Sample DEM at point locations from GeoDataFrame using Xarray Advanced Interpolation (bilinear interpolation)
    If the GeoDataFrame column name is provided (e.g. h_li for IS2 ATL06), also calculate
    difference by subtracting the sampled raster elevation.

    Parameters
    ----------
    da_dem: xr.DataArray
        DEM as xarray DataArray with only 'x' and 'y' spatial dimensions
    gdf_points: gpd.GeoDataFrame
        GeoDataFrame with point geometries to sample DEM at
    diff_col: Optional[str]
        Column in GeoDataFrame to difference [gdf_points[diff_col] - DEM]

    Returns
    -------
    gpd.pd.GeoDataFrame
        GeoDataFrame with sampled DEM elevations and optional 'elev_diff' column
    """
    assert isinstance(da_dem, xr.DataArray), "da_dem must be an xarray.DataArray"
    da_points = gdf_points.get_coordinates().to_xarray()
    # Ensure we don't return a multi-index
    da_dem = da_dem.squeeze()
    samples = (
        da_dem.drop_vars(["band", "spatial_ref"])
        .interp(da_points)
        .to_dataframe(name="dem_elevation")
    )

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
