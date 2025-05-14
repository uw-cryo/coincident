from __future__ import annotations

import warnings
from typing import Any

import geopandas as gpd
import numpy as np
import xarray as xr

try:
    from osgeo import gdal, gdal_array

    gdal.UseExceptions()
except ImportError:
    warnings.warn(
        "gdal python bindings not found. `pip install gdal` for gdaldem functions",
        stacklevel=2,
    )

from coincident._utils import depends_on_optional


# TODO: better to split to separate functions, each with unique return type
def get_elev_diff(
    source: xr.Dataset | gpd.GeoDataFrame,
    reference: xr.DataArray,
    source_col: str | None = None,
    diff_col_name: str | None = "elev_diff",
) -> xr.Dataset | gpd.GeoDataFrame:
    """
    Calculate elevation differences between source and reference datasets.
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

    Notes
    -----
    This function does not reproject data (it assumes source and reference have the same CRS).
    """
    if isinstance(source, xr.Dataset):
        source_elev_diff = source.elevation - reference.elevation
        source[diff_col_name] = source_elev_diff

        return source

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

    msg_type_error = "source must be xr.Dataset or gpd.GeoDataFrame"
    raise TypeError(msg_type_error)


@depends_on_optional("osgeo")
def gdaldem(
    da: xr.DataArray, subcommand: str = "hillshade", **kwargs: Any
) -> xr.DataArray:
    """
    Use GDAL python bindings to perform gdaldem operations

    Parameters
    ----------
    da: xarray.DataArray
        Dataarray containing elevation values
    subcommand: str
        'hillshade' (default), 'aspect', 'slope'
    kwargs: dict
        https://gdal.org/en/stable/api/python/utilities.html#osgeo.gdal.DEMProcessingOptions

    Returns
    -------
    da: xarray.DataArray
        XArray dataset with result of running GDAL algorithm
    """
    # Ensure we're working with a 2D Array
    elevation = da.squeeze()
    src = gdal_array.OpenArray(elevation.to_numpy())
    geotransform = da.rio.transform().to_gdal()
    src.SetGeoTransform(geotransform)

    ds_gdal = gdal.DEMProcessing("", src, subcommand, format="MEM", **kwargs)
    numpy_result = ds_gdal.ReadAsArray()

    if subcommand != "hillshade":
        numpy_result[numpy_result == -9999] = np.nan

    del ds_gdal

    return elevation.copy(data=numpy_result)
