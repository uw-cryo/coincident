from __future__ import annotations

import geopandas as gpd
import numpy as np
import xarray as xr


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


def hillshade(
    da: xr.DataArray, azi: int = 45, alt: int = 45
) -> tuple[tuple[str, ...], np.ndarray]:
    """
    Create hillshade array from DEM.
    Pass with dem['hillshade'] = coincident.plot.utils.hillshade(dem.elevation)

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
