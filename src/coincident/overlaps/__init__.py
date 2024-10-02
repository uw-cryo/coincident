"""
Functions for refining search results based on spatial and temporal overlaps
"""

from __future__ import annotations

from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry
from shapely.geometry.polygon import orient as _orient


def subset_by_maximum_duration(
    gf: gpd.GeoDataFrame,
    max_duration: int = 60,  # days
) -> gpd.GeoDataFrame:
    """
    Subset a GeoDataFrame by a maximum duration.

    Parameters
    ----------
    gf : _gpd.GeoDataFrame
        The input GeoDataFrame containing a 'duration' column.
    max_duration : int, optional
        The maximum duration in days to filter the GeoDataFrame. Default is 60 days.

    Returns
    -------
    _gpd.GeoDataFrame
        A GeoDataFrame filtered to include only rows where the 'duration' is less than or equal to the specified maximum duration.
    """
    max_duration = pd.Timedelta(days=max_duration)
    # NOTE: keep indexes?
    return gf.loc[gf.duration <= max_duration, :].reset_index()


def subset_by_temporal_overlap(
    gf: gpd.GeoDataFrame,
    start_datetime: pd.Timestamp,
    end_datetime: pd.Timestamp,
    temporal_buffer: int = 14,
) -> gpd.GeoDataFrame:
    """
    Subset a DataFrame based on temporal overlap.

    Parameters:
    gf (DataFrame): The input DataFrame with a datetime column.
    start_datetime (Timestamp): The start datetime for the overlap.
    end_datetime (Timestamp): The end datetime for the overlap.
    temporal_buffer (int): The buffer to apply to the start and end datetimes.

    Returns:
    DataFrame: A subset of the input DataFrame with rows that overlap temporally.
    """
    temporal_buffer = pd.Timedelta(days=temporal_buffer)
    # Apply temporal buffer to collection start and end dates
    buffered_start = gf.start_datetime - temporal_buffer
    buffered_end = gf.end_datetime + temporal_buffer

    results_intervals = gpd.pd.IntervalIndex.from_arrays(
        buffered_start, buffered_end, closed="both"
    )
    search_interval = gpd.pd.Interval(start_datetime, end_datetime, closed="both")
    keep = results_intervals.overlaps(search_interval)

    return gf.loc[keep, :].reset_index()


def geographic_area(gf: gpd.GeoDataFrame) -> gpd.pd.Series:
    """
    Calculate the geographic area of each polygon in a GeoDataFrame.
    This function computes the area of polygons in a GeoDataFrame that has a geographic coordinate system.
    It handles both single polygons and multipolygons by summing the areas of individual polygons.
    Parameters
    ----------
    gf : gpd.GeoDataFrame
        A GeoDataFrame containing the geometries for which the area needs to be calculated. The GeoDataFrame
        must have a geographic coordinate system (latitude and longitude).
    Returns
    -------
    pd.Series
        A Pandas Series containing the area of each polygon in the input GeoDataFrame.
    Raises
    ------
    TypeError
        If the GeoDataFrame does not have a geographic coordinate system.
    References
    ----------
    - https://gis.stackexchange.com/questions/413349/calculating-area-of-lat-lon-polygons-without-transformation-using-geopandas
    - https://pyproj4.github.io/pyproj/stable/api/geod.html
    """
    if not gf.crs and gf.crs.is_geographic:
        msg = "geodataframe should have geographic coordinate system"
        raise TypeError(msg)

    geod = gf.crs.get_geod()

    # TODO: sort out numpy typig https://github.com/python/mypy/issues/5480
    def area_calc(geom: shapely.geometry) -> Any:
        if geom.geom_type not in ["MultiPolygon", "Polygon"]:
            return np.nan

        # For MultiPolygon do each separately
        if geom.geom_type == "MultiPolygon":
            return np.sum([area_calc(p) for p in geom.geoms])

        # orient to ensure a counter-clockwise traversal.
        # geometry_area_perimeter returns (area, perimeter)
        return geod.geometry_area_perimeter(_orient(geom, 1))[0]

    return gf.geometry.apply(area_calc)
