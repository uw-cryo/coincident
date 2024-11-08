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
    gpd.GeoDataFrame
        A GeoDataFrame filtered to include only rows where the 'duration' is less than or equal to the specified maximum duration.
    """
    max_duration = pd.Timedelta(days=max_duration)
    # NOTE: keep indexes?
    return gf.loc[gf.duration <= max_duration, :]


def subset_by_temporal_overlap(
    gf: gpd.GeoDataFrame,
    start_datetime: pd.Timestamp,
    end_datetime: pd.Timestamp,
    temporal_buffer: int = 14,
) -> gpd.GeoDataFrame:
    """
    Subset GeoDataFrame with time intervals by a single temporal span.

    Parameters
    ----------
    gf : gpd.GeoDataFrame
        The input GeoDataFrame with 'start_datetime' and 'end_datetime' columns
    start_datetime : pd.Timestamp
        The start datetime for the overlap.
    end_datetime : pd.Timestamp
        The end datetime for the overlap.
    temporal_buffer : int, optional
        The buffer in days to apply to the start and end datetimes (default is 14).

    Returns
    -------
    gpd.GeoDataFrame
        A subset of the input GeoDataFrame with rows that overlap temporally with the specified span.
    """
    temporal_buffer = pd.Timedelta(days=temporal_buffer)
    # Apply temporal buffer to collection start and end dates
    buffered_start = gf.start_datetime - temporal_buffer
    buffered_end = gf.end_datetime + temporal_buffer

    results_intervals = gpd.pd.IntervalIndex.from_arrays(
        buffered_start, buffered_end, closed="both"
    )
    search_interval = gpd.pd.Interval(start_datetime, end_datetime, closed="both")
    keep = results_intervals.overlaps(search_interval)  # pylint: disable=no-member

    return gf.loc[keep, :]


def geographic_area(gf: gpd.GeoDataFrame) -> gpd.pd.Series:
    """
    Calculate the geographic area of each polygon in a GeoDataFrame in KM^2

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

    # TODO: sort out numpy typing https://github.com/python/mypy/issues/5480
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


def subset_by_minimum_area(
    gf: gpd.GeoDataFrame,
    min_area: float = 20.0,  # square kilometers
) -> gpd.GeoDataFrame:
    """
    Subset a GeoDataFrame by a minimum area threshold.

    Parameters
    ----------
    gf : gpd.GeoDataFrame
        The input GeoDataFrame containing geographic features.
    min_area : int, optional
        The minimum area threshold in square kilometers. Features with an area
        less than or equal to this value will be included in the subset. Default is 20.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing only the features with an area less than or equal
        to the specified minimum area.
    """
    areas = geographic_area(gf) * 1e-6

    return gf.loc[areas >= min_area, :]
