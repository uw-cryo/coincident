"""
Functions for refining search results based on spatial and temporal overlaps
"""

from __future__ import annotations

import geopandas as _gpd
import pandas as _pd  # noqa: ICN001


def subset_by_maximum_duration(
    gf: _gpd.GeoDataFrame,
    max_duration: int = 60,  # days
) -> _gpd.GeoDataFrame:
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
    max_duration = _pd.Timedelta(max_duration)
    # NOTE: keep indexes?
    return gf.loc[gf.duration <= max_duration, :].reset_index()


def subset_by_temporal_overlap(
    gf: _gpd.GeoDataFrame,
    start_datetime: _pd.Timestamp,
    end_datetime: _pd.Timestamp,
    temporal_buffer: int = 14,
) -> _gpd.GeoDataFrame:
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
    temporal_buffer = _pd.Timedelta(temporal_buffer)
    # Apply temporal buffer to collection start and end dates
    buffered_start = gf.start_datetime - temporal_buffer
    buffered_end = gf.end_datetime + temporal_buffer

    results_intervals = _gpd.pd.IntervalIndex.from_arrays(
        buffered_start, buffered_end, closed="both"
    )
    search_interval = _gpd.pd.Interval(start_datetime, end_datetime, closed="both")
    keep = results_intervals.overlaps(search_interval)

    return gf.loc[keep, :].reset_index()
