"""
USGS 3DEP Search via WESM.gpkg
"""

from __future__ import annotations

from typing import Any

import fsspec
import geopandas as gpd

from coincident.datasets import usgs
from coincident.overlaps import subset_by_temporal_overlap

defaults = usgs.ThreeDEP()
wesm_gpkg_url = defaults.search


def simplify_footprint(gf: gpd.GeoDataFrame, tolerance: float = 0.1) -> None:
    """complexity of WESM footprints causes slow STAC searches"""
    # parts of a simplified geometry will be no more than `tolerance` distance from the original
    gf = gf.simplify(tolerance)


def stacify_column_names(gf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """match STAC items"""
    name_map = {
        "collect_start": "start_datetime",
        "collect_end": "end_datetime",
    }
    gf = gf.rename(columns=name_map)
    duration = gf.end_datetime - gf.start_datetime
    gf["datetime"] = gf.start_datetime + duration / 2
    gf["duration"] = duration.dt.days

    return gf


def search_convex_hulls(
    url: str = "https://github.com/uw-cryo/stv-aux-data/releases/download/v0.1/WESM-chulls.geoparquet",
    intersects: gpd.GeoDataFrame | None = None,
    search_start: gpd.pd.Timestamp | None = None,
    search_end: gpd.pd.Timestamp | None = None,
    # **kwargs: dict[str, Any] | None,
) -> gpd.GeoDataFrame:
    with fsspec.open(url) as f:
        gf = gpd.read_parquet(f)

    gf = stacify_column_names(gf)
    if intersects is not None:
        gf = gf[gf.intersects(intersects.geometry.iloc[0])]

    if search_start is not None and search_end is not None:
        gf = subset_by_temporal_overlap(gf, search_start, search_end, temporal_buffer=0)

    return gf


def load_wesm_by_fid(
    fids: list[int], url: str = wesm_gpkg_url, **kwargs: dict[str, Any] | None
) -> gpd.GeoDataFrame:
    """
    Extracts full geometry from remote GPKG based on feature id
    Parameters
    ----------
    fids : list(int)
        List of Feature IDs (fids) to extract from GPKG
    url : str
        The URL or file path to the GeoPackage (GPKG) file.

    **kwargs : dict[str, Any], optional
        Additional keyword arguments to pass to `gpd.read_file`.
    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the data from the GeoPackage file with column names
        modified by `stacify_column_names`.
    Notes
    -----
    - The function sets the coordinate reference system (CRS) to "EPSG:4326" and allows
      overriding the existing CRS.
    - Consider using the 'pyogrio' engine and 'pyarrow' for faster reading if the GPKG
      file is local. Refer to https://github.com/geopandas/pyogrio/issues/252 for more
      details.
    """
    # Format SQL: # special case for (a) not (a,)
    query = f"fid in ({fids[0]})" if len(fids) == 1 else f"fid in {*fids,}"

    gf = gpd.read_file(
        url,
        where=query,
        **kwargs,
        # mask=mask, # spatial subset intolerably slow for remote GPKG...
        # NOTE: worth additional dependencies for speed?
        # Only faster I think if GPKG is local https://github.com/geopandas/pyogrio/issues/252
        # engine='pyogrio',
        # pyarrow=True,
    )

    # For the purposes of footprint polygon search, just ignore datum (NAD83)
    gf = gf.set_crs("EPSG:4326", allow_override=True)

    return stacify_column_names(gf)


def swathtime_to_datetime(
    gf: gpd.GeoDataFrame, start_col: str = "START_TIME", end_col: str = "END_TIME"
) -> gpd.GeoDataFrame:
    """
    Convert swath time columns to datetime in a GeoDataFrame.
    This function converts the swath time columns (start and end times) from GPS time to datetime format.
    The conversion is accurate to within tens of seconds (ignores leap seconds).
    Parameters
    ----------
    gf : gpd.GeoDataFrame
        The GeoDataFrame containing the swath time columns.
    start_col : str, optional
        The name of the column containing the start times, by default 'START_TIME'.
    end_col : str, optional
        The name of the column containing the end times, by default 'END_TIME'.
    Returns
    -------
    gpd.GeoDataFrame
        The GeoDataFrame with the converted datetime columns and an additional 'dayofyear' column.
    """
    gf[start_col] = gpd.pd.Timestamp("1980-01-06") + gf[start_col].apply(
        lambda x: gpd.pd.Timedelta(seconds=x + 1e9)
    )
    gf[end_col] = gpd.pd.Timestamp("1980-01-06") + gf[end_col].apply(
        lambda x: gpd.pd.Timedelta(seconds=x + 1e9)
    )
    # Assuming no flights spanning multiple days / happening at midnight :)
    gf["dayofyear"] = gf[start_col].dt.dayofyear
    return gf


def get_swath_polygons(row: gpd.GeoSeries) -> gpd.GeoDataFrame:
    """
    Retrieve swath polygons from a remote URL and convert swath time to datetime.
    Parameters
    ----------
    row : gpd.GeoSeries
        A GeoSeries containing the project and workunit information needed to construct the URL.
    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the swath polygons with swath time converted to datetime.
    """
    base_url = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/metadata"
    # NOTE: does naming convention depend on 'spec' e.g. (USGS Lidar Base Specification 2020 Revision A)?
    if row.start_datetime < gpd.pd.Timestamp("2020"):
        url = f"{base_url}/{row.project}/{row.workunit}/spatial_metadata/USGS/USGS_{row.workunit}_SWATH_POLYGON.shp"
    else:
        url = f"{base_url}/{row.project}/{row.workunit}/spatial_metadata/USGS/USGS_{row.workunit}_SwathPolygon.shp"
    print("Attempting to retrieve:", url)  # noqa: T201
    gf = gpd.read_file(url)

    gf = swathtime_to_datetime(gf)

    # Swath polygons likely have different CRS (EPSG:6350), so reproject
    return gf.to_crs("EPSG:4326")
