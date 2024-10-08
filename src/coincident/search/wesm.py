"""
Spatial metadata for the 3D Elevation Program (3DEP) is available by work unit in the Work Unit Extent Spatial Metadata (WESM).
https://www.usgs.gov/ngp-standards-and-specifications/wesm-data-dictionary
"""

from __future__ import annotations

from typing import Any

import fsspec
import s3fs
from geopandas import GeoDataFrame, read_file, read_parquet
from pandas import Series, Timedelta, Timestamp

from coincident.datasets import usgs
from coincident.overlaps import subset_by_temporal_overlap

defaults = usgs.ThreeDEP()
wesm_gpkg_url = defaults.search


def simplify_footprint(gf: GeoDataFrame, tolerance: float = 0.1) -> None:
    """complexity of WESM footprints causes slow STAC searches"""
    # parts of a simplified geometry will be no more than `tolerance` distance from the original
    gf = gf.simplify(tolerance)


def stacify_column_names(gf: GeoDataFrame) -> GeoDataFrame:
    """
    Renames WESM column names to match STAC (SpatioTemporal Asset Catalog) conventions

    Parameters
    ----------
    gf : GeoDataFrame
        The input GeoDataFrame with columns 'collect_start' and 'collect_end'.

    Returns
    -------
    GeoDataFrame
        The modified GeoDataFrame with renamed columns and additional 'datetime' and 'duration' columns.

    Notes
    -----
        - 'start_datetime': Renamed from 'collect_start'.
        - 'end_datetime': Renamed from 'collect_end'.
        - 'datetime': Midpoint of 'start_datetime' and 'end_datetime'.
        - 'duration': Duration in days between 'start_datetime' and 'end_datetime'.
    """

    name_map = {
        "collect_start": "start_datetime",
        "collect_end": "end_datetime",
    }
    gf = gf.rename(columns=name_map)
    duration = gf.end_datetime - gf.start_datetime
    gf["datetime"] = gf.start_datetime + duration / 2
    gf["dayofyear"] = gf.datetime.dt.dayofyear
    gf["duration"] = duration.dt.days

    return gf


def search_convex_hulls(
    url: str = "https://github.com/uw-cryo/stv-aux-data/releases/download/v0.1/WESM-chulls.geoparquet",
    intersects: GeoDataFrame | None = None,
    search_start: Timestamp | None = None,
    search_end: Timestamp | None = None,
    # **kwargs: dict[str, Any] | None,
) -> GeoDataFrame:
    """
    Search GeoParquet WESM index of convex hulls within a time range and spatial intersection.

    Parameters
    ----------
    url : str, optional
        URL to the GeoParquet file containing the convex hulls, by default "https://github.com/uw-cryo/stv-aux-data/releases/download/v0.1/WESM-chulls.geoparquet".
    intersects : GeoDataFrame, optional
        A GeoDataFrame to spatially intersect with the convex hulls, by default None.
    search_start : Timestamp, optional
        The start time for the temporal search, by default None.
    search_end : Timestamp, optional
        The end time for the temporal search, by default None.

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing the convex hull geometries and FIDs in EPSG:4326
    """
    with fsspec.open(url) as f:
        gf = read_parquet(f)

    gf = stacify_column_names(gf)
    if intersects is not None:
        gf = gf[gf.intersects(intersects.geometry.iloc[0])]

    if search_start is not None and search_end is not None:
        gf = subset_by_temporal_overlap(gf, search_start, search_end, temporal_buffer=0)

    return gf


def load_by_fid(
    fids: list[int], url: str = wesm_gpkg_url, **kwargs: dict[str, Any] | None
) -> GeoDataFrame:
    """
    Extracts full geometry from remote WESM GPKG based on Feature ID (FID)

    Parameters
    ----------
    fids : list(int)
        List of Feature IDs (fids) to extract from GPKG
    url : str
        The URL or file path to the GeoPackage (GPKG) file.

    **kwargs : dict[str, Any], optional
        Additional keyword arguments to pass to :meth:`~geopandas.read_file`

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing the data from the GeoPackage file with column names
        modified by :meth:`stacify_column_names`.
    """
    # Format SQL: # special case for (a) not (a,)
    query = f"fid in ({fids[0]})" if len(fids) == 1 else f"fid in {*fids,}"

    gf = read_file(
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
    gf: GeoDataFrame, start_col: str = "START_TIME", end_col: str = "END_TIME"
) -> GeoDataFrame:
    """
    Convert swath GPS time columns to datetimes in a GeoDataFrame

    Parameters
    ----------
    gf : GeoDataFrame
        The GeoDataFrame containing the swath time columns.
    start_col : str, optional
        The name of the column containing the start times, by default 'START_TIME'.
    end_col : str, optional
        The name of the column containing the end times, by default 'END_TIME'.

    Returns
    -------
    GeoDataFrame
        The GeoDataFrame with additional 'start_ 'dayofyear' column.
    """
    gf["collect_start"] = Timestamp("1980-01-06") + gf[start_col].apply(
        lambda x: Timedelta(seconds=x + 1e9)
    )
    gf["collect_end"] = Timestamp("1980-01-06") + gf[end_col].apply(
        lambda x: Timedelta(seconds=x + 1e9)
    )
    # Assuming no flights spanning multiple days / happening at midnight :)

    return stacify_column_names(gf)


def get_swath_polygons(row: Series) -> GeoDataFrame:
    """
    Retrieve swath polygons from a remote URL and convert swath time to datetime.

    Parameters
    ----------
    row : pd.Series
        A pandas Series containing the project and workunit information needed to construct the URL.

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing the swath polygons with swath time converted to datetime.

    Notes
    -----
        - Swath polygons only available for data collected after 2020
        - Usually takes ~10 to 60s to load one of these shapefiles
    """
    s3 = s3fs.S3FileSystem(anon=True)
    url = f"prd-tnm/StagedProducts/Elevation/metadata/{row.project}/{row.workunit}/spatial_metadata/USGS/**/*.shp"
    try:
        files = s3.glob(url)
        swath_poly_key = next(x for x in files if "swath" in x.lower())
    except StopIteration:
        missing_link = f"http://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/metadata/{row.project}/{row.workunit}"
        message = f"Unable to find swath polygon shapefile in {missing_link}"
        raise FileNotFoundError(message) from None

    uri = swath_poly_key.replace("prd-tnm", "https://prd-tnm.s3.amazonaws.com")
    gf = read_file(uri)

    gf = swathtime_to_datetime(gf)

    # Swath polygons likely have different CRS (EPSG:6350), so reproject
    return gf.to_crs("EPSG:4326")
