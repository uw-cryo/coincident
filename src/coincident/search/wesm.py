"""
Spatial metadata for the 3D Elevation Program (3DEP) is available by work unit in the Work Unit Extent Spatial Metadata (WESM).
https://www.usgs.gov/ngp-standards-and-specifications/wesm-data-dictionary
"""

from __future__ import annotations

from importlib import resources
from typing import Any

import pandas as pd
import pyogrio
from geopandas import GeoDataFrame, read_file
from pandas import Timedelta, Timestamp
from shapely.geometry import box

from coincident.datasets import usgs
from coincident.overlaps import subset_by_temporal_overlap

# Geopandas S3 Client
pyogrio.set_gdal_config_options(
    {
        "AWS_NO_SIGN_REQUEST": True,
        "GDAL_PAM_ENABLED": False,
        "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".gpkg .geojson .json .tif .shp .shx .dbf .prj .cpg",
    }
)

swath_polygon_csv = resources.files("coincident.search") / "swath_polygons.csv"

defaults = usgs.ThreeDEP()
wesm_gpkg_url = defaults.search


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
    # Add collection name to match STAC (used to identify geodataframe contents in other functions)
    gf["collection"] = "3DEP"
    gf["start_datetime"] = pd.to_datetime(
        gf["start_datetime"], format="ISO8601", errors="coerce"
    )
    gf["end_datetime"] = pd.to_datetime(
        gf["end_datetime"], format="ISO8601", errors="coerce"
    )
    duration = gf.end_datetime - gf.start_datetime
    gf["datetime"] = gf.start_datetime + duration / 2
    gf["dayofyear"] = gf.datetime.dt.dayofyear
    gf["duration"] = duration.dt.days

    return gf


def read_wesm_csv(url: str = wesm_gpkg_url) -> GeoDataFrame:
    """
    Read WESM metadata from a remote CSV file.

    The CSV contains up to date metadata and a mapping of workunit to feature ID (FID)
    in the main GPKG file.

    Parameters
    ----------
    url : str, optional
        The URL or file path to the WESM  index file

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing the metadata from the CSV file.
    """
    # TODO: Cache CSV locally, so subsequent reads are faster
    wesm_csv = url.replace(".gpkg", ".csv")
    df = pd.read_csv(wesm_csv)
    df.index += 1
    df.index.name = "fid"
    return df


def search_bboxes(
    url: str = wesm_gpkg_url,
    intersects: GeoDataFrame | None = None,
    search_start: Timestamp | None = None,
    search_end: Timestamp | None = None,
    # **kwargs: dict[str, Any] | None,
) -> GeoDataFrame:
    """
    Search WESM.gpkg bounding boxes for spatial and temporal overlap

    Parameters
    ----------
    url : str, optional
        The URL or file path to the GeoPackage (GPKG) file.
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
    # NOTE: much faster to JUST read bboxes, not full geometry or other columns
    sql = "select * from rtree_WESM_geometry"
    df = pyogrio.read_dataframe(
        url, sql=sql
    )  # , use_arrow=True... arrow probably doesn;t matter for <10000 rows?
    bboxes = df.apply(lambda x: box(x.minx, x.miny, x.maxx, x.maxy), axis=1)
    gf = (
        GeoDataFrame(df, geometry=bboxes, crs="EPSG:4326")
        .rename(columns={"id": "fid"})
        .set_index("fid")
    )

    df = read_wesm_csv(url)
    gf = gf.merge(df, right_index=True, left_index=True)

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
    # Reading a remote WESM by specific FIDs is fast
    query = f"fid in ({fids[0]})" if len(fids) == 1 else f"fid in {(*fids,)}"

    gf = read_file(
        url,
        where=query,
        **kwargs,
        # mask=mask, # spatial subset intolerably slow for remote GPKG...
        # NOTE: worth additional dependencies for speed?
        # Only faster I think if GPKG is local & large https://github.com/geopandas/pyogrio/issues/252
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


def get_swath_polygons(
    workunit: str,
) -> GeoDataFrame:
    """
    Retrieve swath polygons from a remote URL and convert swath time to datetime.

    Parameters
    ----------
    fid : int
        A pandas Series containing WESM feature ID (FID)

    Returns
    -------
    GeoDataFrame
        A GeoDataFrame containing the swath polygons with swath time converted to datetime.

    Notes
    -----
        - Swath polygons not available for all workunits
    """
    # NOTE: this CSV comes from crawling WESM USGS/spatial_metadata in S3 bucket
    df = pd.read_csv(swath_polygon_csv)

    try:
        url = df[df.workunit == workunit].swathpolygon_link.iloc[0]
    except Exception as e:
        message = f"No swath polygons found for workunit={workunit}"
        raise ValueError(message) from e

    # Actually read from S3!
    gf = read_file(url)

    gf = swathtime_to_datetime(gf)
    # Swath polygons likely have different CRS (EPSG:6350), so reproject
    return gf.to_crs("EPSG:4326")


# NOTE it says that 3DEP elevation products were "lastUpdatedDate": "Sep 11, 2019"?
# https://tnmaccess.nationalmap.gov/api/v1/datasets?
def query_tnm_api(
    polygon_str: str,
    tnmdataset: str = "Digital Elevation Model (DEM) 1 meter",
    prodFormats: str = "GeoTIFF",
) -> list[dict[str, Any]]:
    """
    Query the TNM API for USGS DEM products that intersect the given polygon.

    Parameters:
      polygon_str (str): Polygon coordinates string in the required API format.
      tnmdataset (str): The dataset name for the TNM API (default is "Digital Elevation Model (DEM) 1 meter").
      prodFormats (str): The product format filetype to search (default is "GeoTIFF").

    Returns:
      list: A list of JSON items returned from the TNM API.
    """
    # params is a dict[str, object], but requests expects Mapping[str, Any]
    import requests  # type: ignore[import-untyped]

    items = []
    offset = 0
    url = "https://tnmaccess.nationalmap.gov/api/v1/products"

    while True:
        params = {
            "datasets": tnmdataset,
            "polygon": polygon_str,
            "prodFormats": prodFormats,
            "outputFormat": "JSON",
            "max": 200,
            "offset": offset,
        }
        response = requests.get(url, params=params)
        if response.status_code != 200:
            msg_status_code = (
                f"TNM API request failed with status code {response.status_code}"
            )
            raise RuntimeError(msg_status_code)

        data = response.json()
        batch = data.get("items", [])
        if not batch:
            # no more products -> exit loop
            break

        items.extend(batch)
        offset += len(batch)

    return items
