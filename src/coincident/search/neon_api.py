"""
Searching for NEON data with spatial and temporal filters is actually quite difficult.
I couldn't find any STAC catalogs or code to make this easier

https://data.neonscience.org/data-products/DP3.30024.001
https://www.neonscience.org/resources/learning-hub/tutorials/neon-api-intro-requests-py

The process goes as follows:
- Use requests to get some basic metadata for all NEON sites that contain discrete return LiDAR point cloud data
    - Here, we grab the site code, site name, available months, and product URLs
    - Note that grabbing the specific collection time for each flight requires
        a lot more key-value pair digging and i do not believe is worth the effort
    - We grab the 'centroid' (top left of the center-most tile of footprint) here as well
        this is because there is no easy way to search on an input geometry and we can
        use the 'centroid' to make that process easier
- Subset the resulting gdf on the user-input start and end times
- Subset that gdf by only including sites that have 'centroids' within 1 degree (~111km)
    we do this because grabbing the actual footprint geometry is time-consuming.
    the only way to do this dynamically, at least barring the use of
    more unnecessary and unique llibraries/packages from what i found,
    is by making a GET request for each site and then reading in the footprint
    hosted on storage.googleapis.com
- Grab the footprint geometry for each site
    These are hosted on storage.googleapis.com in the form of a shapefile
    that has all of the complex tile footprints as their own row geometries
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
import pyogrio
import requests  # type: ignore[import-untyped]
from pandas import Timestamp
from shapely.geometry import Point


def build_neon_point_gf(sites_url: str) -> gpd.GeoDataFrame:
    """
    Build a GeoDataFrame containing NEON site information, including product URLs and available months.

    Parameters
    ----------
    sites_url : str
        The URL to the NEON API for retrieving site data.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame with site information and 'centroid' geometries. Note that these centroids
        are the top left of the center-most tile for each NEON flight. Also note that there is no
        easy way to pull the flight footprint geometry from NEON.
    """
    sites_request = requests.get(sites_url, timeout=30)
    sites_json = sites_request.json()

    df_neon = pd.concat(
        [
            pd.DataFrame(
                {
                    "id": site["siteCode"],
                    "title": site["siteName"],
                    "start_datetime": dict_lidar["availableMonths"],
                    "end_datetime": dict_lidar["availableMonths"],
                    "product_url": dict_lidar["availableDataUrls"],
                    "latitude": site["siteLatitude"],
                    "longitude": site["siteLongitude"],
                }
            )
            for site in sites_json["data"]
            for dict_lidar in [
                next(
                    (
                        item
                        for item in site["dataProducts"]
                        if item["dataProductCode"] == "DP3.30024.001"
                    ),
                    None,
                )
            ]
            if dict_lidar is not None
        ],
        ignore_index=True,
    )

    df_neon["geometry"] = df_neon.apply(
        lambda row: Point(row["longitude"], row["latitude"]), axis=1
    )
    gf_neon = gpd.GeoDataFrame(df_neon, geometry="geometry", crs="EPSG:4326")

    return gf_neon.drop(columns=["latitude", "longitude"])


def temporal_filter_neon(
    gf_neon: gpd.GeoDataFrame, search_start: Timestamp, search_end: Timestamp
) -> gpd.GeoDataFrame:
    """
    Filter NEON GeoDataFrame by temporal range.

    Parameters
    ----------
    gf_neon : gpd.GeoDataFrame
        The GeoDataFrame containing NEON site data from build_neon_point_gf().
    search_start : pd.Timestamp
        The start of the time range to filter by.
    search_end : pd.Timestamp
        The end of the time range to filter by.

    Returns
    -------
    gpd.GeoDataFrame
        A filtered GeoDataFrame based on the temporal range.
    """
    return gf_neon[
        (
            gf_neon["start_datetime"].apply(lambda x: pd.to_datetime(x + "-01"))
            >= search_start
        )
        & (
            gf_neon["start_datetime"].apply(lambda x: pd.to_datetime(x + "-01"))
            <= search_end
        )
    ].reset_index(drop=True)


def get_neon_bboxes(url: str, fallback_geometry: gpd.GeoSeries) -> gpd.GeoSeries:
    """
    Fetch and return bounding boxes for NEON data products.

    Parameters
    ----------
    url : str
        The URL to the NEON product data.
    fallback_geometry : gpd.GeoSeries
        A fallback geometry to return if footprint cannot be retrieved.

    Returns
    -------
    gpd.GeoSeries
        The bounding box geometry for the NEON data product.
    """
    response = requests.get(url, timeout=30)
    if response.status_code != 200:
        return fallback_geometry
    data = response.json().get("data", {})
    files = data.get("files", [])
    # Find the file that ends with 'merged_tiles.shp' dynamically
    # I really can't think of a better, more-efficient way to do this after dissecting the NEON API
    # The merged tiles shapefiles have similar filepaths, but we're missing some extra metadata
    # that contains another ID needed to hardcode those paths effectively (see 'D02' in commented link below)
    shp_file_url = next(
        (f["url"] for f in files if f["name"].endswith("merged_tiles.shp")), None
    )

    try:
        df_tiles = pyogrio.read_dataframe(shp_file_url).to_crs(4326)
    except:  # noqa: E722
        # TODO: fix provisional footprint grab
        # the provisional (usually synonymous with new) NEON flights don't typically have merged_tiles.shp
        # i tried to iterate over all the different tiles and grab bboxes but too inefficient
        # we'll just return the point geometry for now
        return fallback_geometry

    try:
        return df_tiles.union_all().envelope
    except:  # noqa: E722
        # see below link for example where make_valid() is necessary
        # https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SCBI_3/Metadata/DiscreteLidar/TileBoundary/shps/2019_SCBI_3_merged_tiles.shp
        return df_tiles.geometry.make_valid().union_all().envelope


def search_bboxes(
    intersects: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: Timestamp,
    search_end: Timestamp,
) -> gpd.GeoDataFrame:
    """
    Perform a search for NEON LiDAR metadata and respective bbox footprints. Note that this search
    will take a while if you denote a large aoi or a large time range.

    Parameters
    ----------
    intersects : gpd.GeoDataFrame | gpd.GeoSeries
        The geometry to restrict the search. By default does a global search.
    search_start : pd.Timestamp
        The start of the time range to filter by. By default searches all dates in the catalog (don't do this)
    search_end : pd.Timestamp
        The end of the time range to filter by. By default searches all dates in the catalog (don't do this)

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame with NEON metadata and respective bbox footprints
    """
    if search_start is None and search_end is None:
        search_start = pd.Timestamp(
            "2013-06-01"
        )  # Starting from June, 2005 (first NEON dataset)
        search_end = pd.Timestamp.today()  # Default to today's date
    # note that the above will result in the search taking a very long time
    # if you're searching a larger area
    if intersects is None:
        # Define global bounding box
        intersects = gpd.GeoSeries([Point(-180, -90).buffer(360)])

    sites_url = "http://data.neonscience.org/api/v0/sites"
    # Build NEON point GeoDataFrame
    gf_neon = build_neon_point_gf(sites_url)

    # Apply temporal filter
    gf_neon = temporal_filter_neon(gf_neon, search_start, search_end)

    # Find nearest footprint 'centroids' within the intersection geometry
    gf_neon = (
        gpd.sjoin_nearest(gf_neon, intersects, max_distance=1)
        .drop(columns=["index_right"])
        .reset_index(drop=True)
    )

    # Apply bounding box function
    gf_neon["geometry"] = gf_neon.apply(
        lambda row: get_neon_bboxes(row["product_url"], row["geometry"]), axis=1
    )

    # sjoin to make sure geometries actually overlap with input aoi
    return (
        gf_neon.sjoin(intersects).drop(columns=["index_right"]).reset_index(drop=True)
    )
