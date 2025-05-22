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

import re  # for NEON vendor products, client-side spatial filtering
from typing import Any

import geopandas as gpd
import pandas as pd
import requests  # type: ignore[import-untyped]
from pandas import Timestamp
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union


def _get_tile_bbox(file_name: str, tile_size: int = 1000) -> Polygon | None:
    """
    HELPER FUNCTION FOR load_neon_dem()
    Extract tile lower-left coordinates from the file name and return a shapely box.
    Assumes file name pattern like: NEON_D17_TEAK_DP3_317000_4104000_CHM.tif
    - ..._<easting:6digits>_<northing:7digits>_...
    - ..._<easting:7digits>_<northing:7digits>_...
    - ..._<easting:7digits>_<northing:6digits>_...

    Note that we have to do client-side spatial filtering for NEON products since
    their API doesn't allow for this. Assumes that all NEON products are in respective
    UTM zones.
    """
    patterns = [r"(\d{6})_(\d{7})", r"(\d{7})_(\d{7})", r"(\d{7})_(\d{6})"]

    for pattern in patterns:
        match = re.search(pattern, file_name)
        if match:
            easting = int(match.group(1))
            northing = int(match.group(2))
            return box(easting, northing, easting + tile_size, northing + tile_size)

    return None  # return None if no match is found


def _build_neon_point_gf(sites_url: str) -> gpd.GeoDataFrame:
    """
    Build a GeoDataFrame containing NEON site information, including product URLs and available months.
    The geometry is defined as the point (site longitude, site latitude).

    Parameters
    ----------
    sites_url : str
        The URL to the NEON API for retrieving site data.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame with site information and point geometries.
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


def _temporal_filter_neon(
    gf_neon: gpd.GeoDataFrame, search_start: Timestamp, search_end: Timestamp
) -> gpd.GeoDataFrame:
    """
    Filter NEON GeoDataFrame by temporal range using the available months information.

    Parameters
    ----------
    gf_neon : gpd.GeoDataFrame
        The GeoDataFrame containing NEON site data.
    search_start : pd.Timestamp
        The start of the time range.
    search_end : pd.Timestamp
        The end of the time range.

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


def _get_neon_flight_geometry(url: str, fallback_geometry: Any) -> Any:
    """
    Query the NEON API at the given URL and use regex to extract bounding boxes from TIFF file names.
    The tile boxes are unioned to form the full flight geometry. If valid tile boxes are found,
    they are assumed to be in the site's UTM zone. We then reproject them to EPSG:4326 using the
    site's fallback geometry to estimate the proper UTM CRS. If extraction fails, the fallback is returned.

    Parameters
    ----------
    url : str
        The URL to query NEON product data.
    fallback_geometry : geometry
        A fallback geometry (e.g., the site point) if extraction fails.

    Returns
    -------
    geometry
        The unioned flight geometry (in EPSG:4326) or the fallback geometry.
    """
    try:
        response = requests.get(url, timeout=30)
    except Exception:
        return fallback_geometry

    if response.status_code != 200:
        return fallback_geometry

    data = response.json().get("data", {})
    files = data.get("files", [])
    tiles = []
    for f in files:
        # Process only TIFF files that follow the naming convention
        if f["name"].endswith(".tif"):
            tile = _get_tile_bbox(f["name"])
            if tile is not None:
                tiles.append(tile)
    if tiles:
        # Create a union of all tile boxes (in their native UTM units)
        flight_union = unary_union(tiles)
        # Use the fallback (site point) to estimate the correct UTM CRS
        site_gdf = gpd.GeoDataFrame(geometry=[fallback_geometry], crs="EPSG:4326")
        target_crs = site_gdf.estimate_utm_crs()
        # Create a GeoSeries from the union, assign the UTM CRS, then convert to EPSG:4326
        return gpd.GeoSeries([flight_union], crs=target_crs).to_crs("EPSG:4326")[0]
    return fallback_geometry


def search_bboxes(
    intersects: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: Timestamp,
    search_end: Timestamp,
) -> gpd.GeoDataFrame:
    """
    Perform an efficient spatiotemporal search for NEON LiDAR flight geometries by:
      1. Querying the NEON sites.
      2. Temporally filtering based on available months.
      3. Applying a nearest spatial join to restrict to the input AOI.
      4. For each site, fetching the product data and using regex to extract tile
         bounding boxes which are unioned to form the full flight geometry. The resulting
         geometry is reprojected to EPSG:4326.
      5. Finally, ensuring that the resulting flight geometry intersects the AOI.

    Parameters
    ----------
    intersects : gpd.GeoDataFrame | gpd.GeoSeries
        The area of interest (AOI) geometry (in EPSG:4326).
    search_start : pd.Timestamp
        The start date of the temporal filter.
    search_end : pd.Timestamp
        The end date of the temporal filter.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing NEON site metadata with updated flight geometries in EPSG:4326.
    """
    # Set default time range if not provided
    if search_start is None and search_end is None:
        search_start = pd.Timestamp("2013-06-01")  # earliest NEON dataset date
        search_end = pd.Timestamp.today()

    # Define a global AOI if none is provided
    if intersects is None:
        intersects = gpd.GeoSeries([Point(-180, -90).buffer(360)], crs="EPSG:4326")

    sites_url = "http://data.neonscience.org/api/v0/sites"

    # Build and temporally filter the NEON site points
    gf_neon = _build_neon_point_gf(sites_url)
    gf_neon = _temporal_filter_neon(gf_neon, search_start, search_end)

    # Spatially restrict to sites near the AOI using a nearest join
    gf_neon = gpd.sjoin_nearest(
        gf_neon, intersects[["geometry"]], max_distance=1
    ).reset_index(drop=True)
    if "index_right" in gf_neon.columns:
        gf_neon = gf_neon.drop(columns=["index_right"], errors="ignore")

    # Update each site's geometry by fetching the flight geometry via regex-based extraction.
    # The fallback (site point) is also used to estimate the appropriate UTM CRS.
    gf_neon["geometry"] = gf_neon.apply(
        lambda row: _get_neon_flight_geometry(row["product_url"], row["geometry"]),
        axis=1,
    )

    # Final spatial join to ensure that the flight geometries (now in EPSG:4326) intersect the input AOI.
    gf_neon = gf_neon.sjoin(intersects[["geometry"]]).reset_index(drop=True)
    if "index_right" in gf_neon.columns:
        gf_neon = gf_neon.drop(columns=["index_right"], errors="ignore")

    # add days in -dd to the start_datetime and end_datetime which are in yyyy-mm
    gf_neon["start_datetime"] = gf_neon["start_datetime"] + "-01"
    gf_neon["end_datetime"] = (
        gf_neon["end_datetime"]
        .astype(str)
        .apply(lambda m: f"{m}-{pd.Period(m, freq='M').days_in_month:02d}")
    )

    return gf_neon


# NOTE: product code DP1.30003.001 is for discrete LiDAR point cloud
def query_neon_data_api(
    site_id: str, month_str: str, product_code: str = "DP3.30024.001"
) -> dict[str, Any]:
    """
    Query the NEON API for LiDAR products for a given site and month.

    Parameters:
      site_id (str): The NEON site ID.
      month_str (str): The month string in the format YYYY-MM.
      product_code (str): The NEON product code (default is 'DP3.30024.001').

    Returns:
      dict: The JSON response from the NEON API.
    """
    SERVER = "http://data.neonscience.org/api/v0/"
    query_url = f"{SERVER}data/{product_code}/{site_id}/{month_str}"
    response = requests.get(query_url)
    if response.status_code != 200:
        msg_neon_fail = (
            f"NEON API request failed with status code {response.status_code}"
        )
        raise RuntimeError(msg_neon_fail)
    json_data: dict[str, Any] = response.json()  # Explicit cast for MyPy
    return json_data
