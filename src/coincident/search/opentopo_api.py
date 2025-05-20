"""
Searching for NCALM and NOAA data with an input search polygon is made easy via the OpenTopography API
https://portal.opentopography.org/apidocs/
There is also a GitHub repo containing flight extent polygons for NOAA
https://github.com/OpenTopography/Data_Catalog_Spatial_Boundaries/tree/main/NOAA_Coastal_Lidar
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import boto3
import botocore
import geopandas as gpd
import pandas as pd
import requests  # type: ignore[import-untyped]
from pandas import Timestamp
from shapely.geometry import shape


def _find_noaa_dem_files(
    dataset_id: str | int, prefix: str = "dem/"
) -> list[dict[str, str]]:
    """
    HELPER FUNCTION FOR coincident.io.xarray.load_noaa_dem()

    Find all TIFF or LAZ files in the S3 directory containing the specified NOAA dataset ID.
    Supports multiple “geoid” subdirectories under `laz/`.

    Parameters
    ----------
    dataset_id : str or int
        NOAA dataset identifier (e.g., "8431", "6260", or "10196").
    prefix : str
        The S3 prefix under which to look. Default is "dem/". Use "laz/" for point clouds.

    Returns
    -------
    List[Dict[str, str]]
        A list of dicts, each with keys:
          - 'url': full HTTPS download URL
          - 'name': basename of the object key

    Raises
    ------
    ValueError
        If no directory is found for dataset_id, or if no matching files are found.
    """
    # Initialize unsigned S3 client
    s3 = boto3.client("s3", region_name="us-east-1")
    s3.meta.events.register("choose-signer.s3.*", botocore.handlers.disable_signing)
    bucket = "noaa-nos-coastal-lidar-pds"
    dataset_str = str(dataset_id)

    # 1) If prefix is 'laz/', we must search one level deeper for geoid subfolders
    search_prefixes = [prefix]
    if prefix.rstrip("/").endswith("laz"):
        # list possible geoid*/ top level under 'laz/'
        resp = s3.list_objects_v2(Bucket=bucket, Prefix=prefix, Delimiter="/")
        geoids = [cp["Prefix"] for cp in resp.get("CommonPrefixes", [])]
        search_prefixes = geoids

    # 2) Find the exact dataset directory under one of the search_prefixes
    dataset_dir = None
    for p in search_prefixes:
        pages = s3.get_paginator("list_objects_v2").paginate(
            Bucket=bucket,
            Prefix=p,
            Delimiter="/",
        )
        for page in pages:
            for cp in page.get("CommonPrefixes", []):
                if dataset_str + "/" in cp["Prefix"]:
                    dataset_dir = cp["Prefix"]
                    break
            if dataset_dir:
                break
        if dataset_dir:
            break

    if not dataset_dir:
        msg_no_dir = f"No directory found for dataset ID {dataset_id} under '{prefix}'"
        raise ValueError(msg_no_dir)

    # 3) List all objects under dataset_dir and filter by extension
    f_ext = ".laz" if prefix.rstrip("/").endswith("laz") else ".tif"
    files: list[dict[str, str]] = []
    pages = s3.get_paginator("list_objects_v2").paginate(
        Bucket=bucket,
        Prefix=dataset_dir,
    )
    for page in pages:
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.lower().endswith(f_ext):
                files.append(
                    {
                        "url": f"https://{bucket}.s3.amazonaws.com/{key}",
                        "name": Path(key).name,
                    }
                )

    if not files:
        msg_no_files = f"No '{f_ext}' files found for dataset ID {dataset_id}"
        raise ValueError(msg_no_files)

    return files


def search_ncalm_noaa(
    aoi: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: Timestamp,
    search_end: Timestamp,
    product: str = "lpc",
    dataset: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Perform a search for NCALM LiDAR or NOAA Coastal LiDAR footprints and metadata via
    the OpenTopography API.

    Parameters
    ----------
    aoi : gpd.GeoDataFrame | gpd.GeoSeries
        A GeoDataFrame or GeoSeries containing a geometry to restrict the search area, by default does a global search.
    search_start : Timestamp, optional
        The start datetime for the search, by default searches the entire catalog defined by 'dataset'.
    search_end : Timestamp, optional
        The end datetime for the search, by default searches the entire catalog defined by 'dataset'.
    product : str
        The product type, must be either "lpc" for LiDAR Point Cloud or "dem" for Digital Elevation Model.
        Defaults to "lpc".
    dataset : str
        The dataset type (either "noaa" or "ncalm").

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the search results.

    Raises
    ------
    ValueError
        If no valid `aoi` is provided or the API request fails.
    """
    if dataset not in ["noaa", "ncalm"]:
        msg = f"Unsupported dataset: {dataset}"
        raise ValueError(msg)
    if product not in ["lpc", "dem"]:
        msg_product_arg = "Unsupported 'product' argument"
        raise ValueError(msg_product_arg)
    # Convert product shorthand to full description
    if product == "lpc":
        product = "PointCloud"
    elif product == "dem":
        product = "Raster"
    # If `aoi` is None, set the entire world as the search area
    if aoi is None:
        globe = shape(
            {
                "type": "Polygon",
                "coordinates": [
                    [
                        [-180, -90],
                        [180, -90],
                        [180, 90],
                        [-180, 90],
                        [-180, -90],
                    ]
                ],
            }
        )
        search_poly = gpd.GeoSeries(globe, crs="EPSG:4326")
        coords = ",".join(
            [f"{x},{y}" for x, y in search_poly.geometry[0].exterior.coords]
        )
    else:
        # convex_hull works better than simplify for more-complex geometries (ie. Louisiana)
        # https://raw.githubusercontent.com/unitedstates/districts/refs/heads/gh-pages/states/LA/shape.geojson
        search_poly = aoi.to_crs(4326).union_all()
        search_poly_chull = search_poly.convex_hull
        coords = ",".join([f"{x},{y}" for x, y in search_poly_chull.exterior.coords])

    # https://requests.readthedocs.io/en/latest/user/quickstart/#passing-parameters-in-urls
    url_api_base = "https://portal.opentopography.org/API/otCatalog"
    params_api = {
        "productFormat": "PointCloud",
        "detail": "true",
        "outputFormat": "json",
        "polygon": coords,
        "include_federated": "true" if dataset == "noaa" else "false",
    }

    response = requests.get(url_api_base, params=params_api, timeout=30)
    if response.status_code != 200:
        msg_response = f"Error querying OpenTopography API: {response.status_code}"
        raise ValueError(msg_response)

    catalog = response.json()
    if dataset == "noaa":
        catalog = [
            ds
            for ds in catalog["Datasets"]
            if "NOAA" in ds["Dataset"]["identifier"]["propertyID"]
            and ds["Dataset"]["temporalCoverage"] != ""
        ]
    else:
        catalog = [
            ds for ds in catalog["Datasets"] if ds["Dataset"]["temporalCoverage"] != ""
        ]

    # temporal filtering
    start_dt, end_dt = search_start.to_pydatetime(), search_end.to_pydatetime()
    filtered_catalog = [
        entry
        for entry in catalog
        if any(
            start_dt <= datetime.fromisoformat(date) <= end_dt
            for date in (entry["Dataset"]["temporalCoverage"].split(" / "))
        )
    ]

    # convert the filtered catalog into a GeoDataFrame
    gf_lidar = gpd.GeoDataFrame(
        [
            {
                "id": entry["Dataset"]["identifier"]["value"],
                "name": entry["Dataset"]["identifier"]["value"]
                if dataset == "noaa"
                else entry["Dataset"].get("alternateName"),
                "title": entry["Dataset"]["name"],
                "start_datetime": datetime.fromisoformat(
                    entry["Dataset"]["temporalCoverage"].split(" / ")[0]
                ),
                "end_datetime": datetime.fromisoformat(
                    entry["Dataset"]["temporalCoverage"].split(" / ")[-1]
                ),
                "geometry": shape(
                    entry["Dataset"]["spatialCoverage"]["geo"]["geojson"]["features"][
                        0
                    ]["geometry"]
                ),
            }
            for entry in filtered_catalog
        ],
        crs="EPSG:4326",
    )

    # filter the GeoDataFrame to only include datasets that intersect with the search polygon
    return gf_lidar[gf_lidar.intersects(search_poly)].reset_index(drop=True)


# TODO: add meaningful error messages when search fails
def search_opentopo(
    dataset: str | None,
    intersects: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: Timestamp | None = None,
    search_end: Timestamp | None = None,
) -> gpd.GeoDataFrame:
    """
    Integrate the OpenTopography dataset search into the main search function.

    Parameters
    ----------
    dataset : str
        Dataset alias ('noaa' or 'ncalm').
    intersects : gpd.GeoDataFrame | gpd.GeoSeries
        The geometry to restrict the search.
    search_start : Timestamp, optional
        The start time for the temporal search, by default None.
    search_end : Timestamp, optional
        The end time for the temporal search, by default None.
    **kwargs : Any
        Additional parameters for the search.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the search results from the OpenTopography API.
    """
    if search_start is None and search_end is None:
        search_start = pd.Timestamp(
            "1996-01-01"
        )  # Starting from January 1, 1996 (first NOAAA dataset)
        search_end = pd.Timestamp.today()  # Default to today's date

    # search using the OpenTopography API
    return search_ncalm_noaa(
        intersects, search_start=search_start, search_end=search_end, dataset=dataset
    )
