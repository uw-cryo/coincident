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


def _find_noaa_dem_files(dataset_id: str | int) -> list[dict[str, str]]:
    """
    HELPER FUNCTION FOR coincident.io.xarray.load_noaa_dem()

    Find all TIFF files in the S3 directory containing the specified NOAA dataset ID.

    Parameters
    ----------
    dataset_id : str or int
        NOAA dataset identifier (e.g., "8431" or "6260").

    Returns
    -------
    list[dict[str, str]]
        A list of dictionaries, each containing 'url' and 'name' keys for each TIFF file.

    Raises
    ------
    ValueError
        If no directory is found for the given dataset ID, or if no TIFF files are found.
    """

    # Initialize the S3 client
    s3_client = boto3.client("s3", region_name="us-east-1")
    s3_client.meta.events.register(
        "choose-signer.s3.*", botocore.handlers.disable_signing
    )

    bucket_name = "noaa-nos-coastal-lidar-pds"
    prefix = "dem/"

    # List all directories under the 'dem/' prefix
    paginator = s3_client.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    dataset_directory = None
    for page in pages:
        for common_prefix in page.get("CommonPrefixes", []):
            folder_name = common_prefix["Prefix"]
            if str(dataset_id) in folder_name:
                dataset_directory = folder_name
                break
        if dataset_directory:
            break

    if not dataset_directory:
        msg_dataset = f"No directory found for dataset ID {dataset_id}"
        raise ValueError(msg_dataset)

    # List all TIFF files in the dataset directory
    tif_files = []
    paginator = s3_client.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket=bucket_name, Prefix=dataset_directory)

    for page in pages:
        for obj in page.get("Contents", []):
            if obj["Key"].endswith(".tif"):
                tif_url = f"https://{bucket_name}.s3.amazonaws.com/{obj['Key']}"
                tif_files.append({"url": tif_url, "name": Path(obj["Key"]).name})

    if not tif_files:
        msg_no_tifs = f"No TIFF files found for dataset ID {dataset_id}"
        raise ValueError(msg_no_tifs)

    return tif_files


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
                else entry["Dataset"]["alternateName"],
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
