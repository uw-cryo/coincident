"""
Searching for NCALM and NOAA data with an input search polygon is made easy via the OpenTopography API
https://portal.opentopography.org/apidocs/
There is also a GitHub repo containing flight extent polygons for NOAA
https://github.com/OpenTopography/Data_Catalog_Spatial_Boundaries/tree/main/NOAA_Coastal_Lidar
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd
import requests  # type: ignore[import-untyped]
from obstore.store import S3Store
from pandas import Timestamp
from shapely.geometry import shape


def _list_subdirs(store: S3Store, prefix: str) -> list[str]:
    """Return immediate subdirectory prefixes under `prefix`, each with a trailing slash."""
    p = prefix if prefix.endswith("/") else prefix + "/"
    result = store.list_with_delimiter(p)
    return [cp.rstrip("/") + "/" for cp in result["common_prefixes"]]


def _find_noaa_dem_files(
    dataset_id: str | int, prefix: str = "dem/"
) -> list[dict[str, str]]:
    """
    HELPER FUNCTION FOR coincident.io.xarray.load_noaa_dem()

    Find all TIFF or LAZ files in the S3 directory containing the specified NOAA dataset ID.
    Supports multiple “geoid” subdirectories under `laz/`.

    Parameters
    ----------
    dataset_id
        NOAA dataset identifier (e.g., "8431", "6260", or "10196").
    prefix : str
        The S3 prefix under which to look. Default is "dem/". Use "laz/" for point clouds.

    Returns
    -------
        A list of dicts, each with keys:
          - 'url': full HTTPS download URL
          - 'name': basename of the object key

    Raises
    ------
    ValueError
        If no directory is found for dataset_id, or if no matching files are found.
    """
    # Initialize unsigned S3 store
    bucket = "noaa-nos-coastal-lidar-pds"
    store = S3Store(bucket, region="us-east-1", skip_signature=True)
    dataset_str = str(dataset_id)

    # 1) If prefix is 'laz/', we must search one level deeper for geoid subfolders
    search_prefixes = (
        _list_subdirs(store, prefix)
        if prefix.rstrip("/").endswith("laz")
        else [prefix if prefix.endswith("/") else prefix + "/"]
    )

    # 2) Find the exact dataset directory under one of the search_prefixes
    # dem/ folders embed the ID as a suffix (e.g. "NGS_South_TampBay_2021_8869/")
    # laz/ folders use the ID directly (e.g. "10149/")
    dataset_dir = None
    for p in search_prefixes:
        for cp in _list_subdirs(store, p):
            if dataset_str in cp.rstrip("/").split("/")[-1]:
                dataset_dir = cp
                break
        if dataset_dir:
            break

    if not dataset_dir:
        msg_no_dir = f"No directory found for dataset ID {dataset_id} under '{prefix}'"
        raise ValueError(msg_no_dir)

    # 3) List all objects under dataset_dir and filter by extension
    f_ext = ".laz" if prefix.rstrip("/").endswith("laz") else ".tif"
    files: list[dict[str, str]] = []
    objects = store.list(dataset_dir).collect()
    for obj in objects:
        key = obj["path"]
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
    dataset: str,
    aoi: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: Timestamp,
    search_end: Timestamp,
    product: str = "lpc",
) -> gpd.GeoDataFrame:
    """
    Perform a search for NCALM LiDAR or NOAA Coastal LiDAR footprints and metadata via
    the OpenTopography API.

    Parameters
    ----------
    dataset
        The dataset type (either "noaa" or "ncalm").
    aoi
        A GeoDataFrame or GeoSeries containing a geometry to restrict the search area, by default does a global search.
    search_start
        The start datetime for the search, by default searches the entire catalog defined by 'dataset'.
    search_end
        The end datetime for the search, by default searches the entire catalog defined by 'dataset'.
    product
        The product type, must be either "lpc" for LiDAR Point Cloud or "dem" for Digital Elevation Model.

    Returns
    -------
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
    dataset: str,
    intersects: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: Timestamp | None = None,
    search_end: Timestamp | None = None,
) -> gpd.GeoDataFrame:
    """
    Integrate the OpenTopography dataset search into the main search function.

    Parameters
    ----------
    dataset
        Dataset alias ('noaa' or 'ncalm').
    intersects
        The geometry to restrict the search.
    search_start
        The start time for the temporal search.
    search_end
        The end time for the temporal search.
    **kwargs
        Additional parameters for the search.

    Returns
    -------
        A GeoDataFrame containing the search results from the OpenTopography API.
    """
    if search_start is None and search_end is None:
        search_start = pd.Timestamp(
            "1996-01-01"
        )  # Starting from January 1, 1996 (first NOAAA dataset)
        search_end = pd.Timestamp.today()  # Default to today's date

    # search using the OpenTopography API
    gf = search_ncalm_noaa(
        dataset=dataset,
        aoi=intersects,
        search_start=search_start,
        search_end=search_end,
    )

    gf["collection"] = dataset.upper()

    return gf
