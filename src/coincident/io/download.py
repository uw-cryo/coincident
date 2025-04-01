"""
Convenience functions for downloading STAC items via stac-asset
"""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import geopandas as gpd
import pdal
import requests  # type: ignore[import-untyped]
import stac_asset
from pyproj import CRS
from pystac import Item
from shapely.geometry import box
from tqdm import tqdm

from coincident.io.xarray import _aoi_to_polygon_string, _filter_items_by_project
from coincident.search.neon_api import (
    _get_tile_bbox,
    _utm_crs_from_lonlat,
    query_neon_data_api,
)
from coincident.search.wesm import query_tnm_api


async def download_item(
    item: Item,
    path: str = "/tmp",
    config: str | None = None,
) -> Item:
    """
    Downloads a STAC item to a specified local path.

    Parameters
    ----------
    item : pystac.Item
        The STAC item to be downloaded.
    path : str, optional
        The local directory path where the item will be downloaded. Default is "/tmp".
    config : str, optional
        If config=='maxar', a MAXAR-API-KEY HTTP Header is used for authentication.

    Returns
    -------
    pystac.Item
        The downloaded STAC item.

    Examples
    --------
    >>> localitem = asyncio.run(download_item(item, config=MAXAR_CONFIG))
    """
    posixpath = Path(path)

    if config == "maxar":
        config = stac_asset.Config(
            http_headers={"MAXAR-API-KEY": os.environ.get("MAXAR_API_KEY")}
        )
    # NOTE: prevent in-place modification of remote hrefs with item.clone()
    await stac_asset.download_item(item.clone(), posixpath, config=config)


# Helper function to download files given a list of file dicts.
def _download_files(
    files: list[dict[str, Any]],
    path: str,
    url_key: str,
    product: str,
    chunk_size: int = 65536,  # preset chunk size (64KB)
) -> None:
    """
    Download files given a list of file dictionaries.
    Helper func for download_usgs_dem and download_neon_dem

    Parameters
    ----------
    files : list
        A list of dictionaries where each dict must contain the URL under the key specified by url_key.
    path : str
        The directory path where files will be saved.
    url_key : str
        The key in each file dict that holds the download URL.
    product : str, optional
        Product name used for customizing the progress description.
    chunk_size : int, optional
        The chunk size to use when streaming downloads (default is 64KB).

    Returns
    -------
    None
        Files are saved to the specified path.
    """
    # Set a general description for the outer progress bar based on product if provided.
    outer_desc = (
        f"Downloading {product.upper()} tiles..." if product else "Downloading tiles..."
    )

    for file in tqdm(files, desc=outer_desc, position=0):
        url = file[url_key]
        filename = Path(url).name
        output_file = Path(path) / filename

        # Download the file with streaming and progress tracking.
        with requests.get(url, stream=True) as response:
            response.raise_for_status()
            total_size = int(response.headers.get("content-length", 0))
            inner_progress = tqdm(
                total=total_size,
                unit="iB",
                unit_scale=True,
                desc=filename,
                position=1,
                leave=False,
            )
            with output_file.open("wb") as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        inner_progress.update(len(chunk))
            inner_progress.close()


# NOTE: i do not believe the TNM API supports LPC access for 3DEP :(
# I think that it supports LPD for IfSAR though (tnmdataset = "Lidar Point Cloud (LPC)")
# Note that tnm_api_items has a field 'downloadLazURL' but this is always None from my tests
# https://tnmaccess.nationalmap.gov/api/v1/datasets?
def download_usgs_dem(
    aoi: gpd.GeoDataFrame,
    project: str,
    path: str = "/tmp",
    save_parquet: bool = True,
    tnmdataset: str = "Digital Elevation Model (DEM) 1 meter",
) -> None:
    """
    Download USGS 1-meter DEM tiles based on an AOI by querying the TNM API.
    https://data.usgs.gov/datacatalog/data/USGS:77ae0551-c61e-4979-aedd-d797abdcde0e

    Steps:
      1. Reproject the AOI to EPSG:4326.
      2. Extract the first geometry from the exploded AOI.
      3. Convert the geometry to a polygon string.
      4. Query the TNM API and filter results.
      5. Optionally save metadata as geoparquet.
      6. Download each tile with progress tracking.

    Parameters
    ----------
    aoi : gpd.GeoDataFrame
        Area of interest geometry.
    project : str
        Project identifier to filter results.
    path : str
        Directory path where tiles will be downloaded.
    save_parquet : bool, optional
        Whether to save metadata as geoparquet (default True).
    tnmdataset : str, optional
        TNM dataset identifier (default "Digital Elevation Model (DEM) 1 meter").

    Returns
    -------
    None
        Files are downloaded to the specified output path.
    """
    # 1: Reproject AOI to EPSG:4326 for the TNM API query
    aoi_in_4326 = aoi.to_crs(epsg=4326)

    # 2: Explode and extract the first geometry
    shapely_geom = aoi_in_4326.geometry.explode().iloc[0]

    # 3: Convert geometry to polygon string
    polygon_string = _aoi_to_polygon_string(shapely_geom)

    # 4: Query TNM API and filter results
    tnm_api_items = query_tnm_api(polygon_string, tnmdataset=tnmdataset)
    if not tnm_api_items:
        msg_no_usgs_dem = "TNM API returned no products for the given AOI and dataset."
        raise ValueError(msg_no_usgs_dem)

    filtered_api_items = _filter_items_by_project(tnm_api_items, project)
    if not filtered_api_items:
        msg_no_dem_overlap = (
            "No products found for input project string intersecting the AOI."
        )
        raise ValueError(msg_no_dem_overlap)

    # 5: Optionally save metadata as geoparquet
    if save_parquet:
        geometries = [
            box(
                item["boundingBox"]["minX"],
                item["boundingBox"]["minY"],
                item["boundingBox"]["maxX"],
                item["boundingBox"]["maxY"],
            )
            for item in filtered_api_items
        ]

        gdf = gpd.GeoDataFrame(filtered_api_items, geometry=geometries, crs="EPSG:4326")
        parquet_path = Path(path) / f"stac_tnm_{project}.parquet"
        gdf.to_parquet(parquet_path)

    # 6: Ensure output directory exists
    Path(path).mkdir(parents=True, exist_ok=True)

    # Use the helper function to download files.
    # the 'product' argument here is just for the TQDM progress bar statement
    _download_files(filtered_api_items, path=path, product="DEM", url_key="downloadURL")


def download_usgs_lpc(
    aoi: gpd.GeoDataFrame,
    workunit: str | None = None,
    output_file: str = "/tmp/subset.laz",
    resolution_value: int = 1,
    https_requests: int = 30,
) -> None:
    """
    Download USGS 3DEP LiDAR point cloud data for the given AOI using PDAL's EPT reader.

    Parameters
    ----------
    aoi : gpd.GeoDataFrame
        Area of interest geometry. Must be in EPSG:4326.
    workunit : str, optional
        The 3DEP workunit identifier (e.g., 'USGS_LPC_CO_SoPlatteRiver_Lot5_2013_LAS_2015').
        If not provided, the function automatically finds a workunit that intersects the AOI.
    output_file : str, optional
        Path to save the output LAZ file. Defaults to "/tmp/subset.laz".
    resolution_value : int, optional
        Resolution value for the EPT reader, controlling the level of detail.
        Higher values (e.g., 2, 3) decrease resolution but speed up downloads.
        Defaults to 1 (highest resolution).
    https_requests : int, optional
        Number of concurrent HTTP requests to the EPT server.
        Higher values may increase download speed but could be throttled.
        Defaults to 30.

    Returns
    -------
    None
        The function saves the point cloud data to the specified output file
        but doesn't return any values.

    Notes
    -----
    - The AOI is simplified with a 10-meter tolerance to improve speed.
    - The point cloud is automatically reprojected to match the AOI's coordinate system.
    - For larger areas, consider increasing the resolution_value parameter to speed up downloads.
    """
    # If no workunit is provided, load the resource GeoJSON and find one that intersects the AOI
    # NOTE: I have not tested this thoroughly
    if not workunit:
        resources_url = "https://raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson"
        resources = gpd.read_file(resources_url).set_crs("EPSG:4326")
        aoi_union = aoi.union_all()
        workunit = None
        for _, row in resources.iterrows():
            if row.geometry.intersects(aoi_union):
                workunit = row["name"]
                break
        if not workunit:
            msg_no_workunit_overlap = (
                "No matching 3DEP workunit found for the given AOI."
            )
            raise ValueError(msg_no_workunit_overlap)

    # Construct the EPT URL based on the workunit
    ept_url = (
        f"https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{workunit}/ept.json"
    )

    # Retrieve the dataset metadata to extract the spatial reference.
    response = requests.get(ept_url)
    response.raise_for_status()
    data = response.json()
    srs_wkt = data.get("srs", {}).get("wkt")
    if not srs_wkt:
        msg_ept_wkt_srs = "Could not retrieve SRS information from the EPT endpoint."
        raise ValueError(msg_ept_wkt_srs)
    pointcloud_input_crs = CRS.from_wkt(srs_wkt)

    # Transform the AOI into the point cloud's CRS.
    # Note: CRS.to_epsg() might return None if not EPSG-coded; adjust as needed.
    target_epsg = pointcloud_input_crs.to_epsg() or 3857  # fallback if undefined
    aoi_input = aoi.to_crs(target_epsg)

    # Prepare AOI geometry strings
    aoi_union = aoi_input.union_all().simplify(
        10
    )  # simplify (meters) to remove extra detail for WKT crop speed
    aoi_wkt = aoi_union.wkt

    # Derive bounds string from AOI envelope (format: "([minx, maxx], [miny, maxy])")
    minx, miny, maxx, maxy = aoi_input.total_bounds
    bounds_str = f"([{minx}, {maxx}], [{miny}, {maxy}])"

    # Build the PDAL pipeline to read from the EPT endpoint, crop to AOI, and write to a LAZ file.
    pipeline_json = {
        "pipeline": [
            {
                "type": "readers.ept",
                "filename": ept_url,
                "bounds": bounds_str,
                "requests": https_requests,  # concurrent http requests
                "resolution": resolution_value,  # optional: get coarser level
                # Removed 'threads' parameter - it's causing conflicts
            },
            {"type": "filters.crop", "polygon": aoi_wkt},
            {"type": "writers.las", "filename": output_file, "compression": "laszip"},
        ]
    }

    pipeline = pdal.Pipeline(json.dumps(pipeline_json))
    pipeline.execute_streaming()


def download_neon_dem(
    aoi: gpd.GeoDataFrame,
    datetime_str: str,
    site_id: str,
    product: str,
    path: str = "/tmp",
) -> None:
    """
    Download NEON LiDAR tiles (DSM, DTM, CHM, or LPC) based on an AOI by querying the NEON API.

    Steps:
      1. Convert the datetime string to a month string in the format YYYY-MM.
      2. Determine appropriate UTM CRS for the AOI for spatial filtering.
      3. Query the NEON API using site ID and month.
      4. Filter the returned files based on product type and spatial intersection.
      5. Download filtered tiles with progress tracking.

    Parameters:
      aoi (gpd.GeoDataFrame): Area of interest geometry to query against.
      datetime_str (str): Date string in YYYY-MM-DD format.
      site_id (str): NEON site identifier.
      product (str): Product type to download ('dsm', 'dtm', 'chm', or 'lpc').
      path (str): Directory path where tiles will be downloaded (default "/tmp").

    Returns:
      returns None after downloading files
    """
    if product not in ["dsm", "dtm", "chm", "lpc"]:
        msg_invalid_prod = (
            "Invalid product type. Choose from 'dsm', 'dtm', 'chm', or 'lpc'."
        )
        raise ValueError(msg_invalid_prod)
    # 1: Convert datetime_str to month string
    # TODO: make this simpler
    dt = datetime.strptime(datetime_str, "%Y-%m-%d")
    month_str = dt.strftime("%Y-%m")

    # 2: Determine appropriate UTM CRS for the AOI
    centroid = aoi.union_all().centroid
    target_crs = _utm_crs_from_lonlat(centroid.x, centroid.y)
    aoi_utm = aoi.to_crs(target_crs)
    aoi_geom_utm = aoi_utm.union_all()  # Combined AOI in UTM coordinates

    # 3: Query the NEON API and based on the product type
    product_code = "DP3.30024.001" if product != "lpc" else "DP1.30003.001"
    data_json = query_neon_data_api(site_id, month_str, product_code=product_code)

    # 4: Filter files based on product type and spatial intersection
    files = data_json["data"]["files"]
    filtered_files = []
    product_filter = (
        f"{product.upper()}.tif"
        if product != "lpc"
        else "classified_point_cloud_colorized.laz"
    )
    for f in files:
        if product_filter in f["name"]:
            tile_bbox = _get_tile_bbox(f["name"])
            if tile_bbox is not None and tile_bbox.intersects(aoi_geom_utm):
                filtered_files.append(f)
    if not filtered_files:
        msg_neon_empty = "No matching LiDAR files found for the given spatial query."
        raise RuntimeError(msg_neon_empty)

    # 5: Use the helper function to download files
    _download_files(filtered_files, path=path, url_key="url", product=product)
