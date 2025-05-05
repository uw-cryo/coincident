"""
Convenience functions for downloading STAC items via stac-asset
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import Any

import boto3
import geopandas as gpd
import rasterio
import requests  # type: ignore[import-untyped]
import rioxarray as rxr
import stac_asset
from botocore import UNSIGNED
from botocore.client import Config
from pystac import Item
from shapely.geometry import Polygon, box
from tqdm import tqdm

from coincident.io.xarray import _aoi_to_polygon_string, _filter_items_by_project
from coincident.search.neon_api import (
    _get_tile_bbox,
    _utm_crs_from_lonlat,
    query_neon_data_api,
)
from coincident.search.opentopo_api import _find_noaa_dem_files
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
    output_dir: str = "/tmp",
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
        parquet_path = Path(output_dir) / f"stac_tnm_{project}.parquet"
        gdf.to_parquet(parquet_path)

    # 6: Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Use the helper function to download files.
    # the 'product' argument here is just for the TQDM progress bar statement
    _download_files(
        filtered_api_items, path=output_dir, product="DEM", url_key="downloadURL"
    )


def download_neon_dem(
    aoi: gpd.GeoDataFrame,
    datetime_str: str,
    site_id: str,
    product: str,
    output_dir: str = "/tmp",
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
      output_dir (str): Directory path where tiles will be downloaded (default "/tmp").

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
    _download_files(filtered_files, path=output_dir, url_key="url", product=product)


def download_noaa_dem(
    aoi: gpd.GeoDataFrame,
    dataset_id: str | int,
    output_dir: str = "/tmp",
    res: int = 1,
    clip: bool = True,
) -> None:
    """
    Download NOAA coastal lidar DEM data from S3 based on a dataset identifier and an AOI.
    For each tile that intersects the AOI, the tile is downloaded, optionally coarsened and clipped,
    then saved as a GeoTIFF file in the specified output directory.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest geometry (assumed to be in EPSG:4326 or any other CRS that will be reprojected).
    dataset_id : str or int
        NOAA dataset identifier (e.g., "8431" or "6260").
    output_dir : str, optional
        Directory to save downloaded output files. Defaults to "/tmp".
    res : int, optional
        Resolution factor. If greater than 1, the DEMs will be coarsened by this factor. Defaults to 1.
    clip : bool, optional
        Whether to clip each DEM tile to the AOI. Defaults to True.

    Returns
    -------
    None
        Processed DEM tiles are saved as GeoTIFF files to the specified output directory.

    Notes
    -----
    The function retrieves only the header information from each file on S3 to determine tile bounds,
    then filters and processes only those tiles intersecting the AOI. Files are processed concurrently
    to reduce overall network latency.
    """
    # Ensure the output directory exists.
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Step 1: List all .tif files within the NOAA dataset directory using a helper function.
    tif_files = _find_noaa_dem_files(dataset_id)
    tile_urls = [file_info["url"] for file_info in tif_files]

    # Helper function to fetch a tile's bounds and retain its URL.
    def fetch_tile_geometry(
        url: str,
    ) -> tuple[str, Polygon, rasterio.crs.CRS]:
        # Open the raster header from S3 using rasterio.Env for proper configuration.
        with rasterio.Env(), rasterio.open(url) as src:
            bounds = src.bounds
            tile_crs = src.crs
        # Create a shapely bounding box polygon from the raster bounds.
        polygon = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
        return url, polygon, tile_crs

    # Step 2: Retrieve header metadata for each tile concurrently.
    with ThreadPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(fetch_tile_geometry, tile_urls))

    # Unpack the results into separate lists (url, geometry, and CRS).
    urls, tile_geometries, crs_list = zip(*results, strict=False)
    common_crs = crs_list[0]  # Assume all tiles share the same CRS.

    # Create a GeoDataFrame to map tile geometries to their source URLs.
    tiles_gdf = gpd.GeoDataFrame(
        {"url": urls, "geometry": tile_geometries}, crs=common_crs
    )

    # Step 3: Reproject the AOI to the common CRS of the raster tiles.
    aoi = aoi.to_crs(common_crs)

    # Step 4: Filter tiles by determining which ones intersect the AOI.
    aoi_union = aoi.union_all()
    selected_tiles = tiles_gdf[tiles_gdf.intersects(aoi_union)]

    if selected_tiles.empty:
        msg = f"No TIFF files intersect with the provided AOI for dataset ID {dataset_id}."
        raise ValueError(msg)

    # Step 5: Process and download each intersecting DEM tile.
    for _, row in selected_tiles.iterrows():
        tile_url = row["url"]
        # Extract the filename from the URL.
        file_name = tile_url.split("/")[-1]
        output_file = (
            output_path / f"clipped_{file_name}" if clip else output_path / file_name
        )

        try:
            # Open the remote DEM using rioxarray (reads only metadata/chunks, not full data into memory).
            dem = rxr.open_rasterio(tile_url, masked=True, chunks="auto")
            # Optionally coarsen the DEM if a resolution factor greater than 1 is provided.
            if res > 1:
                dem = dem.coarsen(x=res, y=res, boundary="trim").mean()
            # If clipping is enabled, reproject the AOI to the DEM's CRS and clip.
            if clip:
                dem_crs = dem.rio.crs
                aoi_tile_crs = aoi.to_crs(dem_crs)
                dem = dem.rio.clip(aoi_tile_crs.geometry, dem_crs, drop=True)
            # Save the processed DEM as a GeoTIFF.
            dem.rio.to_raster(output_file)
        except Exception as error:
            msg = f"Unable to process or save tile from URL {tile_url}. Error: {error}"
            raise RuntimeError(msg) from error


def download_ncalm_dem(
    aoi: gpd.GeoDataFrame,
    dataset_id: str | int,
    output_dir: str = "/tmp",
) -> None:
    """
    Download NCALM (or OpenTopo user-hosted) data from OpenTopography S3 based on a dataset identifier and an AOI.
    Your dataset identifier will be the 'name' column from coincident.search.search(dataset="ncalm")
    NOTE: this function will save individual tiles to the output directory of interest, this is due
    to the messier structure of the NCALM OpenTopo catalog (especially confusing for DEMs)

    For DEM  data, the function lists all available DEM tiles under the
    dataset prefix, downloads each tile, reprojects it (if necessary) to the AOI's local UTM
    CRS, clips it to the AOI, and then saves the clipped DEM as a GeoTIFF file.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest geometry (assumed to be in EPSG:4326).
    dataset_id : str or int
        Dataset identifier (e.g., "WA18_Wall" or "OTLAS.072019.6339.1"). This is used as the
        prefix for S3 object keys.
    output_dir : str, optional
        Directory to save output files. Defaults to "/tmp".

    Returns
    -------
    None
        The function writes the cropped/clipped output files to the specified output directory.
    """
    bucket = "raster"  # DEM (raster) bucket
    file_ext = ".tif"  # expecting GeoTIFF files for DEMs

    # Ensure the output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Set the endpoint URL for OpenTopo S3
    endpoint_url = "https://opentopography.s3.sdsc.edu"

    # Create an S3 client configured for unsigned access
    s3 = boto3.client(
        "s3", config=Config(signature_version=UNSIGNED), endpoint_url=endpoint_url
    )

    # List objects under the dataset prefix
    prefix = f"{dataset_id}/"
    response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    if "Contents" not in response:
        msg_empty_contents = f"No objects found for prefix {prefix} in bucket {bucket}."
        raise ValueError(msg_empty_contents)

    # Estimate the AOI's local UTM CRS and reproject the AOI, most NCALM flight are in local UTMs
    # we add another check for this later
    aoi_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(aoi_crs)

    # Process DEM files: collect all keys ending with the DEM extension
    dem_keys = [
        obj["Key"]
        for obj in response.get("Contents", [])
        if obj["Key"].endswith(file_ext)
    ]
    for dem_key in dem_keys:
        dem_url = f"{endpoint_url}/{bucket}/{dem_key}"
        dem = rxr.open_rasterio(dem_url, masked=True)
        if dem.rio.crs != aoi_crs:
            dem = dem.rio.reproject(aoi_crs)
        dem_clipped = dem.rio.clip(aoi_utm.geometry, aoi_utm.crs)
        clipped_output_path = Path(output_dir) / f"clipped_{Path(dem_key).name}"
        dem_clipped.rio.to_raster(clipped_output_path)
