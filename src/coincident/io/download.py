"""
Convenience functions for downloading STAC items via stac-asset
"""

from __future__ import annotations

import json
import logging
import os
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import Any

import boto3
import geopandas as gpd
import requests  # type: ignore[import-untyped]
import rioxarray as rxr
import stac_asset
from botocore import UNSIGNED
from botocore.client import Config
from pyproj import CRS
from pystac import Item
from shapely.geometry import box
from tqdm.autonotebook import tqdm

from coincident.datasets.nasa import GLiHT
from coincident.io.xarray import (
    _aoi_to_polygon_string,
    _fetch_noaa_tile_geometry,
    _filter_items_by_project,
    _translate_gliht_id,
)
from coincident.search import search
from coincident.search.neon_api import (
    _get_tile_bbox,
    query_neon_data_api,
)
from coincident.search.opentopo_api import _find_noaa_dem_files
from coincident.search.wesm import query_tnm_api

logger = logging.getLogger(__name__)


def _set_config(
    item: Item,
    config: dict[str, Any] | None = None,
) -> stac_asset.Config:
    """Set stac_asset config. Automatically add Maxar API key if needed"""
    if item.properties.get("constellation") == "maxar" and config is None:
        config = stac_asset.Config(
            http_headers={"MAXAR-API-KEY": os.environ.get("MAXAR_API_KEY")}
        )
    else:
        config = stac_asset.Config(**config)

    return config


async def read_href(
    item: Item,
    asset: str,
    config: dict[str, Any] | None = None,
) -> Item:
    """
    Open and read a STAC asset href into local memory

    Parameters
    ----------
    item : pystac.Item
        The STAC item to be downloaded.
    asset : str
        The asset name to open (e.g. "cloud-cover" for maxar).
    config : dict, optional
        dictionary of options for :class:`~stac_asset.Config`

    Returns
    -------
    bytes
        Bytes read from file at corresponding href

    Examples
    --------
    Read cloud-cover MultiPolygon estimate from Maxar API
    >>> bytes = asyncio.run(read_href(item, "cloud-cover" config=MAXAR_CONFIG))
    """
    href = item.assets.get(asset).href

    config = _set_config(item, config)

    return await stac_asset.read_href(href, config=config)


async def download_item(
    item: Item,
    path: str = "/tmp",
    config: dict[str, Any] | None = None,
) -> Item:
    """
    Downloads a STAC item to a specified local path.

    Parameters
    ----------
    item : pystac.Item
        The STAC item to be downloaded.
    path : str, optional
        The local directory path where the item will be downloaded. Default is "/tmp".
    config : dict, optional
        dictionary of options for :doc:`stac_asset.Config <stac_asset:stac_asset.Config>`

    Returns
    -------
    pystac.Item
        The downloaded STAC item.

    Examples
    --------
    Download all assets for given item to /tmp directory
    >>> localitem = asyncio.run(download_item(item, config=MAXAR_CONFIG))
    """
    posixpath = Path(path)

    config = _set_config(item, config)

    # NOTE: prevent in-place modification of remote hrefs with item.clone()
    return await stac_asset.download_item(item.clone(), posixpath, config=config)


def download_files(
    files: list[str],
    output_dir: str,
    product: str | None = "",
    chunk_size: int = 65536,  # preset chunk size (64KB)
    earthdata_auth: bool = False,
) -> None:
    """
    Download files from remote URLs (cloud storage, web servers, etc.)
    Used in the download_dem and fetch_lpc_tiles functions

    Parameters
    ----------
    files : list
        A list of file urls to download
    output_dir : str
        The directory path where files will be saved.
    product : str, optional
        Product name used for customizing the progress description.
    chunk_size : int, optional
        The chunk size to use when streaming downloads (default is 64KB).
    earthdata_auth : bool, optional
        Whether to use NASA Earthdata authentication (default False).

    Returns
    -------
    None
        Files are saved to the specified path.
    """
    # Handle Earthdata authentication if needed
    auth = None
    if earthdata_auth:
        username = os.getenv("EARTHDATA_USER")
        password = os.getenv("EARTHDATA_PASS")
        if not username or not password:
            msg_no_ed_cres = """
            EARTHDATA_USER and EARTHDATA_PASS must be set in the environment.
            os.environ["EARTHDATA_USER"] = "username"
            os.environ["EARTHDATA_PASS"] ="password"
            """
            raise OSError(msg_no_ed_cres)
        auth = (username, password)

    outer_desc = (
        f"Downloading {product.upper()} tiles..." if product else "Downloading tiles..."
    )
    for url in tqdm(files, desc=outer_desc, position=0):
        filename = Path(url).name
        dest = Path(output_dir) / filename

        # Skip if file already exists
        if dest.exists():
            continue

        resp = requests.get(url, stream=True, auth=auth)
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        inner = tqdm(
            total=total,
            unit="iB",
            unit_scale=True,
            desc=filename,
            position=1,
            leave=False,
        )
        with dest.open("wb") as f:
            for chunk in resp.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    inner.update(len(chunk))
        inner.close()


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
    download_urls = [item["downloadURL"] for item in filtered_api_items]
    download_files(download_urls, output_dir=output_dir, product="DEM")


def _fetch_usgs_lpc_tiles(
    aoi: gpd.GeoDataFrame,
    project: str,
    output_dir: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Query USGS TNMAccess for LPC products intersecting an AOI, assemble their
    bounding boxes into a GeoDataFrame (EPSG:4326), and optionally write to disk.

    Steps:
      1. Reproject AOI to EPSG:4326 and explode to a single geometry.
      2. Convert that geometry into the TNM APIs `polygon` string format.
      3. Call `query_tnm_api` to retrieve all LPC products (LAS/LAZ) in the AOI.
      4. Filter the returned items by matching `project` via `_filter_items_by_project`.
      5. Build a Shapely box from each items `boundingBox` coordinates.
      6. Assemble into a GeoDataFrame in EPSG:4326.
      7. If `output_dir` is provided, save as `usgs_lpc_products.geojson`.

    Parameters:
      aoi (gpd.GeoDataFrame): Area of interest.
      project (str): A substring or identifier used to filter product titles.
      output_dir (Optional[str]): Directory to write out the GeoJSON. If None,
                                  no file is written.

    Returns:
      gpd.GeoDataFrame: Contains columns from the API plus a `geometry` column
                        with each products bounding box in WGS84.
    """
    # 1. AOI in lon/lat for TNM query
    aoi_4326 = aoi.to_crs(epsg=4326)
    single_geom = aoi_4326.geometry.explode(index_parts=False).iloc[0]

    # 2. Build polygon string for the API
    polygon_string = _aoi_to_polygon_string(single_geom)

    # 3. Query TNM API (only LAS/LAZ for LPC)
    items = query_tnm_api(
        polygon_string, tnmdataset="Lidar Point Cloud (LPC)", prodFormats="LAS,LAZ"
    )

    if not items:
        msg_empty_api = "TNM API returned no products for the given AOI/dataset."
        raise ValueError(msg_empty_api)

    # 4. Filter by project identifier
    filtered = _filter_items_by_project(items, project)
    if not filtered:
        msg_no_lpc = "No LPC products match the project filter and AOI."
        raise ValueError(msg_no_lpc)

    # 5. Build bbox geometries
    geometries = []
    for itm in filtered:
        bb = itm.get("boundingBox", {})
        geom = box(bb["minX"], bb["minY"], bb["maxX"], bb["maxY"])
        geometries.append(geom)

    # 6. Assemble GeoDataFrame in WGS84
    gdf = gpd.GeoDataFrame(filtered, geometry=geometries, crs="EPSG:4326")
    gdf = gdf[["sourceId", "downloadURL", "geometry"]].rename(
        columns={"sourceId": "name", "downloadURL": "url"}
    )

    # 7. Optionally write out
    if output_dir:
        out_path = Path(output_dir) / f"USGS_{project}_lpc_tiles.geojson"
        gdf.to_file(out_path, driver="GeoJSON")

    return gdf


def build_usgs_ept_pipeline(
    aoi: gpd.GeoDataFrame,
    workunit: str | None = None,
    output_dir: str | None = None,
) -> dict[str, Any]:
    """
    Finds the relevant USGS 3DEP EPT dataset for an AOI and returns a dictionary
    structured as a PDAL pipeline JSON to subset the data based on the AOI,
    without downloading point data.

    The returned pipeline includes an 'readers.ept' stage and a 'filters.crop'
    stage configured with the EPT URL, the AOI's bounds, and polygon WKT,
    all in the EPT's spatial reference system (SRS). Users can add their
    custom parameters (like writers and additional filters) to this pipeline
    before executing it with PDAL.

    Parameters:
        aoi (gpd.GeoDataFrame): Area of interest (can be in any CRS).
        workunit (Optional[str]): The specific 3DEP workunit name from coincident.search.search(dataset="3dep")
                                  (e.g., 'CO_CentralEasternPlains_1_2020').
                                  If None, the function will attempt to find a
                                  workunit intersecting the AOI using the USGS
                                  resources GeoJSON.
        output_dir (Optional[str]): Directory to write out the generated PDAL
                                    pipeline JSON file. If None, no file is written.

    Returns:
        dict: A dictionary representing the base PDAL pipeline JSON for subsetting.

    Raises:
        ValueError: If the AOI geometry is not a Polygon/MultiPolygon, or if
                    no matching 3DEP workunit is found for the given AOI.
        requests.exceptions.RequestException: If there's an issue fetching the
                                              EPT metadata.
        RuntimeError: If essential data like SRS WKT is missing/cannot be parsed,
                      or if the resources GeoJSON cannot be read.
    """
    if not aoi.geom_type.isin(["Polygon", "MultiPolygon"]).all():
        aoi_geom_msg = (
            "AOI geometry must be a Polygon or MultiPolygon for PDAL filters.crop."
        )
        raise ValueError(aoi_geom_msg)

    # --- 1. Find the relevant 3DEP workunit for the AOI ---
    if workunit is None:
        resources_url = "https://raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson"
        try:
            resources = gpd.read_file(resources_url).set_crs("EPSG:4326")
        except Exception as e:
            err_msg = f"Could not read USGS resources GeoJSON from {resources_url}: {e}"
            raise RuntimeError(err_msg) from e

        # 4326 for spatial intersect with resources
        aoi_4326 = aoi.to_crs("EPSG:4326")
        aoi_union_4326 = aoi_4326.unary_union

        found_workunit = None
        aoi_bounds_4326 = aoi_4326.total_bounds
        for _, row in resources.iterrows():
            # nested if statement for ruff SIM102
            if (
                row.geometry
                and row.geometry.bounds.intersects(aoi_bounds_4326)
                and row.geometry.intersects(aoi_union_4326)
            ):
                found_workunit = row["name"]
                break  # Found the first intersecting workunit

        if found_workunit is None:
            err_msg = "No matching 3DEP workunit found for the given AOI."
            raise ValueError(err_msg)
        workunit = found_workunit

    # --- 2. Construct and fetch the EPT metadata (ept.json) ---
    ept_url = (
        f"https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{workunit}/ept.json"
    )
    try:
        response = requests.get(ept_url)
        response.raise_for_status()
        ept_data = response.json()
    except requests.exceptions.RequestException as e:
        err_msg = f"Failed to fetch EPT metadata from {ept_url}: {e}"
        raise requests.exceptions.RequestException(err_msg) from e
    except json.JSONDecodeError as e:
        err_msg = f"Failed to decode JSON from EPT metadata at {ept_url}: {e}"
        raise RuntimeError(
            err_msg
        ) from e  # Changed from json.JSONDecodeError to RuntimeError for mypy

    # --- 3. Determine EPT CRS and transform AOI ---
    srs_info = ept_data.get("srs", {})
    srs_wkt = srs_info.get("wkt")

    if srs_wkt is None:
        # Check for horizontal/vertical if wkt is missing, though WKT is preferred
        horizontal_epsg = srs_info.get("horizontal")
        # vertical_epsg = srs_info.get("vertical") # Not needed for horizontal CRS
        if horizontal_epsg:
            try:
                pointcloud_input_crs = CRS.from_epsg(int(horizontal_epsg))
            except Exception as e:
                err_msg = f"Could not create CRS from EPSG:{horizontal_epsg} found in EPT metadata: {e}"
                raise RuntimeError(err_msg) from e
        else:
            err_msg = "Could not retrieve SRS WKT or horizontal EPSG from the EPT endpoint metadata."
            raise RuntimeError(err_msg)
    else:
        try:
            pointcloud_input_crs = CRS.from_wkt(srs_wkt)
        except Exception as e:
            err_msg = f"Could not parse SRS WKT from EPT metadata: {e}"
            raise RuntimeError(err_msg) from e

    aoi_input = aoi.to_crs(pointcloud_input_crs)

    # --- 4. Prepare AOI geometry strings/bounds in EPT CRS for PDAL pipeline ---
    aoi_union = aoi_input.union_all()
    aoi_wkt_geom = aoi_union.simplify(10)  # simplify in meters (UTM zones)
    aoi_wkt = aoi_wkt_geom.wkt

    # Derive bounds from AOI envelope (for readers.ept 'bounds' parameter)
    # PDAL expects bounds in the format "(([minx, maxx], [miny, maxy]))"
    minx, miny, maxx, maxy = aoi_input.total_bounds
    # f-string for UP031
    bounds_str_pdal = f"(([{minx:.10f}, {maxx:.10f}], [{miny:.10f}, {maxy:.10f}]))"

    # --- 5. Assemble the base PDAL pipeline dictionary ---
    # This pipeline defines the input EPT source and the initial spatial subsetting (bounds and crop)
    # Users will add their own writers, filters, etc., to this pipeline.
    pdal_pipeline = {
        "pipeline": [
            {
                "type": "readers.ept",
                "filename": ept_url,
                "bounds": bounds_str_pdal,
                # Optional: Add concurrent requests parameter? This depends on user environment.
                # Let's omit it for a minimal template.
                # "requests": 8,
                # Optional: Add threads? Omit for minimal template.
                # "threads": 2,
            },
            {
                "type": "filters.crop",
                "polygon": aoi_wkt,
                # Optional: Add a name to the filter? Omit for minimal template.
                # "tag": "aoi_crop"
            },
            # User adds their custom stages (filters, writers) here
            {
                "type": "writers.las",
                "filename": f"{workunit}_EPT_subset_pipeline.laz",
                "compression": "laszip",
            },
        ]
    }

    # --- 6. Optionally write out the pipeline dictionary to JSON ---
    if output_dir:
        out_path = Path(output_dir) / f"{workunit}_EPT_subset_pipeline.json"
        # out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w") as f:
            json.dump(pdal_pipeline, f, indent=4)

    # --- 7. Return the pipeline dictionary ---
    return pdal_pipeline


def download_neon_dem(
    aoi: gpd.GeoDataFrame,
    datetime_str: str,
    site_id: str,
    product: str,
    output_dir: str = "/tmp",
) -> None:
    """
    Download NEON LiDAR tiles (DSM, DTM, or CHM) based on an AOI by querying the NEON API.

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
      product (str): Product type to download ('dsm', 'dtm', or 'chm').
      output_dir (str): Directory path where tiles will be downloaded (default "/tmp").

    Returns:
      returns None after downloading files
    """
    if product not in ["dsm", "dtm", "chm"]:
        msg_invalid_prod = "Invalid product type. Choose from 'dsm', 'dtm', 'chm'."
        raise ValueError(msg_invalid_prod)
    # 1: Convert datetime_str to month string
    # TODO: make this simpler
    dt = datetime.strptime(datetime_str, "%Y-%m-%d")
    month_str = dt.strftime("%Y-%m")

    # 2: Determine appropriate UTM CRS for the AOI
    target_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(target_crs)
    aoi_geom_utm = aoi_utm.union_all()  # Combined AOI in UTM coordinates

    # 3: Query the NEON API and based on the product type
    data_json = query_neon_data_api(site_id, month_str, product_code="DP3.30024.001")

    # 4: Filter files based on product type and spatial intersection
    files = data_json["data"]["files"]
    filtered_files = []
    product_filter = f"{product.upper()}.tif"
    for f in files:
        if product_filter in f["name"]:
            tile_bbox = _get_tile_bbox(f["name"])
            if tile_bbox is not None and tile_bbox.intersects(aoi_geom_utm):
                filtered_files.append(f)
    if not filtered_files:
        msg_neon_empty = "No matching LiDAR files found for the given spatial query."
        raise RuntimeError(msg_neon_empty)

    # 5: Use the helper function to download files
    download_urls = [item["url"] for item in filtered_files]
    download_files(download_urls, output_dir=output_dir, product=product)


def _fetch_neon_lpc_tiles(
    aoi: gpd.GeoDataFrame,
    datetime_str: str,
    site_id: str,
    output_dir: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Fetch NEON LPC tile bboxes as a GeoDataFrame (in EPSG:4326),
    optionally writing to disk as a GeoJSON.
    This function uses regex and assumes a tile size of 1km for all NEON tiles
    given the necessity for client-side spatial filtering and no PDAL

    Steps:
      1. Convert the datetime string to YYYY-MM for the API query.
      2. Determine the appropriate UTM CRS for spatial filtering.
      3. Query the NEON API for LPC files.
      4. Filter files whose 1km tile intersects the AOI in UTM.
      5. Parse each filename for its lower left UTM corner, build a 1km box.
      6. Assemble into a GeoDataFrame and reproject to EPSG:4326.
      7. If `output_dir` is provided, save the GeoDataFrame as `lpc_tiles.geojson`.

    Parameters:
      aoi (gpd.GeoDataFrame): AOI geometries.
      datetime_str (str): Date in “YYYY-MM-DD” format.
      site_id (str): NEON site identifier.
      output_dir (Optional[str]): Directory to write `lpc_tiles.geojson`.
                                  If None, nothing is written.

    Returns:
      gpd.GeoDataFrame: Tile extent polygons in lon/lat.
    """
    # 1. Month string
    dt = datetime.strptime(datetime_str, "%Y-%m-%d")
    month_str = dt.strftime("%Y-%m")

    # 2. AOI in UTM for spatial filtering
    utm_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(utm_crs)
    aoi_u = aoi_utm.union_all()

    # 3. Query NEON API for LPC product
    data = query_neon_data_api(site_id, month_str, product_code="DP1.30003.001")
    files = data["data"]["files"]

    # 4. Filter by filename and AOI intersection
    pattern = re.compile(r"_DP1_(\d+)_(\d+)_classified_point_cloud_colorized\.laz$")
    records = []
    for f in files:
        if "classified_point_cloud_colorized.laz" not in f["name"]:
            continue

        m = pattern.search(f["name"])
        if not m:
            # skip any unexpected naming
            continue

        easting, northing = map(float, m.groups())
        # assumes 1km x 1km tiles
        tile = box(easting, northing, easting + 1000, northing + 1000)

        if not tile.intersects(aoi_u):
            continue

        records.append(
            {
                "name": f["name"],
                "url": f["url"],
                "geometry": tile,
            }
        )

    if not records:
        msg_no_lpc_tiles = "No LPC tiles intersect your AOI for the given date."
        raise RuntimeError(msg_no_lpc_tiles)

    # 5. Build GeoDataFrame & reproject
    gdf = gpd.GeoDataFrame(records, crs=utm_crs)
    gdf = gdf.to_crs(epsg=4326)

    # 6. Optionally write out
    if output_dir:
        out_path = Path(output_dir) / f"NEON_{site_id}_{datetime_str}_lpc_tiles.geojson"
        gdf.to_file(out_path, driver="GeoJSON")

    return gdf


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
    output_path = Path(output_dir)
    # output_path.mkdir(parents=True, exist_ok=True)

    # Step 1: List all .tif files within the NOAA dataset directory using a helper function.
    tif_files = _find_noaa_dem_files(dataset_id)
    tile_urls = [file_info["url"] for file_info in tif_files]

    # Step 2: Retrieve header metadata for each tile concurrently.
    with ThreadPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(_fetch_noaa_tile_geometry, tile_urls))

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


def _fetch_noaa_lpc_tiles(
    aoi: gpd.GeoDataFrame,
    dataset_id: str,
    output_dir: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Fetch NOAA DEM tile geometries as a GeoDataFrame (in EPSG:4326),
    optionally writing to disk as a GeoJSON.

    This function processes NOAA DEM LAZ files, extracts tile geometries based on filename patterns,
    filters them by intersection with the provided AOI, and returns a GeoDataFrame in EPSG:4326.

    NOTE: regex was necessary here for client-side spatial filtering, otherwise you'd have to make network
    requests to thousands of laz files

    Steps:
      1. Estimate the appropriate UTM CRS for spatial filtering.
      2. Query the NOAA dataset for LAZ files.
      3. Extract easting and northing from filenames using regex.
      4. Determine tile dimensions by analyzing neighboring tiles.
      5. Construct tile geometries and filter by AOI intersection.
      6. Assemble into a GeoDataFrame and reproject to EPSG:4326.
      7. If output_dir is provided, save the GeoDataFrame as noaa_dem_tiles.geojson.

    Parameters:
      aoi (gpd.GeoDataFrame): AOI geometries.
      dataset_id (str): NOAA dataset identifier.
      output_dir (Optional[str]): Directory to write noaa_dem_tiles.geojson.
                                  If None, nothing is written.

    Returns:
      gpd.GeoDataFrame: Tile extent polygons in lon/lat.
    """
    # 1. Estimate UTM CRS for AOI
    utm_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(utm_crs)
    aoi_union = aoi_utm.union_all()

    # 2. Query NOAA dataset for LAZ files
    laz_files = _find_noaa_dem_files(dataset_id, prefix="laz/")

    # 3. Define patterns to match filenames (NOAA has a lot of different formats)
    patterns = [
        re.compile(
            r"(\d+)_(\d+)e_(\d+)"
        ),  # Matches patterns like '20221027_252000e_3268500n'
        re.compile(
            r"(\d{4})_(\d{4})_(\d{2})_(\d{2})"
        ),  # Matches patterns like '2022_1027_25_15'
        re.compile(
            r"_\w+_(\d+)_(\d+)_\d+\."
        ),  # Matches patterns like '_abc_1234_5678_90.'
        re.compile(
            r"(\w+)_lpc_\d+ft_\d+$"
        ),  # Matches patterns like 'abc_lpc_10ft_1234'
        re.compile(r"(\w+)_lpc$"),  # Matches patterns like 'abc_lpc'
        re.compile(r"_(\d+)e_(\d+)n"),  # Matches patterns like '_252000e_3268500n'
    ]

    # 4. Extract easting and northing from filenames
    tile_coords = []
    row_tiles = defaultdict(list)
    col_tiles = defaultdict(list)

    for f in laz_files:
        name = f["name"]
        matched = False

        for pattern in patterns:
            match = pattern.search(name)
            if match:
                try:
                    # Attempt to extract easting and northing, considering edge cases
                    easting = int(match.group(1))
                    northing = int(match.group(2))

                    # Sanity check for values (e.g., avoid interpreting a date as a coordinate)
                    if (
                        easting < 100000 or northing < 100000
                    ):  # Assuming reasonable limits for easting/northing
                        continue

                    tile_coords.append((easting, northing))
                    row_tiles[northing].append(easting)
                    col_tiles[easting].append(northing)
                    matched = True
                except ValueError:
                    # Handle cases where conversion to int fails
                    continue

        # If no match, continue with the next file
        if not matched:
            continue

    # 5. Sort the rows and columns
    for northing in row_tiles:
        row_tiles[northing].sort()
    for easting in col_tiles:
        col_tiles[easting].sort()

    # 6. Calculate tile dimensions
    tile_sizes = {}
    for easting, northing in tile_coords:
        # Calculate width
        eastings = row_tiles[northing]
        e_index = eastings.index(easting)
        if e_index < len(eastings) - 1:
            width = eastings[e_index + 1] - easting
        elif e_index > 0:
            width = easting - eastings[e_index - 1]
        else:
            width = 0  # Default or handle as needed

        # Calculate height
        northings = col_tiles[easting]
        n_index = northings.index(northing)
        if n_index < len(northings) - 1:
            height = northings[n_index + 1] - northing
        elif n_index > 0:
            height = northing - northings[n_index - 1]
        else:
            height = 0  # Default or handle as needed

        tile_sizes[(easting, northing)] = (width, height)

    # 7. Construct tile geometries and filter by AOI intersection
    records = []
    for f in laz_files:
        name = f["name"]
        matched = False

        for pattern in patterns:
            match = pattern.search(name)
            if match:
                try:
                    # Attempt to extract easting and northing, considering edge cases
                    easting = int(match.group(1))
                    northing = int(match.group(2))

                    # Sanity check for values
                    if (
                        easting < 100000 or northing < 100000
                    ):  # Skip invalid coordinates
                        continue

                    width, height = tile_sizes.get((easting, northing), (0, 0))
                    tile = box(easting, northing, easting + width, northing + height)
                    if tile.intersects(aoi_union):
                        records.append(
                            {"name": name, "url": f["url"], "geometry": tile}
                        )
                    matched = True
                except ValueError:
                    # Handle cases where conversion to int fails
                    continue

        # If no match, continue with the next file
        if not matched:
            continue

    if not records:
        msg_no_lpc_hits = "No LPC tiles intersect your AOI."
        raise RuntimeError(msg_no_lpc_hits)

    # 8. Build GeoDataFrame & reproject
    gdf = gpd.GeoDataFrame(records, crs=utm_crs)
    gdf_wgs84 = gdf.to_crs(epsg=4326)

    # 9. Optionally write out
    if output_dir:
        out_path = Path(output_dir) / f"NOAA_{dataset_id}_lpc_tiles.geojson"
        gdf_wgs84.to_file(out_path, driver="GeoJSON")

    return gdf_wgs84


def download_ncalm_dem(
    aoi: gpd.GeoDataFrame,
    dataset_id: str | int,
    product: str,
    output_dir: str = "/tmp",
) -> None:
    """
    Download NCALM (or OpenTopo user-hosted) data from OpenTopography S3 based on a dataset identifier and an AOI.
    Your dataset identifier will be the 'name' column from coincident.search.search(dataset="ncalm")

    NOTE: NCALM provides both DSMs and DTMs
    e.g. WA18_Wall/WA18_Wall_be/WALL_GEG_1M.tif vs WA18_Wall/WA18_Wall_hh/WALL_GEF_1M.tif
    where the be suffix stands for "bare earth" and the hh suffix stands for "highest hits" (or first return)
    GEG = Grid Elevation (Ground) or the the bare earth elevation grid
    GEF = Grid Elevation (First return) or the the first return elevation grid

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
    product : str
        'dtm' to download the bare-earth DEM (GEG),
        'dsm' to download the first-return DEM (GEF).
    output_dir : str, optional
        Directory to save output files. Defaults to "/tmp".

    Returns
    -------
    None
        The function writes the cropped/clipped output files to the specified output directory.
    """
    product = product.lower()
    if product not in ("dtm", "dsm"):
        msg_invalid_product = "product must be either 'dtm' or 'dsm'"
        raise ValueError(msg_invalid_product)

    # Determine folder and filename suffix
    folder_suffix = "_be" if product == "dtm" else "_hh"
    # file_code = "GEG" if product == "dtm" else "GEF"
    bucket = "raster"
    file_ext = ".tif"

    # Ensure the output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # S3 client for unsigned access
    endpoint_url = "https://opentopography.s3.sdsc.edu"
    s3 = boto3.client(
        "s3",
        config=Config(signature_version=UNSIGNED),
        endpoint_url=endpoint_url,
    )

    # List objects under the dataset prefix + return folder
    prefix = f"{dataset_id}/{dataset_id}{folder_suffix}/"
    response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    if "Contents" not in response:
        msg_empty_project = f"No objects found for prefix {prefix} in bucket {bucket}."
        raise ValueError(msg_empty_project)

    # Estimate the AOI's UTM and reproject
    aoi_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(aoi_crs)

    # Filter keys by suffix
    dem_keys = [
        obj["Key"] for obj in response["Contents"] if obj["Key"].endswith(file_ext)
    ]

    for dem_key in dem_keys:
        dem_url = f"{endpoint_url}/{bucket}/{dem_key}"
        dem = rxr.open_rasterio(dem_url, masked=True)

        # Reproject if needed
        if dem.rio.crs != aoi_crs:
            dem = dem.rio.reproject(aoi_crs)

        # Clip to AOI
        dem_clipped = dem.rio.clip(aoi_utm.geometry, aoi_utm.crs)

        # Save output
        output_name = f"clipped_{Path(dem_key).stem}_{product}.tif"
        clipped_output_path = Path(output_dir) / output_name
        dem_clipped.rio.to_raster(clipped_output_path)


def _fetch_ncalm_lpc_tiles(
    aoi: gpd.GeoDataFrame,
    dataset_name: str,
    output_dir: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Fetch NCALM LPC (1 km LAZ) tiles intersecting an AOI from the OpenTopography 'pc-bulk' bucket,
    optionally downloading them to disk.

    Steps:
      1. Reproject AOI to its local UTM CRS for accurate spatial filtering.
      2. List all .laz keys under pc-bulk/{dataset_id}/.
      3. Client-side filter: extract easting/northing from each filename, build a 1 km box,
         keep only tiles that intersect the AOI in UTM.
      4. Build a GeoDataFrame with columns ['key','url','geometry'] and reproject back to EPSG:4326.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        AOI geometries (any CRS, but assumed geographic or projected).
    dataset_name : str
        NCALM project shortname (e.g. "WA18_Wall").
    output_dir : str | None
        If provided, directory to download .laz files into.

    Returns
    -------
    geopandas.GeoDataFrame
        Columns:
          - key: S3 object key (e.g. "WA18_Wall/617000_5148000.laz")
          - url: full http URL to the LAZ tile
          - geometry: shapely Polygon of the 1 km tile, in EPSG:4326
    """
    # 1. Determine AOI's UTM CRS and reproject
    utm_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(utm_crs)
    union_utm = aoi_utm.union_all()

    # 2. List all .laz keys in pc-bulk
    bucket = "pc-bulk"
    endpoint_url = "https://opentopography.s3.sdsc.edu"
    s3 = boto3.client(
        "s3",
        config=Config(signature_version=UNSIGNED),
        endpoint_url=endpoint_url,
    )
    prefix = f"{dataset_name}/"
    resp = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    if "Contents" not in resp:
        msg_empty_contents = f"No objects found for prefix {prefix} in bucket {bucket}."
        raise ValueError(msg_empty_contents)

    # 3. Filter client-side by 1 km tile intersection
    pattern = re.compile(r"(\d+)\_(\d+)\.laz$")
    records = []
    for obj in resp["Contents"]:
        key = obj["Key"]
        if not key.endswith(".laz"):
            continue
        m = pattern.search(Path(key).name)
        if not m:
            continue
        easting, northing = map(int, m.groups())
        # 1000 for 1km x 1km tile shape (consistent for all NCALM datasets)
        tile = box(easting, northing, easting + 1000, northing + 1000)
        if not tile.intersects(union_utm):
            continue
        url = f"{endpoint_url}/{bucket}/{key}"
        records.append({"name": key, "url": url, "geometry": tile})

    if not records:
        msg_no_overlay = "No LPC tiles intersect your AOI for this dataset."
        raise RuntimeError(msg_no_overlay)

    # 4. Build GeoDataFrame and reproject to 4326
    gdf = gpd.GeoDataFrame(records, crs=utm_crs)
    gdf = gdf.to_crs(epsg=4326)

    # 5. Optionally download
    if output_dir:
        # Path(output_dir).mkdir(parents=True, exist_ok=True)
        for key in gdf["key"]:
            local_path = Path(output_dir) / Path(key).name
            s3.download_file(bucket, key, str(local_path))

    # 6. Optionally save GeoJSON
    if output_dir:
        geojson_path = Path(output_dir) / f"{dataset_name}_lpc_tiles.geojson"
        gdf.to_file(geojson_path, driver="GeoJSON")

    return gdf


def fetch_lpc_tiles(
    aoi: gpd.GeoDataFrame,
    dataset_id: str,
    provider: str,
    datetime_str: str | None = None,
    output_dir: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Unified fetch for LPC tiles across providers: USGS, NOAA, NCALM, NEON.

    Parameters
    ----------
    aoi : gpd.GeoDataFrame
        Area of interest.
    dataset_id : str
        Identifier string for the dataset of interest
          - NOAA: dataset ID, all digits length 1-5. e.g. 10149
          - NEON: dataset ID,4-character uppercase string. e.g. 'ARIK'
          - USGS: 3DEP project name, not to be confused with the workunit, e.g. 'CO_CentralEasternPlains_2020_D20'
          - NCALM: dataset name, not to be confused with the dataset ID, e.g. 'WA18_Wall'
    provider : str
        Provider of LPC data. One of 'usgs', 'noaa', 'ncalm', 'neon'.
    datetime_str : str, optional
        For NEON provider only, a date in 'YYYY-MM-DD' format.
    output_dir : str, optional
        Directory to write provider-specific GeoJSON output.

    Returns
    -------
    gpd.GeoDataFrame
        Tile geometries with 'name', 'url', and 'geometry'.

    Notes
    -----
    If provider is 'neon', you must supply datetime_str.
    """
    # provider
    prov = provider.lower()

    # given the similarities between usgs 3dep and ncalm ids (hard to regex their differences)
    # so we'll just try both ncalm and usgs if not noaa or neon
    if prov == "noaa":
        return _fetch_noaa_lpc_tiles(aoi, dataset_id=dataset_id, output_dir=output_dir)
    if prov == "neon":
        if not datetime_str:
            msg_no_neon_date = "datetime_str is required for NEON provider"
            raise ValueError(msg_no_neon_date)
        return _fetch_neon_lpc_tiles(
            aoi, datetime_str=datetime_str, site_id=dataset_id, output_dir=output_dir
        )
    if prov == "ncalm":
        return _fetch_ncalm_lpc_tiles(
            aoi, dataset_name=dataset_id, output_dir=output_dir
        )
    if prov == "usgs":
        return _fetch_usgs_lpc_tiles(aoi, project=dataset_id, output_dir=output_dir)

    msg_invalid_search = "Invalid search, please double check your parameters and function documentation: 'coincident.io.download.fetch_lpc_tiles?'"
    raise ValueError(msg_invalid_search)


# TODO: add similar logic to the download_urls func instead so every download func is generalizeable?
def _process_gliht_files(
    urls: list[str],
    output_dir: str,
    products: list[str],
    aoi: gpd.GeoDataFrame,
    clip: bool,
    res: int,
) -> None:
    """Helper function to downloads G-LiHT files with optional clipping and coarsening."""
    import tempfile

    # If no processing needed, just download normally
    if not clip and res == 1:
        download_files(
            urls, output_dir, f"G-LiHT {', '.join(products)}", earthdata_auth=True
        )
        return

    # Get Earthdata credentials
    username = os.getenv("EARTHDATA_USER")
    password = os.getenv("EARTHDATA_PASS")
    if not username or not password:
        msg_no_ed_cres = """
        EARTHDATA_USER and EARTHDATA_PASS must be set in the environment.
        os.environ["EARTHDATA_USER"] = "username"
        os.environ["EARTHDATA_PASS"] ="password"
        """
        raise OSError(msg_no_ed_cres)

    for url in tqdm(urls, desc=f"Processing G-LiHT {', '.join(products)} tiles..."):
        filename = Path(url).name
        prefix = []
        if clip:
            prefix.append("clipped")
        if res > 1:
            prefix.append(f"coarsened_{res}x")
        output_filename = f"{'_'.join(prefix) + '_' if prefix else ''}{filename}"
        output_file = Path(output_dir) / output_filename
        if output_file.exists():
            continue

        temp_path = None
        try:
            # Download to temporary file first
            with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp_file:
                temp_path = tmp_file.name

            resp = requests.get(url, stream=True, auth=(username, password))
            resp.raise_for_status()

            with Path(temp_path).open("wb") as f:
                for chunk in resp.iter_content(chunk_size=65536):
                    if chunk:
                        f.write(chunk)

            # Process the local file
            data = rxr.open_rasterio(temp_path, masked=True, chunks="auto")
            # Clip to AOI geometry if requested
            if clip:
                aoi_clipped = (
                    aoi.to_crs(data.rio.crs) if aoi.crs != data.rio.crs else aoi
                )
                data = data.rio.clip(aoi_clipped.geometry, data.rio.crs, drop=True)
            # Optionally coarsen the data
            if res > 1:
                data = data.coarsen(x=res, y=res, boundary="trim").mean()
            data.rio.to_raster(output_file)

        except Exception as error:
            msg = f"Unable to process tile from {url}: {error}"
            logger.error(msg)
            raise RuntimeError(msg) from error
        finally:
            # Clean up temp file
            if temp_path and Path(temp_path).exists():
                Path(temp_path).unlink()


def download_gliht_raster(
    aoi: gpd.GeoDataFrame,
    dataset_id: str,
    products: list[str] | None = None,
    output_dir: str = "/tmp",
    save_parquet: bool = True,
    clip: bool = True,
    res: int = 1,
) -> None:
    """
    Download G-LiHT raster products (e.g., CHM, DTM, DSM) for a given AOI and date.
    This function queries NASA's STAC endpoint for the appropriate G-LiHT product(s),
    downloads the raster COG files via Earthdata authentication, and saves them
    to the specified output directory.

    Supported products:
        - 'ortho': Orthorectified aerial imagery
        - 'chm': Canopy height model
        - 'dsm': Digital surface model
        - 'dtm': Digital terrain model
        - 'hyperspectral_ancillary': Ancillary HSI data
        - 'radiance': Hyperspectral aradiance
        - 'reflectance': Hyperspectral surface reflectance
        - 'hyperspectral_vegetation': HSI-derived veg indices

    NOTE: Not all G-LiHT flights will contain every single product listed
    NOTE: Lidar point cloud data is not supported due to coincident's lack of PDAL support

    Parameters
    ----------
    aoi : gpd.GeoDataFrame
        Area of interest to intersect with G-LiHT tiles.
    datetime_str : str
        Desired acquisition date (format: YYYY or YYYY-MM-DD).
        Provided in the 'datetime' field returned by coincident.search.search
    product_key : str | list[str]
        Key or list of keys to select products from GLIHT_COLLECTIONS_MAP.
    output_dir : str, optional
        Directory path where tiles will be downloaded (default "/tmp").
    save_parquet : bool, optional
        Whether to save metadata as geoparquet (default True).
    clip: bool
        Whether to clip result to the AOI geometry.
    res: int
        Optional coarsening factor for x/y resolution.

    Returns
    -------
    None
        Files are downloaded to the specified output path.
    """
    GLIHT_COLLECTIONS_MAP = {
        "ortho": "GLORTHO_001",  # Orthorectified aerial imagery
        "chm": "GLCHMT_001",  # Tiled version of CHM
        "dsm": "GLDSMT_001",  # Digital surface model (tiled)
        "dtm": "GLDTMT_001",  # Tiled version of DTM
        "hyperspectral_ancillary": "GLHYANC_001",  # Ancillary HSI data
        "radiance": "GLRADS_001",  # Hyperspectral aradiance
        "reflectance": "GLREFL_001",  # Hyperspectral surface reflectance
        "hyperspectral_vegetation": "GLHYVI_001",  # HSI-derived veg indices
    }
    if products is None:
        products = list(GLIHT_COLLECTIONS_MAP.keys())  # Use all available products
        log_all_keys = "No G-LiHT products specified... Downloading all products in 'supported products'"
        logger.warning(log_all_keys)
    # Validate each product in the list
    for product in products:
        if product not in GLIHT_COLLECTIONS_MAP:
            msg_product_key = f"'{product}' is not a valid G-LiHT product key."
            raise ValueError(msg_product_key)

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Process each valid product key
    all_download_urls = []
    all_metadata_items = []

    for key in products:
        collection_id = GLIHT_COLLECTIONS_MAP[key]
        gliht_dataset = GLiHT(collections=[collection_id])

        # Search for intersecting items
        results = search(dataset=gliht_dataset, intersects=aoi)
        possible_ids = _translate_gliht_id(dataset_id, key, GLIHT_COLLECTIONS_MAP)

        # Try to find matching item with any of the possible IDs
        matching_results = None
        for translated_id in possible_ids:
            temp_results = results[results.id == translated_id]
            if not temp_results.empty:
                matching_results = temp_results
                break

        if matching_results is None or matching_results.empty:
            log_empty_results = "No G-LiHT raster tiles found for product '%s' with the given AOI and any of these IDs: %s. Skipping..."
            logger.warning(
                log_empty_results,
                key,
                possible_ids,
            )
            continue

        # Use the matching results
        results = matching_results

        # Extract download URLs and metadata for this product
        for _, item in results.iterrows():
            asset_items = item.assets
            # Find the COG (TIF) asset
            cog_key = next(
                k for k, v in asset_items.items() if v["href"].endswith(".tif")
            )
            cog_href = asset_items[cog_key]["href"]
            all_download_urls.append(cog_href)

            # Collect metadata for parquet file
            all_metadata_items.append(
                {
                    "id": item.id,
                    "collection": item.collection,
                    "datetime": item.datetime,
                    "product_key": key,
                    "asset_key": cog_key,
                    "download_url": cog_href,
                    "bbox": item.bbox,
                }
            )

    # Check if we found any tiles to download
    if not all_download_urls:
        msg_empty_all = (
            "No G-LiHT raster tiles found for any of the valid product keys."
        )
        raise ValueError(msg_empty_all)

    # Optionally save metadata as geoparquet
    if save_parquet:
        geometries = [
            box(
                item["bbox"]["xmin"],
                item["bbox"]["ymin"],
                item["bbox"]["xmax"],
                item["bbox"]["ymax"],
            )
            for item in all_metadata_items
        ]

        gdf = gpd.GeoDataFrame(all_metadata_items, geometry=geometries, crs="EPSG:4326")
        product_keys_str = "_".join(products)
        parquet_path = Path(output_dir) / f"stac_gliht_{product_keys_str}.parquet"
        gdf.to_parquet(parquet_path)

    # Process files with optional clipping and coarsening
    _process_gliht_files(all_download_urls, output_dir, products, aoi, clip, res)
