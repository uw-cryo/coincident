"""
Subset and Sample using Xarray
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import Any

import geopandas as gpd
import odc.stac
import pystac
import rasterio

# NOTE: must import for odc.stac outputs to have .rio accessor
import rioxarray
import xarray as xr
from shapely.geometry import MultiPolygon, Polygon, box

from coincident.datasets.nasa import GLiHT
from coincident.search import search
from coincident.search.neon_api import (
    _get_tile_bbox,
    query_neon_data_api,
)
from coincident.search.opentopo_api import _find_noaa_dem_files
from coincident.search.stac import to_pystac_items
from coincident.search.wesm import query_tnm_api

# Sets GDAL_DISABLE_READDIR_ON_OPEN to 'EMPTY_DIR' etc.
odc.stac.configure_rio(cloud_defaults=True, VSICURL_PC_URL_SIGNING="YES")

# TODO: investigate performance of GTI versus VRT
cop_to_7912 = "https://raw.githubusercontent.com/uw-cryo/3D_CRS_Transformation_Resources/refs/heads/main/globaldems/COP30_hh_7912.vrt"
nasadem_to_7912 = "https://raw.githubusercontent.com/uw-cryo/3D_CRS_Transformation_Resources/refs/heads/main/globaldems/nasadem_7912.vrt"
threedep_to_7912 = "https://raw.githubusercontent.com/uw-cryo/3D_CRS_Transformation_Resources/refs/heads/main/3dep/USGS_Seamless_DEM_13_7912.vrt"


def load_dem_7912(dataset: str, aoi: gpd.GeoDataFrame) -> xr.DataArray:
    """
    Load a Digital Elevation Model (DEM) dataset and clip it to the area of interest (AOI).
    Reprojection-on-reading to EPSG:7912 is performed by GDAL.

    Parameters
    ----------
    dataset : str
        The name of the DEM dataset to load. Must be one of 'cop30', 'nasadem', or '3dep'.
    aoi : gpd.GeoDataFrame
        Polygon with lon,lat coordinates to which the DEM will be clipped.

    Returns
    -------
    xr.DataArray
        The clipped DEM data as an xarray DataArray.

    Notes
    -----
    - This function currently eagerly pulls all data into memory, so may fail if aoi is large
    """
    # NOTE: Currently only getting data from OpenTopography VRTs, but leaving this in for future
    # if dataset in ["cop30", "nasadem"]:
    #     r = requests.get(
    #         "https://planetarycomputer.microsoft.com/api/sas/v1/token/pcstacitems/items",
    #         timeout=10,
    #     )
    #     token = r.json()["token"]
    #     if dataset == "cop30":
    #         uri = cop_to_7912
    #     elif dataset == "nasadem":
    #         uri = nasadem_to_7912
    #     ENV = rasterio.Env(
    #         GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    #         AZURE_STORAGE_ACCOUNT="pcstacitems",
    #         AZURE_STORAGE_SAS_TOKEN=token,
    #         VSICURL_PC_URL_SIGNING="YES",
    #     )
    if dataset == "cop30":
        uri = cop_to_7912
    elif dataset == "nasadem":
        uri = nasadem_to_7912
    elif dataset == "3dep":
        uri = threedep_to_7912
    else:
        msg = f"Unknown dataset: {dataset}, must be 'cop30', 'nasadem', or '3dep'"
        raise ValueError(msg)

    ENV = rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", AWS_NO_SIGN_REQUEST="YES"
    )

    # Pull all data into memory
    with ENV:
        return (
            xr.open_dataarray(uri, engine="rasterio")
            .rio.clip_box(**aoi.bounds)
            .squeeze()
            .load()
        )


def to_dataset(
    gf: gpd.GeoDataFrame,
    bands: list[str] | None = None,
    aoi: gpd.GeoDataFrame | None = None,
    mask: bool = False,
    **kwargs: Any,
) -> xr.DataArray:
    """
    Convert a GeoDataFrame to an xarray DataArray using odc.stac

    Parameters
    ----------
    gf : gpd.GeoDataFrame
        The GeoDataFrame containing the geospatial data.
    aoi : gpd.GeoDataFrame, optional
        GeoDataFrame containing area of interest (e.g. search polygon)
    bands : str, optional
        The asset keys to extract from the STAC items, by default ["data"].
    **kwargs : dict
        Additional keyword arguments passed to `odc.stac.load`.

    Returns
    -------
    xr.Dataset
        The resulting stacked xarray Dataset backed by dask arrays
    """
    if bands is None:
        bands = ["data"]

    items = to_pystac_items(gf)
    ds = odc.stac.load(items, bands=bands, geopolygon=aoi, chunks={}, **kwargs)
    if mask and aoi is not None:
        ds = ds.rio.clip(aoi.geometry)

    return ds


# NOTE: make this more general to open any single STAC asset?
def open_maxar_browse(item: pystac.Item, overview_level: int = 0) -> xr.DataArray:
    """
    Open a browse image from a STAC item using the specified overview level.

    Parameters
    ----------
    item : pystac.Item
        The STAC item containing the browse image asset.
    overview_level : int, optional
        The overview level to use when opening the image, by default 0.

    Returns
    -------
    xr.DataArray
        The opened browse image as an xarray DataArray.

    Notes
    -----
    The function uses the `rasterio` engine to open the image and sets the
    `GDAL_DISABLE_READDIR_ON_OPEN` and `GDAL_HTTP_HEADERS` environment variables
    for optimized reading and authentication, respectively.
    """

    # href = item.assets['browse']['href'] # accessed via geopandas
    href = item.assets["browse"].href  # accessed via pystac

    env = rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        GDAL_HTTP_HEADERS=f"MAXAR-API-KEY:{os.environ['MAXAR_API_KEY']}",
    )
    with env:
        return xr.open_dataarray(
            href,
            engine="rasterio",
            mask_and_scale=False,  # otherwise uint8 -> float32!
            backend_kwargs={"open_kwargs": {"overview_level": overview_level}},
        )


def _aoi_to_polygon_string(geom: Polygon | MultiPolygon) -> str:
    """
    HELPER FUNCTION FOR load_usgs_dem()
    Convert a shapely geometry into a polygon string required by the TNM API.

    This updated helper now accepts a shapely geometry (as produced in search.search).

    Parameters:
      geom (Polygon or MultiPolygon): Geometry defining the AOI.

    Returns:
      str: A string in the format "x1 y1, x2 y2, ..." representing the polygon coordinates.
    """
    if geom.geom_type == "Polygon":
        coords = list(geom.exterior.coords)
    elif geom.geom_type == "MultiPolygon":
        coords = list(geom.convex_hull.exterior.coords)
    else:
        msg_aoi_err = "Unsupported geometry type for AOI"
        raise ValueError(msg_aoi_err)
    return ", ".join(f"{x} {y}" for x, y in coords)


def _filter_items_by_project(
    api_items: list[dict[str, Any]], project_identifier: str
) -> list[dict[str, Any]]:
    """
    HELPER FUNCTION FOR load_usgs_dem()
    Filter TNM API items to only include those corresponding to a specific project.

    Parameters:
      api_items (list): List of items (dictionaries) returned by the TNM API.
      project_identifier (str): The project identifier to filter by.

    Returns:
      list: A filtered list of API items matching the specified project.
    """
    filtered_items = []
    for item in api_items:
        vendor_meta_url = item.get("vendorMetaUrl", "")
        if vendor_meta_url:
            try:
                remainder = vendor_meta_url.split("metadata/")[-1]
                parts = remainder.split("/")
                if len(parts) >= 2:
                    candidate_project = parts[0]
                    if candidate_project == project_identifier:
                        filtered_items.append(item)
            except Exception as e:
                msg_vendor_url = (
                    f"Unable to parse vendorMetaUrl '{vendor_meta_url}'. Error: {e}"
                )
                raise RuntimeError(msg_vendor_url) from e

    return filtered_items


# NOTE: see https://apps.nationalmap.gov/tnmaccess/#/product for more datasets we can load in
# LPCs and 5m IFSAR DSMs might be of interest
# but we'll have to do something with this wonky "_filter_items_by_project" function i made if
# we add datasets
def load_usgs_dem(
    aoi: gpd.GeoDataFrame,
    project: str,
    tnmdataset: str = "Digital Elevation Model (DEM) 1 meter",
    res: int = 1,
    clip: bool = True,
) -> xr.DataArray:
    """
    Load and merge USGS 1-meter DEM tiles based on an AOI by querying the TNM API.

    Steps:
      1. Reproject the AOI to EPSG:4326.
      2. Extract the first geometry from the exploded AOI (to match the geometry type used in search.search).
      3. Convert the geometry to a polygon string using a private helper.
      4. Query the TNM API via the moved function in coincident.search.wesm.
      5. Filter the API items using a private helper.
      6. Load and optionally coarsen each GeoTIFF tile.
      7. Merge the DEM tiles and optionally clip the mosaic to the AOI.

    Parameters:
        aoi (gpd.GeoDataFrame): Area of interest geometry to query against
        project (str): Project identifier to filter results
        tnmdataset (str): TNM dataset identifier (default "Digital Elevation Model (DEM) 1 meter")
        res (int): Resolution factor to coarsen DEM by (default 1)
        clip (bool): Whether to clip final mosaic to AOI (default True)

    Returns:
      xr.DataArray: The merged (and optionally clipped) DEM mosaic.
    """
    # 1: Reproject AOI to EPSG:4326 for the TNM API query.
    aoi_in_4326 = aoi.to_crs(epsg=4326)

    # 2: Explode and extract the first geometry.
    shapely_geom = aoi_in_4326.geometry.explode().iloc[0]

    # 3: Convert geometry to polygon string using the updated helper function.
    polygon_string = _aoi_to_polygon_string(shapely_geom)

    # 4: Query the TNM API (function moved to coincident.search.wesm).
    tnm_api_items = query_tnm_api(polygon_string, tnmdataset=tnmdataset)
    if not tnm_api_items:
        msg_empty_tnm = "TNM API returned no products for the given AOI and dataset."
        raise ValueError(msg_empty_tnm)

    # 5: Filter the returned items using the private helper.
    filtered_api_items = _filter_items_by_project(tnm_api_items, project)
    if not filtered_api_items:
        msg_no_project_id = (
            "No products found for input project string intersecting the AOI."
        )
        raise ValueError(msg_no_project_id)

    # 6: Extract TIFF URLs and open each GeoTIFF lazily.
    tiff_urls = [
        item["urls"]["TIFF"]
        for item in filtered_api_items
        if "urls" in item and "TIFF" in item["urls"]
    ]
    dem_dataarrays = []
    for url in tiff_urls:
        try:
            dem_tile = rioxarray.open_rasterio(url, masked=True, chunks="auto")
            if res > 1:
                dem_tile = dem_tile.coarsen(x=res, y=res, boundary="trim").mean()
            dem_dataarrays.append(dem_tile)
        except Exception as error:
            msg_vendor_err = f"Unable to parse vendorMetaUrl. Error: {error}"
            raise RuntimeError(msg_vendor_err) from error
    if not dem_dataarrays:
        msg_tiff_open = "None of the TIFF files could be opened."
        raise RuntimeError(msg_tiff_open)

    # 7: Merge DEM tiles and clip if required.
    tile_crs = dem_dataarrays[0].rio.crs
    merged_dem = (
        xr.concat(dem_dataarrays, dim="band").max(dim="band").rename("elevation")
    )
    if clip:
        aoi_tile_crs = aoi.to_crs(tile_crs)
        merged_dem = merged_dem.rio.clip(aoi_tile_crs.geometry, tile_crs, drop=True)

    return merged_dem


# TODO: look into warning below that prints for every "for f in filtered_files:"
# RuntimeWarning: TIFFReadDirectoryCheckOrder:Invalid TIFF directory; tags are not sorted in ascending order if riods.subdatasets:
# NOTE: this function requires client-side spatial querying. I went with regex for this because it's significantly faster
# than reading in metadata headers via rasterio as seen in load_noaa_dem but we can always switch
def load_neon_dem(
    aoi: gpd.GeoDataFrame,
    datetime_str: str,
    site_id: str,
    product: str,
    res: int = 1,
    clip: bool = True,
) -> xr.DataArray:
    """
    Load and merge NEON LiDAR tiles (DSM, DTM, or CHM) based on an AOI by querying the NEON API.

    Steps:
        1. Convert the datetime string to a month string in the format YYYY-MM.
        2. Determine appropriate UTM CRS for the AOI.
        3. Query the NEON API using a preset product code.
        4. Filter the returned files based on product type and spatial intersection.
        5. Load and optionally coarsen each GeoTIFF tile.
        6. Merge the tiles and optionally clip the mosaic to the AOI.

    Parameters:
        aoi (gpd.GeoDataFrame): Area of interest geometry to query against
        datetime_str (str): Date string in YYYY-MM-DD format
        site_id (str): NEON site identifier
        product (str): Product type to load ('dsm', 'dtm', or 'chm')
        res (int): Resolution factor to coarsen DEM by (default 1)
        clip (bool): Whether to clip final mosaic to AOI (default True)

    Returns:
        xr.DataArray: The merged (and optionally clipped) LiDAR mosaic
    """
    # 1: Convert datetime_str to month string
    dt = datetime.strptime(datetime_str, "%Y-%m-%d")
    month_str = dt.strftime("%Y-%m")

    # 2: Determine appropriate UTM CRS for the AOI
    target_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(target_crs)
    aoi_geom_utm = aoi_utm.union_all()

    # 3: Query the NEON API
    data_json = query_neon_data_api(site_id, month_str)
    if (
        data_json.get("data", {}).get("release") == "PROVISIONAL"
        and len(data_json.get("data", {}).get("files", [])) == 0
    ):
        msg_provisional = (
            f"Data for {site_id} in {month_str} has PROVISIONAL status with no files."
        )
        raise ValueError(msg_provisional)

    # 4: Spatial filtering: filter for product TIFF files and intersect with AOI
    files = data_json["data"]["files"]
    filtered_files = []
    # Filter based on product (e.g., 'DSM.tif', 'DTM.tif', or 'CHM.tif')
    product_filter = f"{product.upper()}.tif"
    for f in files:
        if product_filter in f["name"]:
            tile_bbox = _get_tile_bbox(f["name"])
            # Ensure we got a valid tile bounding box and that it intersects the AOI
            if tile_bbox is not None and tile_bbox.intersects(aoi_geom_utm):
                filtered_files.append(f)
    if not filtered_files:
        msg_neon_empty = "No matching LiDAR files found for the given spatial query."
        raise RuntimeError(msg_neon_empty)

    # 5: Set an optional resolution factor for coarsening. If res > 1, the tiles are coarsened.
    dem_dataarrays = []
    for f in filtered_files:
        try:
            da = rioxarray.open_rasterio(f["url"], masked=True, chunks="auto")
            if res > 1:
                da = da.coarsen(x=res, y=res, boundary="trim").mean()
            dem_dataarrays.append(da)
        except Exception as error:
            msg_neon_open = f"Unable to open file {f['name']}. Error: {error}"
            raise RuntimeError(msg_neon_open) from error
    if not dem_dataarrays:
        msg_neon_tif_open = "None of the TIFF files could be opened."
        raise RuntimeError(msg_neon_tif_open)

    # 6: Merge the tiles:
    # Concatenate along a new "band" dimension and then take the max across bands.
    tile_crs = dem_dataarrays[0].rio.crs
    merged = xr.concat(dem_dataarrays, dim="band").max(dim="band").rename("elevation")
    # Reproject the AOI to the tile CRS and clip the mosaic
    if clip:
        aoi_tile_crs = aoi.to_crs(tile_crs)
        merged = merged.rio.clip(aoi_tile_crs.geometry, tile_crs, drop=True)

    return merged


def load_ncalm_dem(
    aoi: gpd.GeoDataFrame,
    dataset_id: str | int,
    product: str,
    res: int = 1,
    clip: bool = True,
) -> xr.DataArray:
    """
    Load NCALM DEM (DTM or DSM) tiles directly from OpenTopography S3 and return as xarray DataArray.

    NCALM provides:
      - DTM (bare-earth) with file code GEG in folder suffix '_be'
      - DSM (first-return) with file code GEF in folder suffix '_hh'

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest geometries (in EPSG:4326).
    dataset_id : str or int
        NCALM dataset identifier (e.g., 'WA18_Wall').
    product : str
        'dtm' for bare-earth, 'dsm' for first-return.
    res : int, optional
        Factor to coarsen the DEM (default 1 = no coarsening).
    clip : bool, optional
        Whether to clip final mosaic to the AOI (default True).

    Returns
    -------
    xarray.DataArray or list of xarray.DataArray
        DEM tile arrays or merged mosaic with name 'elevation'.
    """
    import boto3
    from botocore import UNSIGNED
    from botocore.client import Config

    # Validate product
    product = product.lower()
    if product not in ("dtm", "dsm"):
        msg_invalid_prod = "product must be 'dtm' or 'dsm'"
        raise ValueError(msg_invalid_prod)

    # Determine folder and file code
    folder_suffix = "_be" if product == "dtm" else "_hh"
    # file_code = "GEG" if product == "dtm" else "GEF"

    # S3 setup
    bucket = "raster"
    endpoint_url = "https://opentopography.s3.sdsc.edu"
    s3 = boto3.client(
        "s3",
        config=Config(signature_version=UNSIGNED),
        endpoint_url=endpoint_url,
    )

    # List DEM tiles
    prefix = f"{dataset_id}/{dataset_id}{folder_suffix}/"
    resp = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    if "Contents" not in resp:
        msg_no_contents = f"No objects found for prefix {prefix} in bucket {bucket}."
        raise ValueError(msg_no_contents)

    # AOI UTM reprojection
    aoi_crs = aoi.estimate_utm_crs()
    aoi_utm = aoi.to_crs(aoi_crs)

    # Filter keys matching DEM code and resolution
    pattern = ".tif"
    # file_ext = ".tif"
    # pattern = f"_{file_code}_1M{file_ext}"
    key = next(
        (obj["Key"] for obj in resp["Contents"] if obj["Key"].endswith(pattern)), None
    )
    if not key:
        msg_no_overlay = f"No DEM files ({pattern}) found for dataset {dataset_id}"
        raise ValueError(msg_no_overlay)
    url = f"{endpoint_url}/{bucket}/{key}"
    da = rioxarray.open_rasterio(url, masked=True)
    # select band1
    da = da.sel(band=1).drop_vars("band")

    # reproject if needed
    if da.rio.crs != aoi_crs:
        da = da.rio.reproject(aoi_crs)

    # coarsen if requested
    if res > 1:
        da = da.coarsen(x=res, y=res, boundary="trim").mean()

    # clip if requested
    if clip:
        da = da.rio.clip(aoi_utm.geometry, aoi_crs, drop=True)

    da.name = "elevation"
    return da


# Helper function to fetch a tile's bounds and retain its URL.
def _fetch_noaa_tile_geometry(
    url: str,
) -> tuple[str, Polygon, rasterio.crs.CRS]:
    # Open the raster header from S3 using rasterio.Env for proper configuration.
    with rasterio.Env(), rasterio.open(url) as src:
        bounds = src.bounds
        tile_crs = src.crs
    # Create a shapely bounding box polygon from the raster bounds.
    polygon = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    return url, polygon, tile_crs


# NOTE: NOAA has separate ids and catalogs for their "lidar" and "imagery and elevation raster" datasets
# https://coast.noaa.gov/htdata/lidar1_z/
# https://coast.noaa.gov/htdata/raster1/
# Our current NOAA search functionality only supports the lidar datasets via opentopo
# Note that we get the dataset IDs from our NOAA search for the lidar datasets, not the dem datasets
# that being said the dem ids = lidar ids + 1
# e.g. "Great Bay NERR UAS Lidar" has id 10175 for lidar and id 10176 for dem
# NOTE: NOAA tiles also don't have a consistent size... see below
# 'https://noaa-nos-coastal-lidar-pds.s3.amazonaws.com/dem/NGS_South_TampBay_Topobathy_2021_9481/2021_324000e_3062000n_dem.tif'
# 'https://noaa-nos-coastal-lidar-pds.s3.amazonaws.com/dem/NGS_South_TampBay_Topobathy_2021_9481/2021_364000e_3082000n_dem.tif'
def load_noaa_dem(
    aoi: gpd.GeoDataFrame, dataset_id: str | int, res: int = 1, clip: bool = True
) -> xr.DataArray:
    """
    Process NOAA coastal lidar DEM data from S3 based on a dataset identifier and an AOI.
    Returns an xarray DataArray with elevation values clipped to the AOI without downloading files locally.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest geometry (assumed to be in EPSG:4326 or any other CRS that will be reprojected).
    dataset_id : str or int
        NOAA dataset identifier (e.g., "8431" or "6260").
    res : int, optional
        Resolution factor. If greater than 1, the data will be coarsened by this factor. Defaults to 1.
    clip : bool, optional
        Whether to clip the output mosaic to the AOI. Defaults to True.

    Returns
    -------
    xarray.DataArray
        An xarray DataArray containing DEM values clipped to the AOI.

    Notes
    -----
    This function queries only the raster headers from S3 to retrieve bounds information quickly.
    The header metadata is processed concurrently to reduce I/O latency.
    Tiles are spatially filtered against the AOI and merged into a mosaic using xarray.
    """
    # Step 1: List all .tif files within the NOAA dataset directory
    tif_files = _find_noaa_dem_files(dataset_id)
    tile_urls = [file_info["url"] for file_info in tif_files]

    # Use ThreadPoolExecutor to concurrently fetch header metadata from each tile URL
    with ThreadPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(_fetch_noaa_tile_geometry, tile_urls))

    # Unpack results: each item is a tuple (url, geometry, crs)
    urls, tile_geometries, crs_list = zip(*results, strict=False)
    common_crs = crs_list[0]  # Assume all tiles share the same CRS

    # Create a GeoDataFrame that includes both geometry and the corresponding URL.
    tiles_gdf = gpd.GeoDataFrame(
        {"url": urls, "geometry": tile_geometries}, crs=common_crs
    )

    # --- Preparing the AOI ---
    # Reproject AOI to the common CRS of the raster tiles.
    aoi = aoi.to_crs(common_crs)

    # --- Spatial Filtering ---
    aoi_union = aoi.union_all()
    selected_tiles = tiles_gdf[tiles_gdf.intersects(aoi_union)]

    if selected_tiles.empty:
        msg = (
            f"No TIFF files intersect with the provided AOI for dataset ID {dataset_id}"
        )
        raise ValueError(msg)

    # Step 5: Process each DEM and store the resulting DataArray in a list.
    dem_dataarrays = []
    for _, row in selected_tiles.iterrows():
        try:
            # Open the remote raster using rioxarray without loading full data (only metadata and chunks)
            da = rioxarray.open_rasterio(row["url"], masked=True, chunks="auto")
            # Coarsen the resolution if a factor greater than 1 is provided.
            if res > 1:
                da = da.coarsen(x=res, y=res, boundary="trim").mean()
            dem_dataarrays.append(da)
        except Exception as error:
            msg = f"Unable to open file from URL {row['url']}. Error: {error}"
            raise RuntimeError(msg) from error

    if not dem_dataarrays:
        msg_empty_arrays = "None of the TIFF files could be opened."
        raise RuntimeError(msg_empty_arrays)

    # Step 6: Merge the DEM tiles.
    # Concatenate along a new "band" dimension and use the maximum value across bands to form a mosaic.
    tile_crs = dem_dataarrays[0].rio.crs
    merged_dem = (
        xr.concat(dem_dataarrays, dim="band").max(dim="band").rename("elevation")
    )

    # Step 7: Reproject the AOI to tile CRS and clip the mosaic if requested.
    if clip:
        aoi_tile_crs = aoi.to_crs(tile_crs)
        merged_dem = merged_dem.rio.clip(aoi_tile_crs.geometry, tile_crs, drop=True)

    return merged_dem


def _translate_gliht_id(
    original_id: str, target_product: str, GLIHT_COLLECTIONS_MAP: dict[str, Any]
) -> tuple[str, ...]:
    """
    Helper function for load_gliht_raster to translate IDs between datasets
    Because most people will input the search ID for the default GLiHT collection (metadata)
    rather than the id for what is needed (e.g. dsm vs dtm cvs chm)
    """
    collection_code = GLIHT_COLLECTIONS_MAP[target_product].replace("_001", "")

    # Check if the target collection code is already in the original ID
    if original_id.startswith(collection_code):
        return (original_id,)

    parts = original_id.split("_", 1)  # Split only on first underscore
    # Return both possible formats, sometimes a tail is added sometimes not (e.g. CHM vs DSM)
    format1 = f"{collection_code}_{parts[1]}_{target_product.upper()}"
    format2 = f"{collection_code}_{parts[1]}"
    return format1, format2


def load_gliht_raster(
    aoi: gpd.GeoDataFrame,
    dataset_id: str,
    product: str,
) -> xr.Dataset:
    """
    Load a G-LiHT raster product (e.g., CHM, DTM, DSM) for a given AOI and date.

    This function queries NASA's STAC endpoint for the appropriate G-LiHT product,
    downloads the raster COG(s) via Earthdata authentication, and returns a
    Dask-backed, clipped xarray.Dataset with renamed bands.

    For products with multiple data assets (e.g., DSM with Aspect, Slope, Rugosity),
    all assets will be loaded as separate bands in the output Dataset.

    NOTE: many DSM products from G-LiHT do not have an elevation band provided,
    users should consider searching for CHM products for replacement
    """
    GLIHT_COLLECTIONS_MAP = {
        "ortho": "GLORTHO_001",  # Orthorectified aerial imagery
        "chm": "GLCHMT_001",  # Canopy height model
        "dsm": "GLDSMT_001",  # Digital surface model
        "dtm": "GLDTMT_001",  # Digital terrain model
        "hyperspectral_ancillary": "GLHYANC_001",  # Ancillary HSI data
        "radiance": "GLRADS_001",  # Hyperspectral radiance
        "reflectance": "GLREFL_001",  # Hyperspectral surface reflectance
        "hyperspectral_vegetation": "GLHYVI_001",  # HSI-derived veg indices
    }
    if product not in GLIHT_COLLECTIONS_MAP:
        msg_product_key = f"'{product}' is not a valid G-LiHT product key."
        raise ValueError(msg_product_key)

    collection_id = GLIHT_COLLECTIONS_MAP[product]
    gliht_dataset = GLiHT(collections=[collection_id])

    # Search for intersecting items and find the specific one by its ID
    results = search(dataset=gliht_dataset, intersects=aoi)
    if results.empty:
        msg_empty_gliht = "No GLiHT raster tiles found for the given AOI and id."
        raise ValueError(msg_empty_gliht)

    possible_ids = _translate_gliht_id(dataset_id, product, GLIHT_COLLECTIONS_MAP)
    df_item = None
    for translated_id in possible_ids:
        matching_items = results[results.id == translated_id]
        if not matching_items.empty:
            df_item = matching_items.iloc[[0]]  # Keep as gdf
            break

    if df_item is None:
        msg_no_match = f"No GLiHT item found with any of these IDs: {possible_ids}"
        raise ValueError(msg_no_match)

    asset_items = df_item.assets.item()

    # Find all data assets that are GeoTIFFs
    data_assets = {
        key: asset
        for key, asset in asset_items.items()
        if "href" in asset
        and asset["href"].endswith(".tif")
        and "data" in asset.get("roles", [])
    }
    if not data_assets:
        msg_no_data_assets = f"No .tif data assets found for item {dataset_id}"
        raise ValueError(msg_no_data_assets)

    asset_keys_to_load = list(data_assets.keys())

    # Check for Earthdata credentials in environment
    username = os.getenv("EARTHDATA_USER")
    password = os.getenv("EARTHDATA_PASS")
    if not username or not password:
        msg_no_ed_cres = """
        EARTHDATA_USER and EARTHDATA_PASS must be set in the environment.
        os.environ["EARTHDATA_USER"] = "username"
        os.environ["EARTHDATA_PASS"] ="password"
        """
        raise OSError(msg_no_ed_cres)

    # Use the href from the first found data asset to create the grid template
    first_href = data_assets[asset_keys_to_load[0]]["href"]

    # Use a rasterio.Env context to handle authentication for remote access
    with rasterio.Env(
        GDAL_HTTP_USERNAME=username,
        GDAL_HTTP_PASSWORD=password,
        GDAL_HTTP_COOKIEFILE=Path("~/cookies.txt").expanduser(),
        GDAL_HTTP_COOKIEJAR=Path("~/cookies.txt").expanduser(),
    ):
        # Create a template from the first asset to define the output grid
        template = rioxarray.open_rasterio(
            first_href, chunks={"x": 512, "y": 512}, masked=True
        ).squeeze(drop=True)

        # Load all selected data assets using their keys and the template grid
        ds = to_dataset(
            df_item,
            bands=asset_keys_to_load,
            mask=True,
            like=template,
        )

    # Create a map to rename the long band names to simpler, more readable names
    rename_map = {band_name: band_name.split("_")[-1] for band_name in ds.data_vars}
    return ds.rename(rename_map)
