"""
Subset and Sample using Xarray
"""

from __future__ import annotations

import os
import re  # for NEON vendor products, client-side spatial filtering
from datetime import datetime
from typing import Any

import geopandas as gpd
import odc.stac
import pystac
import rasterio

# NOTE: must import for odc.stac outputs to have .rio accessor
import rioxarray
import xarray as xr
from shapely.geometry import MultiPolygon, Polygon, box

from coincident.search.neon_api import query_neon_data_api
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
    if mask:
        ds = ds.rio.clip(aoi)

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
        GDAL_HTTP_HEADERS=f'MAXAR-API-KEY:{os.environ["MAXAR_API_KEY"]}',
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


def _utm_crs_from_lonlat(lon: float, lat: float) -> str:
    """
    HELPER FUNCTION FOR load_neon_dem()
    Compute the UTM CRS for a given longitude and latitude.
    Assumes northern hemisphere if lat >= 0, otherwise southern.
    Returns the CRS as a string, e.g., "EPSG:326XX" or "EPSG:327XX".
    """
    zone = int((lon + 180) / 6) + 1
    epsg_code = 32600 + zone if lat >= 0 else 32700 + zone
    return f"EPSG:{epsg_code}"


def _get_tile_bbox(file_name: str, tile_size: int = 1000) -> Polygon | None:
    """
    HELPER FUNCTION FOR load_neon_dem()
    Extract tile lower-left coordinates from the file name and return a shapely box.
    Assumes file name pattern like: NEON_D17_TEAK_DP3_317000_4104000_CHM.tif

    Note that we have to do client-side spatial filtering for NEON products since
    their API doesn't allow for this. Assumes that all NEON products are in respective
    UTM zones.
    """
    match = re.search(r"_(\d{6})_(\d{7})_", file_name)
    if match:
        easting = int(match.group(1))
        northing = int(match.group(2))
        return box(easting, northing, easting + tile_size, northing + tile_size)

    return None  # return None if no match is found


# TODO: look into warning below that prints for every "for f in filtered_files:"
# RuntimeWarning: TIFFReadDirectoryCheckOrder:Invalid TIFF directory; tags are not sorted in ascending order if riods.subdatasets:
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

    Returns:
      xr.DataArray: The merged (and optionally clipped) LiDAR mosaic.
    """
    # 1: Convert datetime_str to month string
    dt = datetime.strptime(datetime_str, "%Y-%m-%d")
    month_str = dt.strftime("%Y-%m")

    # 2: Determine appropriate UTM CRS for the AOI
    # we need to do this for client-side spatial filtering
    # Get the centroid of the AOI (in EPSG:4326)
    centroid = aoi.union_all().centroid
    target_crs = _utm_crs_from_lonlat(centroid.x, centroid.y)
    # Reproject AOI to target CRS
    aoi_utm = aoi.to_crs(target_crs)
    aoi_geom_utm = aoi_utm.union_all()  # Combined AOI in UTM coordinates

    # 3: Query the NEON API
    data_json = query_neon_data_api(site_id, month_str)

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
