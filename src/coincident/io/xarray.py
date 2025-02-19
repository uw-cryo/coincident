"""
Subset and Sample using Xarray
"""

from __future__ import annotations

import os
from typing import Any

import geopandas as gpd
import odc.stac
import pystac
import rasterio

# NOTE: must import for odc.stac outputs to have .rio accessor
import rioxarray  # noqa: F401
import xarray as xr

from coincident.search.stac import to_pystac_items

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
