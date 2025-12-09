"""
GDAL Python bindings + PROJ for handling correct 3D Reprojections
"""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

try:
    from osgeo import gdal, gdal_array

    gdal.UseExceptions()
except ImportError:
    warnings.warn(
        "gdal python bindings not found. `pip install gdal` for gdaldem functions",
        stacklevel=2,
    )


def gdaldem(
    da: xr.DataArray, subcommand: str = "hillshade", **kwargs: Any
) -> xr.DataArray:
    """
    Use GDAL python bindings to perform gdaldem operations

    Parameters
    ----------
    da: xarray.DataArray
        Dataarray containing elevation values
    subcommand: str
        'hillshade' (default), 'aspect', 'slope'
    kwargs: dict
        https://gdal.org/en/stable/api/python/utilities.html#osgeo.gdal.DEMProcessingOptions

    Returns
    -------
    da: xarray.DataArray
        XArray dataset with result of running GDAL algorithm
    """
    # Ensure we're working with a 2D Array
    elevation = da.squeeze()
    src = gdal_array.OpenArray(elevation.to_numpy())
    geotransform = da.rio.transform().to_gdal()
    src.SetGeoTransform(geotransform)

    # Automatically set vertical scale
    scale = 111120 if da.rio.crs.is_geographic else 1

    ds_gdal = gdal.DEMProcessing(
        "", src, subcommand, format="MEM", scale=scale, **kwargs
    )
    numpy_result = ds_gdal.ReadAsArray()

    if subcommand != "hillshade":
        numpy_result[numpy_result == -9999] = np.nan

    del ds_gdal

    return elevation.copy(data=numpy_result)


def create_vrt(
    url_list: list[str],
    outputBounds: tuple[float, float, float, float] | None = None,
    resolution: str = "highest",
    custom_res: float | None = None,
    tap: bool = False,
    resampling: str = "nearest",
    a_srs: str | None = None,  # override SRS
    prepend_vsicurl: bool = True,
    ignore_overviews: bool = True,
) -> xr.DataArray:
    """
    Create a VRT (Virtual Raster) from a list of URLs using GDAL.
    https://gdal.org/en/stable/api/python/utilities.html#osgeo.gdal.BuildVRT

    Parameters:
    -----------
    url_list : list of str
        List of URLs pointing to raster files

    outputBounds : tuple
        (minX, minY, maxX, maxY) in target SRS coordinates

    resolution : str
        GDAL VRT resolution option. One of {'highest', 'lowest', 'average', 'user'}

    custom_res : float | None
        User-specified resolution in target SRS units

    tap : bool
        Whether to use target aligned pixels

    resampling : str
        Resampling algorithm to use. One of
        {nearest|bilinear|cubic|cubicspline|lanczos|average|mode}

    a_srs : str
        Override the SRS of the input files with this SRS

    prepend_vsicurl : bool
        Whether to prefix URLs with /vsicurl/ for remote access. If files are not COG
        and high resolution, set to False to use HTTP driver and download entire file up front.

    ignore_overviews : bool
        If True, ignore existing overviews in the source datasets.

    **kwargs : dict
        Options passed to gdal.BuildVRTOptions()
    Returns:
    --------
    str
        Path to the temporary VRT file
    """
    gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")

    # Prefix URLs with vsicurl if they start with http
    url_list = [
        f"/vsicurl/{url}" if url.startswith("http") else url for url in url_list
    ]

    # TAP only used with specified resolution
    if custom_res:
        resolution = "user"
        tap = True

    # Create a named temporary file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".vrt", delete=False) as temp_vrt:
        vrt_path = temp_vrt.name

    # Build VRT from the list of URLs
    # for USGS seamless, 'nearest' resampling should be fine
    vrt_options = gdal.BuildVRTOptions(
        resampleAlg=resampling,
        targetAlignedPixels=tap,
        resolution=resolution,
        xRes=custom_res,
        yRes=custom_res,
        outputBounds=outputBounds,
        outputSRS=a_srs,
        creationOptions=["BLOCKYSIZE=512", "BLOCKXSIZE=512"],
    )
    vrt_ds = gdal.BuildVRT(vrt_path, url_list, options=vrt_options)
    # Close the dataset to flush to disk
    vrt_ds = None  # noqa: F841

    # A bit hacky, but we want to use /vsicurl/ when *constructing* the VRT, but not reading tiff tiles later on
    # Read VRT content once if modifications are needed
    # TODO: alternatively use gdal bindings to modify VRT in memory
    if not prepend_vsicurl or ignore_overviews:
        # print("Modifying VRT content...")
        with Path(vrt_path).open("rt") as fp:
            vrt_content = fp.readlines()

        if not prepend_vsicurl:
            # print("removing /vsicurl/ from VRT...")
            vrt_content = [line.replace("/vsicurl/", "") for line in vrt_content]

        if ignore_overviews:
            # print("removing OverviewList from VRT...")
            vrt_content = [
                line
                for line in vrt_content
                if not line.strip().startswith("<OverviewList")
            ]

        with Path(vrt_path).open("w") as fp:
            fp.writelines(vrt_content)

    return vrt_path


# NOTE: create mapping of common transforms to simplify user input?
# NAD83/NAVD88 <-> ITRF2014 & ITRF2020
# WGS84+EGM2008 <-> ITRF2014 & ITRF2020
def warp_with_pipeline(
    inputData: str,
    srcSRS: str,
    dstSRS: str,
    pipeline: str,
    resampling: str = "bilinear",
) -> xr.DataArray:
    """
    Use a PROJ pipeline to ensure 3D reprojection happens

    Parameters
    ----------
    inputData : str
        Input raster file path or URL
    srcSRS : str
        Source spatial reference system. e.g. 'EPSG:9055+3855' or custom WKT2 string
    dstSRS : str
        Destination spatial reference system. e.g. 'EPSG:7912' or custom WKT2 string
    pipeline : str
        PROJ pipeline string defining the transformation
    resampling : str
        Resampling algorithm to use. One of
        {nearest|bilinear|cubic|cubicspline|lanczos|average|mode}

    Returns
    -------
    xr.DataArray
        DataArray containing the warped raster


    Examples
    -------
    # (projinfo -s EPSG:9055+3855 -t EPSG:7661 -o PROJ --single-line -q)
    PROJ_PIPELINE='+proj=pipeline +step +proj=axisswap +order=2,1 +step +proj=unitconvert +xy_in=deg +xy_out=rad +step +proj=vgridshift +grids=us_nga_egm08_25.tif +multiplier=1 +step +proj=unitconvert +xy_in=rad +xy_out=deg +step +proj=axisswap +order=2,1'
    warp_with_pipeline(input='cop30.tif',
                       srcSRS='EPSG:9055+3855',
                       dstSRS='EPSG:7661',
                       pipeline=PROJ_PIPELINE)
    """
    gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    inputData = f"/vsicurl/{inputData}" if inputData.startswith("http") else inputData

    # Create a named temporary file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".vrt", delete=False) as temp_vrt:
        vrt_path = temp_vrt.name

    # Get input raster properties
    with gdal.Open(inputData) as input_ds:
        geotransform = input_ds.GetGeoTransform()
        xRes = abs(geotransform[1])
        yRes = abs(geotransform[5])

    # 6.71089e+07
    warp_options = gdal.WarpOptions(
        format="VRT",
        srcSRS=srcSRS,
        dstSRS=dstSRS,
        outputType=gdal.GDT_Float32,
        resampleAlg=resampling,
        coordinateOperation=pipeline,
        warpMemoryLimit=1000,  # in MB
        targetAlignedPixels=True,
        xRes=xRes,  # maintain resolution of input
        yRes=yRes,
        multithread=True,  # I don't think this gets stored in VRT?,
        creationOptions=["BLOCKYSIZE=512", "BLOCKXSIZE=512"],
    )
    warped_ds = gdal.Warp(vrt_path, inputData, options=warp_options)  # noqa: F841
    # numpy_result = warped_ds.ReadAsArray()
    # del warped_ds
    # return xr.DataArray(numpy_result)
    return vrt_path


# def load_dem_ellipsoid(dataset: str, aoi: gpd.GeoDataFrame, dstSRS: str = "EPSG:7912") -> xr.DataArray:
#     """
#     Load a Digital Elevation Model (DEM) dataset reprojected to ITRF ellipsoid heights.

#     1. Search for DEM tiles covering the area of interest (AOI).
#     2. Create a VRT from the DEM tiles.
#     3. Reproject the VRT to the specified destination spatial reference system (SRS)
#        using a PROJ pipeline that includes a vertical datum transformation to ellipsoid heights.
#     4. Clip the reprojected DEM to the AOI.

#     Parameters
#     ----------
#     dataset : str
#         Name of the DEM dataset to load. Supported datasets include 'cop30',
#         'nasadem', '3dep_10m', and '3dep_1m'
#     aoi : gpd.GeoDataFrame
#         Area of interest to clip the DEM to.
#     dstSRS : str
#         Destination spatial reference system. Supported = 'EPSG:7912', 'EPSG:9989'

#     Returns
#     -------
#     xr.DataArray
#         DataArray containing the DEM with ellipsoid heights.
#     """
#     print('TODO')
