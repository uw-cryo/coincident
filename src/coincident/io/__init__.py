"""
coincident: Search and analysis of STV Precursor Coincident Datasets

I/O module is organized into files by libraries used to read data. In
general, a datatype from that library is returned:

download.*   Read bytes into local memory or copy remote files to local disk
xarray.*     Client-side reading of Xarray DataArrays or Datasets
geopandas.*  Client-side reading of Geopandas GeoDataFrames
sliderule.*  Server-side sliderule processing. Return GeoDataFrames or directly write to GeoParquet file
"""

from __future__ import annotations

import os
import warnings

if not os.environ.get("EARTHDATA_TOKEN"):
    warnings.warn(
        "The 'EARTHDATA_TOKEN' environment variable is not set. "
        "Some download functionality may not work. "
        "See https://urs.earthdata.nasa.gov/documentation/for_users/user_token to obtain a token.",
        UserWarning,
        stacklevel=2,
    )

from coincident.io import download, gdal, proj, sliderule, xarray

__all__ = ["download", "gdal", "proj", "sliderule", "xarray"]
