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

from coincident.io import download, sliderule, xarray

__all__ = ["download", "sliderule", "xarray"]
