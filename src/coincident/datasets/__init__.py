"""
Supported datasets

Files are organized by data provider, each of which can contain multiple datasets

ATL03_v6: https://cmr.earthdata.nasa.gov/search/concepts/C2596864127-NSIDC_CPRD
ATL06_v6: https://cmr.earthdata.nasa.gov/search/concepts/C2564427300-NSIDC_ECS
GEDI_2a_v2: https://cmr.earthdata.nasa.gov/search/concepts/C1908348134-LPDAAC_ECS # did this ID change?!
USGS_3DEP: https://cmr.earthdata.nasa.gov/search/concepts/C2021957295-LPCLOUD
WESM:
"""

from __future__ import annotations

from coincident.datasets import (
    csda,
    maxar,
    nasa,
    neon,
    opentopo,
    planetary_computer,
    usgs,
)
from coincident.datasets.general import Dataset

# Convenience mapping of string aliases to supported dataset classes
_datasets = [
    maxar.Stereo(),
    usgs.ThreeDEP(),
    nasa.ICESat2(),
    nasa.GEDI(),
    nasa.GLiHT(),
    nasa.ABLVIS2_1(),
    nasa.AFLVIS2_1(),
    planetary_computer.COP30(),
    planetary_computer.WorldCover(),
    csda.TDX(),
    opentopo.NOAA(),
    opentopo.NCALM(),
    neon.NEON(),
]

aliases = [x.alias for x in _datasets]
_alias_to_Dataset = dict(zip(aliases, _datasets, strict=False))

__all__ = ["Dataset", "csda", "maxar", "nasa", "opentopo", "planetary_computer", "usgs"]
