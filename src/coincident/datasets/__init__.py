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

from coincident.datasets import csda, maxar, nasa, planetary_computer, usgs
from coincident.datasets.general import Dataset

# Convenience mapping of string aliases to supported dataset classes
_datasets = [
    maxar.Stereo(),
    usgs.ThreeDEP(),
    nasa.ICESat2(),
    nasa.GEDI(),
    planetary_computer.COP30(),
    planetary_computer.WorldCover(),
    csda.TDX(),
]

aliases = [x.alias for x in _datasets]
_alias_to_Dataset = dict(zip(aliases, _datasets, strict=False))

__all__ = ["Dataset", "usgs", "maxar", "nasa", "planetary_computer", "csda"]
