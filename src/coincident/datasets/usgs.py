"""
USGS 3DEP LiDAR

https://www.usgs.gov/3d-elevation-program
"""

# from pydantic.dataclasses import dataclass
from __future__ import annotations

from dataclasses import dataclass

from coincident.datasets.general import Dataset


@dataclass
class ThreeDEP(Dataset):
    """Essential metadata for USGS 3DEP"""

    has_stac_api: bool = False
    search: str = (
        # "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/metadata/WESM.gpkg"
        "s3://prd-tnm/StagedProducts/Elevation/metadata/WESM.gpkg"
    )
    start: str = "2000-12-01"
    end: str | None = None
    type: str = "lidar"
    alias: str = "3dep"
    provider: str = "usgs"
