"""
NEON LiDAR data

https://data.neonscience.org/data-products/DP3.30024.001
NEON (National Ecological Observatory Network). Elevation - LiDAR (DP3.30024.001), provisional data. Dataset accessed from https://data.neonscience.org/data-products/DP3.30024.001 on September 11, 2024. Data archived at [your DOI].
"""

# from pydantic.dataclasses import dataclass
from __future__ import annotations

from dataclasses import dataclass

from coincident.datasets.general import Dataset


@dataclass
class NEON(Dataset):
    """Essential metadata for NEON"""

    has_stac_api: bool = False
    search: str = "https://data.neonscience.org/data-api"
    start: str = "2013-06-01"
    end: str = "2023-08-01"
    type: str = "lidar"
    alias: str = "neon"
