"""
NEON Discrete LiDAR Datasets

https://data.neonscience.org/data-products/DP3.30024.001
"""

from __future__ import annotations

from dataclasses import dataclass

from coincident.datasets.general import Dataset


@dataclass
class NEON(Dataset):
    """Essential metadata for NEON LiDAR"""

    has_stac_api: bool = False
    search: str = "http://data.neonscience.org/api/v0/sites"
    start: str = "2005-06-01"  # Or 2013-06-01?
    end: str | None = None
    type: str = "lidar"
    alias: str = "neon"
    provider: str = "neon"
