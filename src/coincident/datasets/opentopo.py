"""
OpenTopography Datasets (NOAA, NCALM)

NOAA: https://www.opentopography.org/
NCALM: https://www.ncalm.org/
"""

from __future__ import annotations

from dataclasses import dataclass

from coincident.datasets.general import Dataset


@dataclass
class NOAA(Dataset):
    """Essential metadata for OpenTopography NOAA LiDAR"""

    has_stac_api: bool = False
    search: str = "https://portal.opentopography.org/API/otCatalog?productFormat=PointCloud&polygon={coords}&detail=true&outputFormat=json&include_federated=true"
    start: str = "1996-10-09"
    end: str | None = None
    type: str = "lidar"
    alias: str = "noaa"
    provider: str = "opentopography"


@dataclass
class NCALM(Dataset):
    """Essential metadata for OpenTopography NCALM LiDAR"""

    has_stac_api: bool = False
    search: str = "https://portal.opentopography.org/API/otCatalog?productFormat=PointCloud&polygon={coords}&detail=true&outputFormat=json&include_federated=false"
    start: str = "2003-05-15"
    end: str | None = None
    type: str = "lidar"
    alias: str = "ncalm"
    provider: str = "opentopography"
