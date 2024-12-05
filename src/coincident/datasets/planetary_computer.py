"""
Microsoft Planetary Computer

https://planetarycomputer.microsoft.com/docs/quickstarts/reading-stac/

- https://planetarycomputer.microsoft.com/dataset/esa-worldcover
- https://planetarycomputer.microsoft.com/dataset/cop-dem-glo-30
"""

from __future__ import annotations

from dataclasses import dataclass, field

from coincident.datasets.general import Dataset

STACAPI = "https://planetarycomputer.microsoft.com/api/stac/v1"


@dataclass
class COP30(Dataset):
    """Essential metadata and data access for Copernicus DEM"""

    alias: str = "cop30"
    has_stac_api: bool = True
    collections: list[str] = field(default_factory=lambda: ["cop-dem-glo-30"])
    search: str = STACAPI
    start: str | None = None  # Copernicus DEM has 'representative' datetime: 2021-04-22
    end: str | None = None
    type: str = "dem"
    provider: str = "microsoft"


@dataclass
class WorldCover(Dataset):
    """Essential metadata for ESA WorldCover"""

    alias: str = "worldcover"
    has_stac_api: bool = True
    collections: list[str] = field(default_factory=lambda: ["esa-worldcover"])
    search: str = STACAPI
    start: str = "2020-01-01"
    end: str = "2021-12-31"
    type: str = "lulc"
    provider: str = "microsoft"
