"""
Microsoft Planetary Computer

https://planetarycomputer.microsoft.com/docs/quickstarts/reading-stac/

- https://planetarycomputer.microsoft.com/dataset/esa-worldcover
- https://planetarycomputer.microsoft.com/dataset/cop-dem-glo-30
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

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
    classmap: dict[int, Any] = field(
        default_factory=lambda: {
            10: {"hex": "#006400", "description": "Tree cover"},
            20: {"hex": "#FFBB22", "description": "Shrubland"},
            30: {"hex": "#FFFF4C", "description": "Grassland"},
            40: {"hex": "#F096FF", "description": "Cropland"},
            50: {"hex": "#FA0000", "description": "Built-up"},
            60: {"hex": "#B4B4B4", "description": "Bare / sparse vegetation"},
            70: {"hex": "#F0F0F0", "description": "Snow and ice"},
            80: {"hex": "#0064C8", "description": "Permanent water bodies"},
            90: {"hex": "#0096A0", "description": "Herbaceous wetland"},
            95: {"hex": "#00CF75", "description": "Mangroves"},
            100: {"hex": "#FAE6A0", "description": "Moss and lichen"},
        }
    )
