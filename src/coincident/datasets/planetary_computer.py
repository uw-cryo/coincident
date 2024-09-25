"""
Microsoft Planetary Computer

https://planetarycomputer.microsoft.com/docs/quickstarts/reading-stac/
"""

from __future__ import annotations

from dataclasses import dataclass, field

from coincident.datasets.general import Dataset


@dataclass
class COP30(Dataset):
    """Essential metadata for Copernicus DEM"""

    alias: str = "cop30"
    has_stac_api: bool = True
    collections: list[str] = field(default_factory=lambda: ["cop-dem-glo-30"])
    search: str = "https://planetarycomputer.microsoft.com/api/stac/v1"
    start: str | None = None  # NOTE: has 'representative' datetime of 2021-04-22
    end: str | None = None
    type: str = "sar"
