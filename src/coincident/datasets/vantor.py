"""
Vantor (formerly Maxar) VHR stereo imagery

https://www.vantor.com
"""

# from pydantic.dataclasses import dataclass, Field # type: ignore[attr-defined]
from __future__ import annotations

from dataclasses import dataclass, field

from coincident.datasets.general import Dataset


@dataclass
class Stereo(Dataset):
    """Essential metadata for Vantor In Track Stereo"""

    alias: str = "vantor"
    has_stac_api: bool = True
    collections: list[str] = field(
        default_factory=lambda: ["wv01", "wv02", "wv03-vnir", "ge01"]
    )
    search: str = "https://api.maxar.com/discovery/v1/search"
    start: str = "2007-01-01"
    end: str | None = None
    type: str = "stereo"
    # Unique to Vantor
    area_based_calc: bool = False
    provider: str = "vantor"
