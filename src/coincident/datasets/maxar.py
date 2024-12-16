"""
Maxar VHR stereo imagery

https://www.maxar.com/maxar-intelligence/products/satellite-imagery
"""

# from pydantic.dataclasses import dataclass, Field # type: ignore[attr-defined]
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from coincident.datasets.general import Dataset


class Collection(str, Enum):
    """wv01,wv02,wv03-vnir,ge01"""

    wv01 = "wv01"
    wv02 = "wv02"
    wv03_vnir = "wv03-vnir"
    ge01 = "ge01"


@dataclass
class Stereo(Dataset):
    """Essential metadata for Maxar In Track Stereo"""

    alias: str = "maxar"
    has_stac_api: bool = True
    collections: list[Collection] = field(
        default_factory=lambda: ["wv01", "wv02", "wv03-vnir", "ge01"]
    )  # type: ignore[assignment]
    search: str = "https://api.maxar.com/discovery/v1/search"
    start: str = "2007-01-01"
    end: str | None = None
    type: str = "stereo"
    # Unique to Maxar
    area_based_calc: bool = False
    provider: str = "maxar"
