"""
Basic dataset metadata structure
"""

# from pydantic.dataclasses import dataclass, Field # type: ignore[attr-defined]
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


# NOTE: default to None for all of these?
@dataclass
class Dataset:
    """Essential metadata for supported datasets"""

    alias: str | None = None  # nickname
    has_stac_api: bool | None = None  # whether or not its a STAC API
    collections: list[str] = field(
        default_factory=list
    )  # STAC collection names of specific datasets
    search: str | None = None  # search API endpoint
    start: str | None = None  # first acquisition date
    end: str | None = None  # last acquisition date (or None if ongoing)
    spatial_ref: str | None = (
        None  # override spatial reference system (EPSG code or WKT2/PROJJSON string)
    )
    type: str | None = None  # lidar | stereo | altimeter | sar
    provider: str | None = None  # usgs, maxar, nasa, microsoft, csda
    provider_docs: str | None = None  # URL to provider documentation
    # Pystac client default limit=100, but seems set by API endpoint as well (nasa cmr-stac=25)
    stac_kwargs: dict[str, Any] = field(default_factory=lambda: {"limit": 1000})
