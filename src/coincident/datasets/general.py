"""
Supported datasets:

ATL03_v6: https://cmr.earthdata.nasa.gov/search/concepts/C2596864127-NSIDC_CPRD
ATL06_v6: https://cmr.earthdata.nasa.gov/search/concepts/C2564427300-NSIDC_ECS
GEDI_2a_v2: https://cmr.earthdata.nasa.gov/search/concepts/C1908348134-LPDAAC_ECS
USGS_3DEP: https://cmr.earthdata.nasa.gov/search/concepts/C2021957295-LPCLOUD
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
    type: str | None = None  # lidar | stereo | altimeter | sar
    provider: str | None = None  # usgs, maxar, nasa, microsoft, csda
    # Pystac client default limit=100, but seems set by API endpoint as well (nasa cmr-stac=25)
    stac_kwargs: dict[str, Any] = field(default_factory=lambda: {"limit": 1000})
