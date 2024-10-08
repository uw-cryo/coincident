"""
NASA Datasets

https://github.com/nasa/cmr-stac

NOTE: Different search endpoints dependng on DAAC
"""

# from pydantic.dataclasses import dataclass
from __future__ import annotations

from dataclasses import dataclass, field

from coincident.datasets.general import Dataset


@dataclass
class ICESat2(Dataset):
    """Essential metadata for ICESat-2 Alitmeter"""

    # https://cmr.earthdata.nasa.gov/stac/NSIDC_ECS/collections/ATL03_006
    has_stac_api: bool = True
    search: str = "https://cmr.earthdata.nasa.gov/stac/NSIDC_ECS"
    start: str | None = "2018-10-13"
    type: str = "altimeter"
    alias: str = "icesat-2"
    collections: list[str] = field(
        default_factory=lambda: ["ATL03_006"]
    )  # ATL08_006 etc.
    provider: str = "nasa"


@dataclass
class GEDI(Dataset):
    """Essential metadata for GEDI Altimeter"""

    # NOTE: parse temporal & bbox from collection metadata?
    # https://cmr.earthdata.nasa.gov/stac/LPCLOUD/collections/GEDI02_A_002
    has_stac_api: bool = True
    search: str = "https://cmr.earthdata.nasa.gov/stac/LPCLOUD"
    start: str = "2019-04-04"
    # https://www.earthdata.nasa.gov/news/nasa-announces-pause-gedi-mission
    end: str = "2023-03-17"
    type: str = "altimeter"
    alias: str = "gedi"
    collections: list[str] = field(default_factory=lambda: ["GEDI02_A_002"])
    provider: str = "nasa"
