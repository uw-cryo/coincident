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


# LVIS
# NOTE: different LVIS catalogs are scattered and only some have STAC support...
# Accessing LVIS data through an LP DAAC via earthaccess would include the below two
# datasets, making them mostly-redundant. Can't find STAC catalogs for G-LiHT either
# but the ASF STAC endpoint is promising and will allow us to easily implement
# some Sentinel-1, ALOS, etc access
# https://cmr.earthdata.nasa.gov/stac/ASF
# =======
@dataclass
class ABLVIS2_1(Dataset):
    """
    Essential metadata for ABoVE LVIS L2 Geolocated Surface Elevation Product, Version 1

    Search extent limited to bbox [-158, 48, -104, 72] and 2017-06-29 to 2017-07-17

    This data set contains surface elevation data over Alaska and Western Canada measured by the
    NASA Land, Vegetation, and Ice Sensor (LVIS), an airborne lidar scanning laser altimeter.
    The data were collected as part of NASA's Terrestrial Ecology Program campaign,
    the Arctic-Boreal Vulnerability Experiment (ABoVE).
    """

    has_stac_api: bool = True
    search: str = "https://cmr.earthdata.nasa.gov/stac/NSIDC_ECS"
    start: str = "2017-06-29"
    end: str = "2017-07-17"
    type: str = "lidar"
    alias: str = "ablvis2_1"
    collections: list[str] = field(default_factory=lambda: ["ABLVIS2_1"])
    provider: str = "nasa"


@dataclass
class AFLVIS2_1(Dataset):
    """
    Essential metadata for AfriSAR LVIS L2 Geolocated Surface Elevation Product, Version 1

    Search extent limited to bbox [8, -2, 12, 1] and 2016-02-20 to 2016-03-08

    This data set contains surface elevation data over Gabon, Africa.
    The measurements were taken by the NASA Land, Vegetation, and Ice Sensor (LVIS),
    an airborne lidar scanning laser altimeter.
    The data were collected as part of a NASA campaign,
    in collaboration with the European Space Agency (ESA) mission AfriSAR.
    """

    has_stac_api: bool = True
    search: str = "https://cmr.earthdata.nasa.gov/stac/NSIDC_ECS"
    start: str = "2016-02-20"
    end: str = "2016-03-08"
    type: str = "lidar"
    alias: str = "aflvis2_1"
    collections: list[str] = field(default_factory=lambda: ["AFLVIS2_1"])
    provider: str = "nasa"
