"""
USGS 3DEP LiDAR

https://www.usgs.gov/3d-elevation-program
"""

# from pydantic.dataclasses import dataclass
from __future__ import annotations

from dataclasses import dataclass

from coincident.datasets.general import Dataset


@dataclass
class ThreeDEP(Dataset):
    """Essential metadata for USGS 3DEP"""

    has_stac_api: bool = False
    search: str = (
        # "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/metadata/WESM.gpkg"
        "s3://prd-tnm/StagedProducts/Elevation/metadata/WESM.gpkg"
    )
    start: str = "2000-12-01"
    end: str | None = None
    spatial_ref: str | None = "EPSG:6318+5703"  # NAD83(2011) + NAVD88 height
    type: str = "lidar"
    alias: str = "3dep"
    provider: str = "usgs"
    provider_docs: str | None = (
        "https://www.usgs.gov/3d-elevation-program/about-3dep-products-services"
    )


@dataclass
class ThreeDEP_1m(Dataset):
    """3DEP 1-meter seamless"""

    has_stac_api: bool = False
    search: str = "https://github.com/uw-cryo/stac-3dep-1m/releases/download/v0.2/catalogv2.parquet"
    start: str = "2000-12-01"
    end: str | None = None
    spatial_ref: str | None = (
        "UTM: NAD83(2011) + NAVD88"  # Note this is not a standard ID, b/c source data has many UTM zones
    )
    type: str = "lidar"
    alias: str = "3dep-1m"
    provider: str = "usgs"
    provider_docs: str | None = (
        "https://portal.opentopography.org/raster?opentopoID=OTNED.012021.4269.3"
    )
