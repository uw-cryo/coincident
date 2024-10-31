"""
NASA Commercial SmallSat Data Acquisition Program (CSDA)

https://www.earthdata.nasa.gov/esds/csda

TSX/TDX = https://csdap.earthdata.nasa.gov/stac/collections/airbus
# TODO: compare this archive against
# https://api.oneatlas.airbus.com/api-catalog-v2/radar/overview/index.html

Notes
https://github.com/stac-extensions/sar?tab=readme-ov-file#product-type
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from coincident.datasets.general import Dataset


@dataclass
class TDX(Dataset):
    """Essential metadata for TanDEM-X / TerraSAR-X
    platform =
     - PAZ (PAX)
     - TerraSAR-X (TSX-1)
     - TanDEM-X (TDX-1)

    sar:product_type =
    - Coregistered Single Look Slant Range Complex (COSSC)
    - Geocoded Ellipsoid Corrected (GEC)
    - Single Look Slant Range Complex (SSC)
    - Enhanced Ellipsoid Corrected (EEC)

    sar:instrument_mode =
     - Staring Spotlight (ST)
     - High Resolution Spotlight (HS)
     - Spotlight (SL)
     - Stripmap (SM)
     - ScanSAR (SC)
     - Wide ScanSAR (WSC)

    Other properties to filter on:
    sar:polarizations = ['HH']
    orbitDirection = ASCENDING / DESCENDING
    """

    alias: str = "tdx"
    has_stac_api: bool = True
    collections: list[str] = field(default_factory=lambda: ["airbus"])
    search: str = "https://csdap.earthdata.nasa.gov/stac"
    start: str = "2007-07-01"  # TODO: check this
    end: str | None = None
    type: str = "sar"
    provider: str = "csda"
    stac_kwargs: dict[str, Any] = field(
        default_factory=lambda: {
            "limit": 1000,
            "filter": {
                "op": "and",
                "args": [
                    # exclude PAZ, only SSC products
                    {
                        "op": "in",
                        "args": [{"property": "platform"}, ["TDX-1", "TSX-1"]],
                    },
                    {"op": "=", "args": [{"property": "sar:product_type"}, "SSC"]},
                ],
            },
        }
    )
