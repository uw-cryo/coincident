"""
Datasets for:
- NGS survey mark control monuments
- Permanent GNSS control position(s)
"""

from __future__ import annotations

from dataclasses import dataclass

from coincident.datasets.general import Dataset


@dataclass
class NgsStations(Dataset):
    """Essential metadata for NGS Survey Stations."""

    has_stac_api: bool = False
    # The search function is not a simple URL, so we point to the API docs
    search: str | None = "https://geodesy.noaa.gov/api/nde/bounds"
    start: str | None = "1983-01-01"  # TODO: double check this
    end: str | None = None
    type: str = "control"
    alias: str = "ngs-stations"

    provider: str = "ngs"


@dataclass
class GageStations(Dataset):
    """Essential metadata for GAGE GPS Survey Stations."""

    has_stac_api: bool = False
    # The search function is not a simple URL, so we point to the API docs
    search: str | None = "https://web-services.unavco.org/gps/metadata/sites/v1"
    start: str | None = "1986-01-01"  # TODO: double check this
    end: str | None = None
    type: str = "control"
    alias: str = "gage-stations"

    provider: str = "gage"
