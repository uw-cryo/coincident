"""
Wrapper for various data provider search APIs

Results returned as GeoDataFrame with consistent STAC column headings and best dtypes
"""

from __future__ import annotations

from coincident.search.main import search
from coincident.search.stac import to_pystac_items

__all__ = ["search", "to_pystac_items"]
