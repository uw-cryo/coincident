"""
Wrapper for various data provider search APIs

Results returned as GeoDataFrame with consistent STAC column headings and best dtypes
"""

from __future__ import annotations

from coincident.search.main import cascading_search, search
from coincident.search.stac import to_pystac_items

__all__ = ["search", "cascading_search", "to_pystac_items"]
