"""
USGS 3DEP Search via WESM.gpkg
"""

from __future__ import annotations

from typing import Any

import geopandas as gpd


def simplify_footprint(gf: gpd.GeoDataFrame, tolerance: float = 0.1) -> None:
    """complexity of WESM footprints causes slow STAC searches"""
    # parts of a simplified geometry will be no more than `tolerance` distance from the original
    gf = gf.simplify(tolerance)


def stacify_column_names(gf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """match STAC items"""
    name_map = {
        "collect_start": "start_datetime",
        "collect_end": "end_datetime",
    }
    gf = gf.rename(columns=name_map)
    duration = gf.end_datetime - gf.start_datetime
    gf["datetime"] = gf.start_datetime + duration / 2
    gf["duration"] = duration.dt.days

    return gf


def search_gpkg(
    url: str, mask: gpd.GeoDataFrame | None = None, **kwargs: dict[str, Any] | None
) -> gpd.GeoDataFrame:
    """pass a geodataframe or shapely geometry to read portion of full WESM.gpkg"""

    gf = gpd.read_file(
        url,
        mask=mask,
        **kwargs,
        # NOTE: worth additional dependencies for speed?
        # Only faster I think if GPKG is local https://github.com/geopandas/pyogrio/issues/252
        # engine='pyogrio',
        # pyarrow=True,
    )

    # For the purposes of search, just ignore datum (NAD83)
    gf = gf.set_crs("EPSG:4326", allow_override=True)

    return stacify_column_names(gf)
