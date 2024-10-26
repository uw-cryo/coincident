"""
Read ICESat-2 & GEDI using Sliderule
"""

from __future__ import annotations

from typing import Any

import geopandas as gpd
import pandas as pd
from sliderule import gedi, icesat2, toregion

# TODO: make these configurable?
default_params = {
    "ats": 20.0,
    "cnt": 5,
    "len": 40.0,
    "res": 20.0,
    "maxi": 6,
}


def load_gedi(
    aoi: gpd.GeoDataFrame,
    start_datetime: pd.Timestamp,
    end_datetime: pd.Timestamp,
    quality_flags: bool = True,
    **kwargs: Any,
) -> gpd.GeoDataFrame:
    # NOTE: allow passing private cluster?
    gedi.init("slideruleearth.io")
    params = default_params.copy()
    if quality_flags:
        params["degrade_flag"] = 0
        params["l2_quality_flag"] = 1
    params["poly"] = toregion(aoi)["poly"]
    params["t0"] = start_datetime.tz_localize(None).isoformat()
    params["t1"] = end_datetime.tz_localize(None).isoformat()
    params.update(kwargs)
    return gedi.gedi02ap(params)


def load_icesat2(
    aoi: gpd.GeoDataFrame,
    start_datetime: pd.Timestamp,
    end_datetime: pd.Timestamp,
    **kwargs: Any,
) -> gpd.GeoDataFrame:
    # NOTE: allow passing private cluster?
    # NOTE: deal with kwargs
    icesat2.init("slideruleearth.io")
    params = default_params.copy()
    params["srt"] = icesat2.SRT_LAND
    params["cnf"] = icesat2.CNF_SURFACE_HIGH
    params["poly"] = toregion(aoi)["poly"]
    params["t0"] = start_datetime.tz_localize(None).isoformat()
    params["t1"] = end_datetime.tz_localize(None).isoformat()
    params.update(kwargs)
    return icesat2.atl06p(params)
