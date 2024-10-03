"""
STAC Search Functions
"""

from __future__ import annotations

import contextlib
import warnings
from typing import Any

import geopandas as gpd
import pystac_client

# maxar_platform tries to authenticate upon import, change error to warning
try:
    import maxar_platform  # noqa: F401
except Exception as e:
    message = f"maxar_plaform authentication error: {e}, setting MAXAR_API_KEY environment variable is recommended."
    warnings.warn(message, stacklevel=2)

# MSPC authenticated access has higher limits
# try:
#     import planetary_computer
#     modifier = planetary_computer.sign_inplace
# except ModuleNotFoundError as e:
#     modifier = None
#     error_msg = str(e)
#     warnings.warn(f'install plantary_computer package for API TOKEN handling: {error_msg}')
modifier = None


# NOTE: consolidate to minimal number of columns at this point
def add_columns(gf: gpd.GeoDataFrame, featureCollection: dict[str, Any]) -> None:
    """Add stac_id, href, browse & metadata links to dataframe in-place"""
    first_asset = featureCollection["features"][0]["assets"]
    # Common to all collections?
    gf["stac_id"] = [x["id"] for x in featureCollection["features"]]
    # NOTE: NASA STAC item asset naming seems inconsistent... 'data' or stac_id for direct download link...
    # ICESAt-2 has 'thumbnail_9' up to 28 and then 'data' is 'thumbnail_28', '03/ATL03_20230103090928_02111806_006_02' with
    # NOTE: MAXAR items have: 'browse', 'cloud-cover', 'sample-point-set'
    # NOTE: COP30 'tilejson', 'rendered_preview', 'data'
    assets = [f["assets"] for f in featureCollection["features"]]
    if "data" in first_asset:
        gf["data"] = [a["data"].get("href") for a in assets]
    else:
        # NOTE: Better to parse by 'role' ('data','browse', etc) since asset keys can be anything?
        # ugh, apparently it's not STAC Core (maxar API doesn't use roles)
        with contextlib.suppress(KeyError):
            gf["data"] = [
                v.get("href")
                for A in assets
                for (k, v) in A.items()
                if v["roles"] == ["data"]
            ]
    if "browse" in first_asset:
        gf["browse"] = [a["browse"].get("href") for a in assets]
    if "metadata" in first_asset:
        gf["metadata"] = [
            a["metadata"].get("href") for a in assets
        ]  # Not sure if any link to this? or add 'self' link?


def to_geopandas(
    search_results: pystac_client.item_search.ItemSearch,
) -> gpd.GeoDataFrame:
    """Convert returned from STAC API to geodataframe"""
    featureCollection = search_results.item_collection_as_dict()
    gf = gpd.GeoDataFrame.from_features(featureCollection, crs="epsg:4326")

    # datetime_cols = ['collect_time_end','collect_time_start']  #MAXAR-specific
    # ValueError: cannot supply both a tz and a timezone-naive dtype (i.e. datetime64[ns]): Error while type casting for column 'datetime'
    # gf.astype({'datetime':'datetime64[ns]'}) #All?
    # gf['datetime'] = gpd.pd.to_datetime(gf.datetime)
    # NOTE: assumption = all timestamps UTC, don't worry about timezone from here on out
    gf["datetime"] = gpd.pd.to_datetime(gf.datetime).dt.tz_localize(None)
    gf["dayofyear"] = gf["datetime"].dt.dayofyear
    if "start_datetime" in gf.columns:
        gf["start_datetime"] = gpd.pd.to_datetime(gf.start_datetime).dt.tz_localize(
            None
        )
        gf["end_datetime"] = gpd.pd.to_datetime(gf.end_datetime).dt.tz_localize(None)
    if "collect_time_end" in gf.columns:
        gf["start_datetime"] = gpd.pd.to_datetime(gf.collect_time_start).dt.tz_localize(
            None
        )
        gf["end_datetime"] = gpd.pd.to_datetime(gf.collect_time_end).dt.tz_localize(
            None
        )

    # Add URLs to specific data formats
    add_columns(gf, featureCollection)

    return gf


def search(
    client: pystac_client.client.Client, **kwargs: dict[str, Any] | None
) -> pystac_client.item_search.ItemSearch:
    """Search any STAC API (e.g. https://github.com/nasa/cmr-stac)"""
    return client.search(
        method="GET",  # default POST but NASA-CMR-STAC seems to require GET
        **kwargs,
    )


def configure_maxar_client(area_based_calc: bool = True) -> pystac_client.client.Client:
    """ """
    import maxar_platform.discovery  # automatically checks authentication

    client = maxar_platform.discovery.open_catalog(catalog="imagery")
    # Custom Maxar setting
    client.area_based_calc = area_based_calc

    return client


def configure_stac_client(url: str) -> pystac_client.client.Client:
    """ """
    return pystac_client.Client.open(url=url)


def configure_mspc_client(url: str) -> pystac_client.client.Client:
    """ """
    return pystac_client.Client.open(
        url=url,
        modifier=modifier,  # planetary_computer.sign_inplace?
    )
