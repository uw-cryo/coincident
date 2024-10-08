"""
STAC Search Functions
"""

from __future__ import annotations

from typing import Any

import geopandas as gpd
import planetary_computer
import pystac
import pystac_client
import stac_geoparquet

# Only needed if planetary_computer SDK is optional dependency...
# try:
#     import planetary_computer
#     modifier = planetary_computer.sign_inplace
# except ModuleNotFoundError as e:
#     modifier = None
#     message = f'{str(e)}: for authenticated requests install plantary_computer package (https://github.com/microsoft/planetary-computer-sdk-for-python)'
#     warnings.warn(message)

# NOTE: add function to keep only single thumbnail in assets from NASA CMR searches?


def to_geopandas(
    collection: pystac.item_collection.ItemCollection,
) -> gpd.GeoDataFrame:
    """Convert returned from STAC API to geodataframe via arrow"""
    record_batch_reader = stac_geoparquet.arrow.parse_stac_items_to_arrow(collection)
    gf = gpd.GeoDataFrame.from_arrow(record_batch_reader)  # doesn't keep arrow dtypes

    # Additional columns for convenience
    # NOTE: these become entries under STAC properties
    gf["dayofyear"] = gf["datetime"].dt.dayofyear

    return gf


def to_pystac_items(gf: gpd.GeoDataFrame) -> list[pystac.Item]:
    """
    Converts a :class:`~geopandas.GeoDataFrame` to a list of PySTAC Items.

    Parameters
    ----------
    gf : GeoDataFrame
        The GeoDataFrame to be converted.

    Returns
    -------
    list[pystac.Item]
        A list of PySTAC Items created from the GeoDataFrame.
    """

    batch = stac_geoparquet.arrow.stac_table_to_items(gf.to_arrow())

    return [pystac.Item.from_dict(x) for x in batch]


def search(
    client: pystac_client.client.Client, **kwargs: dict[str, Any] | None
) -> pystac_client.item_search.ItemSearch:
    """Search any STAC API (e.g. https://github.com/nasa/cmr-stac)"""
    results = client.search(
        **kwargs,
    )
    return results.item_collection()  # actually run the query


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
    return pystac_client.Client.open(url=url, modifier=planetary_computer.sign_inplace)
