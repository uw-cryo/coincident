"""
STAC Search Functions
"""

from __future__ import annotations

import warnings
from typing import Any

import geopandas as gpd
import planetary_computer
import pystac
import pystac_client
import stac_geoparquet

# Any import that requires auth to work will be optional
try:
    import maxar_platform
except ImportError:
    message = "'maxar-platform' package not found. Install for maxar functionality: https://pypi.org/project/maxar-platform/"
    warnings.warn(message, stacklevel=2)

try:
    import maxar_platform.discovery
except maxar_platform.session.NoSessionCredentials:
    message = "Unable to authenticate with Maxar API. Please set MAXAR_API_KEY environment variable."
    warnings.warn(message, stacklevel=2)


def to_geopandas(
    collection: pystac.item_collection.ItemCollection,
) -> gpd.GeoDataFrame:
    """
    Convert a STAC ItemCollection to a GeoDataFrame.
    This function converts a given STAC ItemCollection to a GeoDataFrame using the
    `stac_geoparquet.arrow.parse_stac_items_to_arrow` method. It also adds an additional
    column 'dayofyear' for convenience.

    Parameters
    ----------
    collection : pystac.item_collection.ItemCollection
        The STAC ItemCollection to be converted.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the data from the STAC ItemCollection.

    Raises
    ------
    ValueError
        If the provided ItemCollection is empty.
    """
    # Catch if no items are passed
    if len(collection) == 0:
        message = "ItemCollection is empty, cannot convert to GeoDataFrame"
        raise ValueError(message)

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
    # NOTE: add logging for kwargs?
    # print(kwargs)
    results = client.search(
        **kwargs,
    )
    return results.item_collection()  # actually run the query


def configure_maxar_client(area_based_calc: bool = True) -> pystac_client.client.Client:
    # automatically checks authentication

    client = maxar_platform.discovery.open_catalog(catalog="imagery")
    # Custom Maxar setting
    client.area_based_calc = area_based_calc

    return client


def configure_stac_client(url: str) -> pystac_client.client.Client:
    return pystac_client.Client.open(url=url)


def configure_mspc_client(url: str) -> pystac_client.client.Client:
    return pystac_client.Client.open(url=url, modifier=planetary_computer.sign_inplace)
