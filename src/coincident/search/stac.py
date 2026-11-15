"""
STAC Search Functions
"""

from __future__ import annotations

from typing import Any

import arro3.core
import cql2
import geopandas as gpd
import pystac
import pystac_client
import rustac
from pystac_client.stac_api_io import StacApiIO

from coincident._utils import VANTOR_API_KEY_HEADER, vantor_api_key


class VantorStacApiIO(StacApiIO):
    """StacApiIO for the Vantor Discovery API

    Adds two Vantor-specific behaviors to standard STAC API requests:

    - CQL2 text filters are converted to CQL2 JSON, since the Discovery API
      does not accept text filters in POST requests
    - Optionally request area-based calculations (e.g. cloud cover within the
      search AOI rather than the entire scene)
      https://developers.maxar.com/docs/discovery/guides/discovery-guide#area-based-calculations
    """

    def __init__(self, area_based_calc: bool = False, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.area_based_calc = area_based_calc

    def request(
        self,
        href: str,
        method: str | None = None,
        headers: dict[str, str] | None = None,
        parameters: dict[str, Any] | None = None,
    ) -> str:
        if parameters:
            parameters = dict(parameters)
            has_aoi = "bbox" in parameters or "intersects" in parameters
            if method == "POST":
                if parameters.get("filter-lang") == "cql2-text":
                    parameters["filter"] = cql2.parse_text(
                        parameters["filter"]
                    ).to_json()
                    parameters["filter-lang"] = "cql2-json"
                if self.area_based_calc and has_aoi:
                    parameters["area-based-calc"] = True
            elif self.area_based_calc and has_aoi:
                parameters["area-based-calc"] = 1
        return str(super().request(href, method, headers, parameters))


def _filter_assets(assets: gpd.GeoDataFrame) -> dict[str, str]:
    """Remove key:None pairs from assets"""
    keep_keys = []
    for k, v in assets.items():
        if v is not None:
            keep_keys.append(k)

    return {key: assets[key] for key in keep_keys}


def to_geopandas(
    collection: list[pystac.Item]
    | pystac.item_collection.ItemCollection
    | arro3.core.Table
    | None,
) -> gpd.GeoDataFrame:
    """
    Convert a STAC ItemCollection to a GeoDataFrame.
    It also adds an additional column 'dayofyear' for convenience.

    Parameters
    ----------
    collection
        The STAC ItemCollection to be converted.

    Returns
    -------
        A GeoDataFrame containing the data from the STAC ItemCollection.

    Raises
    ------
    ValueError
        If the provided ItemCollection is empty.
    """
    # Catch if no items are passed
    if collection is None or (hasattr(collection, "__len__") and len(collection) == 0):
        message = "ItemCollection is empty, cannot convert to GeoDataFrame"
        raise ValueError(message)

    if isinstance(collection, pystac.item_collection.ItemCollection):
        collection = rustac.to_arrow(collection.to_dict())

    # NOTE: I believe all STAC footprints are expected as EPSG:4326 (proj:geometry is optional)
    gf = gpd.GeoDataFrame.from_arrow(collection).set_crs(epsg=4326)

    # Workaround stac-geoparquet limitation https://github.com/stac-utils/stac-geoparquet/issues/82
    gf["assets"] = gf["assets"].apply(_filter_assets)

    # Additional columns for convenience
    # NOTE: these become entries under STAC properties
    gf["dayofyear"] = gf["datetime"].dt.day_of_year

    return gf


def to_pystac_items(gf: gpd.GeoDataFrame) -> list[pystac.Item]:
    """
    Converts a GeoDataFrame to a list of PySTAC Items.

    Parameters
    ----------
    gf
        The GeoDataFrame to be converted.

    Returns
    -------
        A list of PySTAC Items created from the GeoDataFrame.
    """

    feature_collection = rustac.from_arrow(gf.to_arrow())
    item_collection = pystac.item_collection.ItemCollection.from_dict(
        feature_collection
    )
    return item_collection.items


def search(
    client: pystac_client.client.Client, **kwargs: Any
) -> pystac.item_collection.ItemCollection:
    """Search any STAC API (e.g. https://github.com/nasa/cmr-stac)"""
    # NOTE: add logging for kwargs?
    # print(kwargs)
    results = client.search(
        **kwargs,
    )
    return results.item_collection()  # actually run the query


def configure_vantor_client(
    url: str,
    area_based_calc: bool = False,
) -> pystac_client.client.Client:
    """STAC client for the Vantor Discovery API, authenticated with VANTOR_API_KEY"""
    stac_io = VantorStacApiIO(
        area_based_calc=area_based_calc,
        headers={VANTOR_API_KEY_HEADER: vantor_api_key()},
        timeout=60,
        max_retries=5,
    )
    return pystac_client.Client.open(url=url, stac_io=stac_io)


def configure_stac_client(url: str) -> pystac_client.client.Client:
    return pystac_client.Client.open(
        url=url, stac_io=StacApiIO(timeout=60, max_retries=5)
    )
