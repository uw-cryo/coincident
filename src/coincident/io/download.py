"""
Convenience functions for downloading STAC items via stac-asset
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import stac_asset
from pystac import Item


def _set_config(
    item: Item,
    config: dict[str, Any] | None = None,
) -> stac_asset.Config:
    """Set stac_asset config. Automatically add Maxar API key if needed"""
    if item.properties.get("constellation") == "maxar" and config is None:
        config = stac_asset.Config(
            http_headers={"MAXAR-API-KEY": os.environ.get("MAXAR_API_KEY")}
        )
    else:
        config = stac_asset.Config(**config)

    return config


async def read_href(
    item: Item,
    asset: str,
    config: dict[str, Any] | None = None,
) -> Item:
    """
    Open and read a STAC asset href into local memory

    Parameters
    ----------
    item : pystac.Item
        The STAC item to be downloaded.
    asset : str
        The asset name to open (e.g. "cloud-cover" for maxar).
    config : dict, optional
        dictionary of options for :class:`~stac_asset.Config`

    Returns
    -------
    bytes
        Bytes read from file at corresponding href

    Examples
    --------
    Read cloud-cover MultiPolygon estimate from Maxar API
    >>> bytes = asyncio.run(read_href(item, "cloud-cover" config=MAXAR_CONFIG))
    """
    href = item.assets.get(asset).href

    config = _set_config(item, config)

    return await stac_asset.read_href(href, config=config)


async def download_item(
    item: Item,
    path: str = "/tmp",
    config: dict[str, Any] | None = None,
) -> Item:
    """
    Downloads a STAC item to a specified local path.

    Parameters
    ----------
    item : pystac.Item
        The STAC item to be downloaded.
    path : str, optional
        The local directory path where the item will be downloaded. Default is "/tmp".
    config : dict, optional
        dictionary of options for :doc:`stac_asset.Config <stac_asset:stac_asset.Config>`

    Returns
    -------
    pystac.Item
        The downloaded STAC item.

    Examples
    --------
    Download all assets for given item to /tmp directory
    >>> localitem = asyncio.run(download_item(item, config=MAXAR_CONFIG))
    """
    posixpath = Path(path)

    config = _set_config(item, config)

    # NOTE: prevent in-place modification of remote hrefs with item.clone()
    return await stac_asset.download_item(item.clone(), posixpath, config=config)
