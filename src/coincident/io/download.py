"""
Convenience functions for downloading STAC items via stac-asset
"""

from __future__ import annotations

import os
from pathlib import Path

import stac_asset
from pystac import Item


async def download_item(
    item: Item,
    path: str = "/tmp",
    config: str | None = None,
) -> Item:
    """
    Downloads a STAC item to a specified local path.

    Parameters
    ----------
    item : pystac.Item
        The STAC item to be downloaded.
    path : str, optional
        The local directory path where the item will be downloaded. Default is "/tmp".
    config : str, optional
        If config=='maxar', a MAXAR-API-KEY HTTP Header is used for authentication.

    Returns
    -------
    pystac.Item
        The downloaded STAC item.

    Examples
    --------
    >>> localitem = asyncio.run(download_item(item, config=MAXAR_CONFIG))
    """
    posixpath = Path(path)

    if config == "maxar":
        config = stac_asset.Config(
            http_headers={"MAXAR-API-KEY": os.environ.get("MAXAR_API_KEY")}
        )
    # NOTE: prevent in-place modification of remote hrefs with item.clone()
    await stac_asset.download_item(item.clone(), posixpath, config=config)
