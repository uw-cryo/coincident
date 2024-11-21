"""
Maxar VHR stereo imagery

https://www.maxar.com/maxar-intelligence/products/satellite-imagery
"""

# from pydantic.dataclasses import dataclass, Field # type: ignore[attr-defined]
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from stac_asset import Config

import os
import warnings
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

import numpy as np
import pystac
import rasterio
import stac_asset
import xarray as xr

from coincident._utils import depends_on_optional
from coincident.datasets.general import Dataset

MAXAR_CONFIG = stac_asset.Config(
    http_headers={"MAXAR-API-KEY": os.environ.get("MAXAR_API_KEY")}
)
try:
    import matplotlib.pyplot as plt
except ImportError:
    warnings.warn(
        "'matplotlib' package not found. Install for plotting functions: https://matplotlib.org/stable/install/index.html",
        stacklevel=2,
    )


class Collection(str, Enum):
    """wv01,wv02,wv03-vnir,ge01"""

    wv01 = "wv01"
    wv02 = "wv02"
    wv03_vnir = "wv03-vnir"
    ge01 = "ge01"


@dataclass
class Stereo(Dataset):
    """Essential metadata for Maxar In Track Stereo"""

    alias: str = "maxar"
    has_stac_api: bool = True
    collections: list[Collection] = field(
        default_factory=lambda: ["wv01", "wv02", "wv03-vnir", "ge01"]
    )  # type: ignore[assignment]
    search: str = "https://api.maxar.com/discovery/v1/search"
    start: str = "2007-01-01"
    end: str | None = None
    type: str = "stereo"
    # Unique to Maxar
    area_based_calc: bool = False
    provider: str = "maxar"


# NOTE: expose as coincident.search.download
# Or just add documentation for how to do this?
async def download_item(
    item: pystac.Item,
    path: str = "/tmp",
    config: Config = MAXAR_CONFIG,
) -> pystac.Item:
    """localitem = asyncio.run(download_item(item, config=MAXAR_CONFIG))"""
    posixpath = Path(path)
    item = await stac_asset.download_item(item, posixpath, config=config)
    return item  # noqa: RET504


def open_browse(item: pystac.Item, overview_level: int = 0) -> xr.DataArray:
    """
    Open a browse image from a STAC item using the specified overview level.

    Parameters
    ----------
    item : pystac.Item
        The STAC item containing the browse image asset.
    overview_level : int, optional
        The overview level to use when opening the image, by default 0.

    Returns
    -------
    xr.DataArray
        The opened browse image as an xarray DataArray.

    Notes
    -----
    The function uses the `rasterio` engine to open the image and sets the
    `GDAL_DISABLE_READDIR_ON_OPEN` and `GDAL_HTTP_HEADERS` environment variables
    for optimized reading and authentication, respectively.
    """

    # href = item.assets['browse']['href'] # accessed via geopandas
    href = item.assets["browse"].href  # accessed via pystac

    env = rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        GDAL_HTTP_HEADERS=f'MAXAR-API-KEY:{os.environ["MAXAR_API_KEY"]}',
    )
    with env:
        return xr.open_dataarray(
            href,
            engine="rasterio",
            mask_and_scale=False,  # otherwise uint8 -> float32!
            backend_kwargs={"open_kwargs": {"overview_level": overview_level}},
        )


@depends_on_optional("matplotlib")
def plot_browse(item: pystac.Item, overview_level: int = 0) -> None:
    """
    Plots a browse image from a STAC item using Matplotlib.

    Parameters
    ----------
    item : pystac.Item
        The STAC item containing the browse image to be plotted.
    overview_level : int, optional
        The overview level of the browse image to be opened, by default 0.

    Returns
    -------
    None
        This function does not return any value
    """

    da = open_browse(item, overview_level=overview_level)
    mid_lat = da.y[int(da.y.size / 2)].to_numpy()  # PD011

    fig, ax = plt.subplots(figsize=(8, 11))  # pylint: disable=unused-variable
    da.plot.imshow(rgb="band", add_labels=False, ax=ax)
    ax.set_aspect(aspect=1 / np.cos(np.deg2rad(mid_lat)))
    ax.set_title(item.id)
