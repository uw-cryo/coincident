from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

import coincident

network = pytest.mark.network
try:
    import maxar_platform.discovery  # noqa: F401

    not_authenticated = False
except:  # noqa: E722
    not_authenticated = True
maxar_authenticated = pytest.mark.skipif(
    not_authenticated, reason="Not authenticated with Maxar API"
)


@network
@maxar_authenticated
@pytest.mark.filterwarnings("ignore:the actual content type does not match")
@pytest.mark.skip(reason="temporarily skip to see if CI failing due to async code")
def test_download_maxar_browse():
    gf = coincident.search.search(dataset="maxar", ids=["102001008EC5AC00"])
    item = coincident.search.stac.to_pystac_items(gf)[0]
    asyncio.run(coincident.io.download.download_item(item, config="maxar"))
    assert Path("/tmp/102001008EC5AC00.browse.tif").exists()
    assert Path("/tmp/102001008EC5AC00.json").exists()


# @network
# def test_download_planetary_computer():
#     print('todo: small cop30 asset?')
