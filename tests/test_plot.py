from __future__ import annotations

import pytest
from matplotlib.collections import QuadMesh

import coincident

network = pytest.mark.network


@network
def test_plot_esa_worldcover_valid(aoi):
    """Test `plot_esa_worldcover` with valid WorldCover dataset."""

    gf_wc = coincident.search.search(
        dataset="worldcover",
        intersects=aoi,
        datetime=["2021"],
    )
    ds = coincident.io.xarray.to_dataset(
        gf_wc,
        bands=["map"],
        aoi=aoi,
        resolution=0.1,  # ~1km
    ).compute()
    ds = ds.rename(map="landcover")
    ax = coincident.plot.plot_esa_worldcover(ds)
    assert ax is not None, "Expected a valid Matplotlib Axes object."
    # https://matplotlib.org/stable/users/prev_whats_new/whats_new_3.4.0.html
    # https://github.com/matplotlib/matplotlib/blob/main/lib/matplotlib/tests/test_contour.py#L146
    assert any(
        isinstance(c, QuadMesh) for c in ax.get_children()
    ), "Expected at least one pcolormesh object in the plot."
