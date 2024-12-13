from __future__ import annotations

import pytest
import xarray as xr

import coincident
from coincident.io.xarray import plot_esa_worldcover, to_dataset

# Decorate tests requiring internet (slow & flaky)
network = pytest.mark.network


@network
def test_to_dataset_with_cop30(aoi):
    """Test `to_dataset` functionality with COP30 dataset."""
    gf_cop30 = coincident.search.search(
        dataset="cop30",
        intersects=aoi,
    )
    ds = to_dataset(
        gf_cop30,
        aoi=aoi,
        resolution=0.1,  # ~1km
    ).compute()
    assert isinstance(ds, xr.Dataset), "Expected output to be an xarray Dataset."
    assert "data" in ds.data_vars, "Expected 'data' variable in the Dataset."


@network
def test_to_dataset_with_worldcover(aoi):
    """Test `to_dataset` functionality with WorldCover dataset."""
    gf_wc = coincident.search.search(
        dataset="worldcover",
        intersects=aoi,
        datetime=["2020"],
    )
    ds = to_dataset(
        gf_wc,
        bands=["map"],
        aoi=aoi,
        resolution=0.1,  # ~1km
    ).compute()
    assert isinstance(ds, xr.Dataset), "Expected output to be an xarray Dataset."
    assert "map" in ds.data_vars, "Expected 'map' variable in the Dataset."


@network
def test_plot_esa_worldcover_valid(aoi):
    """Test `plot_esa_worldcover` with valid WorldCover dataset."""
    import matplotlib.collections

    gf_wc = coincident.search.search(
        dataset="worldcover",
        intersects=aoi,
        datetime=["2021"],
    )
    ds = to_dataset(
        gf_wc,
        bands=["map"],
        aoi=aoi,
        resolution=0.1,  # ~1km
    ).compute()
    ds = ds.rename(map="landcover")
    ax = plot_esa_worldcover(ds)
    assert ax is not None, "Expected a valid Matplotlib Axes object."
    # https://matplotlib.org/stable/users/prev_whats_new/whats_new_3.4.0.html
    # https://github.com/matplotlib/matplotlib/blob/main/lib/matplotlib/tests/test_contour.py#L146
    assert any(
        isinstance(c, matplotlib.collections.QuadMesh) for c in ax.get_children()
    ), "Expected at least one pcolormesh object in the plot."
