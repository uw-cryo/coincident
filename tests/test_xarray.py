from __future__ import annotations

import pytest
import xarray as xr

import coincident
from coincident.io.xarray import to_dataset

network = pytest.mark.network


@network
def test_to_dataset_with_cop30_lazy(aoi):
    """Test `to_dataset` functionality with COP30 dataset."""
    gf_cop30 = coincident.search.search(
        dataset="cop30",
        intersects=aoi,
    )
    ds = to_dataset(
        gf_cop30,
        aoi=aoi,
        resolution=0.1,  # ~1km
    )
    assert isinstance(ds, xr.Dataset), "Expected output to be an xarray Dataset."
    assert "data" in ds.data_vars, "Expected 'data' variable in the Dataset."


@network
def test_to_dataset_with_worldcover_lazy(aoi):
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
    )
    assert isinstance(ds, xr.Dataset), "Expected output to be an xarray Dataset."
    assert "map" in ds.data_vars, "Expected 'map' variable in the Dataset."
