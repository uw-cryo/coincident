from __future__ import annotations

import pytest
import xarray as xr

import coincident
from coincident.io.xarray import plot_esa_worldcover, to_dataset

# These tests are based on the workflow from https://coincident.readthedocs.io/en/latest/examples/contextual_data.html
# I couldn't think of a way to generalized these tests since the
# xarray funcs are pretty tailored to COP30 and ESA WC

# mypy doesn't like the resolution argument in to_dataset
# Argument "resolution" to "to_dataset" has incompatible type "float"; expected "dict[str, Any]"  [arg-type]

@pytest.fixture
def search_aoi():
    """Fixture to load a simplified AOI from realistic lidar data."""
    workunit = "CO_WestCentral_2019"
    df_wesm = coincident.search.wesm.read_wesm_csv()
    gf_lidar = coincident.search.wesm.load_by_fid(
        df_wesm[df_wesm.workunit == workunit].index
    )
    return gf_lidar.simplify(0.01)


@pytest.mark.network
def test_to_dataset_with_cop30(search_aoi):
    """Test `to_dataset` functionality with COP30 dataset."""
    gf_cop30 = coincident.search.search(
        dataset="cop30",
        intersects=search_aoi,
    )
    ds = to_dataset(
        gf_cop30,
        aoi=search_aoi,
        resolution=0.00081,  # ~90m
        mask=True,
    ).compute()
    assert isinstance(ds, xr.Dataset), "Expected output to be an xarray Dataset."
    assert "data" in ds.data_vars, "Expected 'data' variable in the Dataset."


@pytest.mark.network
def test_to_dataset_with_worldcover(search_aoi):
    """Test `to_dataset` functionality with WorldCover dataset."""
    gf_wc = coincident.search.search(
        dataset="worldcover",
        intersects=search_aoi,
        datetime=["2020"],
    )
    ds = to_dataset(
        gf_wc,
        bands=["map"],
        aoi=search_aoi,
        resolution=0.00081,  # ~90m
        mask=True,
    ).compute()
    ds = ds.rename(map="landcover")
    assert isinstance(ds, xr.Dataset), "Expected output to be an xarray Dataset."
    assert "landcover" in ds.data_vars, "Expected 'landcover' variable in the Dataset."


# Tests for `plot_esa_worldcover`
@pytest.mark.network
def test_plot_esa_worldcover_valid(search_aoi):
    """Test `plot_esa_worldcover` with valid WorldCover dataset."""
    gf_wc = coincident.search.search(
        dataset="worldcover",
        intersects=search_aoi,
        datetime=["2020"],
    )
    ds = to_dataset(
        gf_wc,
        bands=["map"],
        aoi=search_aoi,
        resolution=0.00081,  # ~90m
        mask=True,
    ).compute()
    ds = ds.rename(map="landcover")
    ax = plot_esa_worldcover(ds)
    assert ax is not None, "Expected a valid Matplotlib Axes object."
    assert len(ax.images) > 0, "Expected at least one image in the plot."
    ax.set_title("ESA WorldCover")
