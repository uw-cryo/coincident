from __future__ import annotations

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr
from matplotlib.collections import QuadMesh

import coincident.io.xarray
import coincident.plot
import coincident.search

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


def test_create_hillshade_tiny(dem_tiny):
    """Test 'create_hillshade' with a tiny DEM (10x10 cop30)"""

    data_elev = dem_tiny.elevation.squeeze()
    # Test cases with different azimuth and altitude combinations
    test_params = [
        (45, 45),  # default values
        (90, 30),
    ]
    for azi, alt in test_params:
        hillshade = coincident.plot.create_hillshade(data_elev, azi=azi, alt=alt)
        assert hillshade.shape == (10, 10), f"Shape incorrect for azi={azi}, alt={alt}"
        assert not np.isnan(
            hillshade
        ).any(), f"NaN values found for azi={azi}, alt={alt}"
        assert hillshade.dtype == np.uint8, f"Wrong dtype for azi={azi}, alt={alt}"
        assert hillshade.min() >= 0, f"Values < 0 found for azi={azi}, alt={alt}"
        assert hillshade.max() <= 255, f"Values > 255 found for azi={azi}, alt={alt}"


def test_clear_labels():
    """Test that clear_labels removes axis labels correctly"""

    fig, ax = plt.subplots()
    coincident.plot.clear_labels(ax)
    assert isinstance(
        ax.xaxis.get_major_formatter(), plt.NullFormatter
    ), "x-axis abels not cleared"
    assert isinstance(
        ax.yaxis.get_major_formatter(), plt.NullFormatter
    ), "y-axis abels not cleared"


def test_plot_dem_no_hillshade(dem_tiny):
    """Test plot_dem with tiny DEM input without hillshade"""

    fig, ax = plt.subplots()
    ax = coincident.plot.plot_dem(dem_tiny.elevation.squeeze(), ax, title="Test DEM")
    # Check basic plot properties
    assert ax.get_title() == "Test DEM", "Plot title does not match expected value"
    assert ax.images[0].get_alpha() == 0.5, "DEM layer alpha should be 0.5"
    assert isinstance(
        ax.xaxis.get_major_formatter(), plt.NullFormatter
    ), "x-axis abels not cleared"
    assert isinstance(
        ax.yaxis.get_major_formatter(), plt.NullFormatter
    ), "y-axis abels not cleared"
    assert isinstance(ax, plt.Axes), "Return value should be a matplotlib Axes object"


def test_plot_dem_with_hillshade(dem_tiny):
    """Test plot_dem with tiny DEM input with hillshade"""

    data_elev = dem_tiny.elevation.squeeze()
    dem_tiny["hillshade"] = (("y", "x"), coincident.plot.create_hillshade(data_elev))
    fig, ax = plt.subplots()
    coincident.plot.plot_dem(data_elev, ax, da_hillshade=dem_tiny.hillshade.squeeze())
    # Check that both hillshade and DEM layers are present
    assert isinstance(ax, plt.Axes), "Return value should be a matplotlib Axes object"
    assert len(ax.images) == 2, "Plot should contain exactly 2 image layers"
    assert ax.images[0].get_alpha() == 1.0, "Hillshade layer alpha should be 1.0"
    assert ax.images[1].get_alpha() == 0.5, "DEM layer alpha should be 0.5"
    assert ax.get_aspect() == 1.0, "Plot aspect ratio should be equal"


def test_plot_altimeter_points_no_hillshade(points_tiny):
    """Test plot_altimeter_points with point data without hillshade"""
    fig, ax = plt.subplots()
    ax = coincident.plot.plot_altimeter_points(
        points_tiny, ax, column="h_li", title="Test Points"
    )

    assert isinstance(ax, plt.Axes), "Return value should be a matplotlib Axes object"
    assert ax.get_title() == "Test Points", "Plot title does not match expected value"
    assert ax.collections[0].get_alpha() == 0.5, "Points layer alpha should be 0.5"
    assert isinstance(
        ax.xaxis.get_major_formatter(), plt.NullFormatter
    ), "x-axis labels not cleared"
    assert isinstance(
        ax.yaxis.get_major_formatter(), plt.NullFormatter
    ), "y-axis labels not cleared"


def test_plot_altimeter_points_with_hillshade(dem_tiny, points_tiny):
    """Test plot_altimeter_points with point data and hillshade background"""
    dem_tiny["hillshade"] = (
        ("y", "x"),
        coincident.plot.create_hillshade(dem_tiny.elevation.squeeze()),
    )

    fig, ax = plt.subplots()
    ax = coincident.plot.plot_altimeter_points(
        points_tiny, ax, column="h_li", da_hillshade=dem_tiny.hillshade
    )
    assert isinstance(ax, plt.Axes), "Return value should be a matplotlib Axes object"
    assert len(ax.images) == 1, "Plot should contain exactly 1 hillshade layer"
    assert len(ax.collections) == 1, "Plot should contain exactly 1 point collection"
    assert ax.images[0].get_alpha() == 1.0, "Hillshade layer alpha should be 1.0"
    assert ax.collections[0].get_alpha() == 0.5, "Points layer alpha should be 0.5"
    assert ax.get_aspect() == 1.0, "Plot aspect ratio should be equal"


def test_get_elev_diff_gf(dem_tiny, points_tiny):
    """Test get_elev_diff with point data source"""
    # TODO: add test to make sure differencing returns expected values
    gf_diff = coincident.plot.get_elev_diff(points_tiny, dem_tiny, source_col="h_li")
    assert type(gf_diff) is gpd.GeoDataFrame, "Returned object should be a GeoDataFrame"
    assert len(gf_diff) == 10, "Returned GeoDataFrame should have 10 rows"
    assert (
        "elev_diff" in gf_diff.columns
    ), "GeoDataFrame should contain elev_diff column"
    assert all(
        geom.geom_type == "Point" for geom in gf_diff.geometry
    ), "All geometries should be Points"
    assert gf_diff["elev_diff"].dtype == float, "elev_diff column should be float type"


def test_get_elev_diff_ds(dem_tiny):
    """Test get_elev_diff with raster data source"""
    # TODO: add test to make sure differencing returns expected values
    dem_tiny_2 = dem_tiny.copy()
    random_elevations = np.random.uniform(2935, 2965, size=(1, 10, 10))  # noqa: NPY002
    dem_tiny_2["elevation"] = (("band", "y", "x"), random_elevations)
    ds_diff = coincident.plot.get_elev_diff(dem_tiny_2, dem_tiny)
    assert type(ds_diff) is xr.Dataset, "Returned object should be an xr dataset"
    assert "elev_diff" in ds_diff.data_vars, "Dataset should contain elev_diff variable"
    assert ds_diff.elev_diff.shape == (1, 10, 10), "elev_diff should have 10x10 grid"
    assert (
        ds_diff["elev_diff"].dtype == float
    ), "elev_diff variable should be float type"


def test_plot_diff_hist_raster(dem_tiny):
    """Test get_elev_diff with both raster data source"""
    dem_tiny_2 = dem_tiny.copy()
    random_elevations = np.random.uniform(2935, 2965, size=(1, 10, 10))  # noqa: NPY002
    dem_tiny_2["elevation"] = (("band", "y", "x"), random_elevations)
    ds_diff = coincident.plot.get_elev_diff(dem_tiny_2, dem_tiny)

    f, ax = plt.subplots()
    ax = coincident.plot.plot_diff_hist(ds_diff.elev_diff, ax=ax)
    assert isinstance(ax, plt.Axes), "Return value should be a matplotlib Axes object"
    assert len(ax.get_children()) > 0, "Figure should exist"
    assert ax.get_legend() is not None, "Legend should exist"


def test_plot_diff_hist_point(points_tiny, dem_tiny):
    """Test get_elev_diff with both point data source"""

    gf_diff = coincident.plot.get_elev_diff(points_tiny, dem_tiny, source_col="h_li")
    f, ax = plt.subplots()
    ax = coincident.plot.plot_diff_hist(gf_diff.elev_diff, ax=ax)
    assert isinstance(ax, plt.Axes), "Return value should be a matplotlib Axes object"
    assert len(ax.get_children()) > 0, "Figure should exist"
    assert ax.get_legend() is not None, "Legend should exist"


# network to grab ESA worldcover
@network
def test_compare_dems(dem_tiny, points_tiny):
    """Test compare_dems with permutations of different numbers of DEMs and point gfs"""
    # Takes ~8secs to run
    # Create multiple DEMs
    dem_tiny["hillshade"] = (
        ("y", "x"),
        coincident.plot.create_hillshade(dem_tiny.elevation.squeeze()),
    )
    dem_tiny_2 = dem_tiny.copy()
    dem_tiny_3 = dem_tiny.copy()

    # Test with 3 DEMs and no GeoDataFrames
    axd_3_dems = coincident.plot.compare_dems(
        [dem_tiny, dem_tiny_2, dem_tiny_3], show=False
    )
    assert len(axd_3_dems) == 9, "Expected 9 axes"
    expected_keys = [
        "dem_0",
        "dem_1",
        "dem_2",
        "worldcover",
        "diff_1",
        "diff_2",
        "wc_legend",
        "hist_1",
        "hist_2",
    ]
    assert all(key in axd_3_dems for key in expected_keys), "Missing expected axes"
    for i in range(3):
        assert len(axd_3_dems[f"dem_{i}"].images) > 0, f"dem_{i} should have image data"
    for i in range(1, 3):
        assert (
            len(axd_3_dems[f"diff_{i}"].images) > 0
        ), f"diff_{i} should have image data"
    for i in range(1, 3):
        assert (
            len(axd_3_dems[f"hist_{i}"].patches) > 0
        ), f"hist_{i} should have histogram bars"

    # Test with 3 DEMs and 2 GeoDataFrames
    axd_full = coincident.plot.compare_dems(
        [dem_tiny, dem_tiny_2, dem_tiny_3],
        {"is2": (points_tiny, "h_li"), "gedi": (points_tiny, "elevation_hr")},
        show=False,
    )
    assert len(axd_full) == 15, "Expected 15 axes"
    expected_keys = [
        "dem_0",
        "dem_1",
        "dem_2",
        "is2",
        "gedi",
        "worldcover",
        "diff_1",
        "diff_2",
        "diff_is2",
        "diff_gedi",
        "wc_legend",
        "hist_1",
        "hist_2",
        "hist_is2",
        "hist_gedi",
    ]
    assert all(key in axd_full for key in expected_keys), "Missing expected axes"
    for i in range(3):
        assert len(axd_full[f"dem_{i}"].images) > 0, f"dem_{i} should have image data"
    for i in range(1, 3):
        assert len(axd_full[f"diff_{i}"].images) > 0, f"diff_{i} should have image data"
    for i in range(1, 3):
        assert (
            len(axd_full[f"hist_{i}"].patches) > 0
        ), f"hist_{i} should have histogram bars"
    for key in ["is2", "gedi"]:
        assert (
            len(axd_full[f"diff_{key}"].collections) > 0
        ), f"diff_{key} should have point data"
    for i in range(1, 3):
        assert (
            len(axd_full[f"hist_{i}"].patches) > 0
        ), f"hist_{i} should have histogram bars"
    for key in ["is2", "gedi"]:
        assert (
            len(axd_full[f"hist_{key}"].patches) > 0
        ), f"hist_{key} should have histogram bars"
