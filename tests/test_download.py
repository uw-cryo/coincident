from __future__ import annotations

import asyncio
from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import box

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
    asyncio.run(coincident.io.download.download_item(item))
    assert Path("/tmp/102001008EC5AC00.browse.tif").exists()
    assert Path("/tmp/102001008EC5AC00.json").exists()


@network
def test_3dep_tile_gdf(usgs_neon_dem_bbox):
    gf_usgs_lpc_tiles = coincident.io.download._fetch_usgs_lpc_tiles(
        aoi=usgs_neon_dem_bbox, project="CO_CentralEasternPlains_2020_D20"
    )
    assert gf_usgs_lpc_tiles.shape == (1, 3)
    expected_columns = {"name", "url", "geometry"}
    assert set(gf_usgs_lpc_tiles.columns) == expected_columns
    for col in expected_columns:
        assert gf_usgs_lpc_tiles[col].notna().all(), (
            f"Column {col} contains null values"
        )
    assert not gf_usgs_lpc_tiles.geometry.is_empty.any()


@network
def test_3dep_ept(usgs_neon_dem_bbox):
    json_pdal_pipeline = coincident.io.download.build_usgs_ept_pipeline(
        usgs_neon_dem_bbox, workunit="CO_CentralEasternPlains_1_2020"
    )
    assert len(json_pdal_pipeline["pipeline"]) == 3
    assert all("type" in stage for stage in json_pdal_pipeline["pipeline"])


@network
def test_neon_tile_gdf(usgs_neon_dem_bbox):
    gf_neon_lpc_tiles = coincident.io.download._fetch_neon_lpc_tiles(
        aoi=usgs_neon_dem_bbox, datetime_str="2020-06-30", site_id="ARIK"
    )
    assert gf_neon_lpc_tiles.shape == (1, 3)
    expected_columns = {"name", "url", "geometry"}
    assert set(gf_neon_lpc_tiles.columns) == expected_columns
    for col in expected_columns:
        assert gf_neon_lpc_tiles[col].notna().all(), (
            f"Column {col} contains null values"
        )
    assert not gf_neon_lpc_tiles.geometry.is_empty.any()


@network
def test_ncalm_tile_gdf(ncalm_dem_bbox):
    gf_ncalm_lpc_tiles = coincident.io.download._fetch_ncalm_lpc_tiles(
        aoi=ncalm_dem_bbox, dataset_name="WA18_Wall"
    )
    assert gf_ncalm_lpc_tiles.shape == (1, 3)
    expected_columns = {"name", "url", "geometry"}
    assert set(gf_ncalm_lpc_tiles.columns) == expected_columns
    for col in expected_columns:
        assert gf_ncalm_lpc_tiles[col].notna().all(), (
            f"Column {col} contains null values"
        )
    assert not gf_ncalm_lpc_tiles.geometry.is_empty.any()


@network
def test_noaa_tile_gdf():
    noaa_dem_bbox = gpd.GeoDataFrame(
        geometry=[box(-83.8468, 29.8067, -83.8308, 29.8227)], crs="EPSG:4326"
    )
    gf_noaa_lpc_tiles = coincident.io.download._fetch_noaa_lpc_tiles(
        aoi=noaa_dem_bbox, dataset_id="10149"
    )
    assert gf_noaa_lpc_tiles.shape == (3, 3), (
        f"Unexpected shape: {gf_noaa_lpc_tiles.shape}, expected (1, 3)"
    )
    expected_columns = {"name", "url", "geometry"}
    assert set(gf_noaa_lpc_tiles.columns) == expected_columns
    for col in expected_columns:
        assert gf_noaa_lpc_tiles[col].notna().all(), (
            f"Column {col} contains null values"
        )
    assert not gf_noaa_lpc_tiles.geometry.is_empty.any()
