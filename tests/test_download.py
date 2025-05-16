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
    asyncio.run(coincident.io.download.download_item(item))
    assert Path("/tmp/102001008EC5AC00.browse.tif").exists()
    assert Path("/tmp/102001008EC5AC00.json").exists()


# @network
# def test_download_planetary_computer():
#     print('todo: small cop30 asset?')


# NOTE: I believe all of the "test_xxxx_dem" functions to be redundant since they rely on the
#       same logic as the "load_xxxx_dem" functions from the xarray module
# # ~30 secs on my laptop regardless of number of tiles downloaded (where I only tested up to 6 tiles, current func just has one tile)
# @network
# def test_3dep_dem(usgs_neon_dem_bbox):
#     coincident.io.download.download_usgs_dem(
#         aoi=usgs_neon_dem_bbox,
#         project="CO_CentralEasternPlains_2020_D20",
#         output_dir="/tmp",
#         save_parquet=True,
#     )
#     assert Path("/tmp/USGS_1M_13_x71y440_CO_CentralEasternPlains_2020_D20.tif").exists()
#     assert Path("/tmp/stac_tnm_CO_CentralEasternPlains_2020_D20.parquet").exists()
#     with rasterio.open(
#         "/tmp/USGS_1M_13_x71y440_CO_CentralEasternPlains_2020_D20.tif"
#     ) as src:
#         assert src.height > 0, "TIF has no height"
#         assert src.width > 0, "TIF has no width"
#         assert src.count > 0, "TIF has no bands"
#     gf_stac = gpd.read_parquet("/tmp/stac_tnm_CO_CentralEasternPlains_2020_D20.parquet")
#     assert gf_stac.shape[0] > 0, "Parquet file has no rows"
#     assert "geometry" in gf_stac.columns, "Parquet file missing geometry column"


@network
def test_3dep_tile_gdf(usgs_neon_dem_bbox):
    gf_usgs_lpc_tiles = coincident.io.download.fetch_usgs_lpc_tiles(
        aoi=usgs_neon_dem_bbox, project="CO_CentralEasternPlains_2020_D20"
    )
    assert gf_usgs_lpc_tiles.shape == (1, 3), (
        f"Unexpected shape: {gf_usgs_lpc_tiles.shape}, expected (1, 3)"
    )
    expected_columns = ["name", "url", "geometry"]
    assert all(col in gf_usgs_lpc_tiles.columns for col in expected_columns), (
        f"Missing columns. Found: {gf_usgs_lpc_tiles.columns}"
    )
    for col in expected_columns:
        assert gf_usgs_lpc_tiles[col].notna().all(), (
            f"Column {col} contains null values"
        )
    assert not gf_usgs_lpc_tiles.geometry.is_empty.any(), "Contains empty geometries"


@network
def test_3dep_ept(usgs_neon_dem_bbox):
    json_pdal_pipeline = coincident.io.download.build_usgs_ept_pipeline(
        usgs_neon_dem_bbox, workunit="CO_CentralEasternPlains_1_2020"
    )
    assert len(json_pdal_pipeline["pipeline"]) == 3
    assert all("type" in stage for stage in json_pdal_pipeline["pipeline"])


# @network
# def test_neon_dem(usgs_neon_dem_bbox):
#     coincident.io.download.download_neon_dem(
#         aoi=usgs_neon_dem_bbox,
#         datetime_str="2020-06-30",
#         site_id="ARIK",
#         product="dsm",
#         output_dir="/tmp",
#     )
#     assert Path("/tmp/NEON_D10_ARIK_DP3_713000_4394000_DSM.tif").exists()
#     with rasterio.open("/tmp/NEON_D10_ARIK_DP3_713000_4394000_DSM.tif") as src:
#         assert src.height > 0, "TIF has no height"
#         assert src.width > 0, "TIF has no width"
#         assert src.count > 0, "TIF has no bands"


@network
def test_neon_tile_gdf(usgs_neon_dem_bbox):
    gf_neon_lpc_tiles = coincident.io.download.fetch_neon_lpc_tiles(
        aoi=usgs_neon_dem_bbox, datetime_str="2020-06-30", site_id="ARIK"
    )
    assert gf_neon_lpc_tiles.shape == (1, 3), (
        f"Unexpected shape: {gf_neon_lpc_tiles.shape}, expected (1, 3)"
    )
    expected_columns = ["name", "url", "geometry"]
    assert all(col in gf_neon_lpc_tiles.columns for col in expected_columns), (
        f"Missing columns. Found: {gf_neon_lpc_tiles.columns}"
    )
    for col in expected_columns:
        assert gf_neon_lpc_tiles[col].notna().all(), (
            f"Column {col} contains null values"
        )
    assert not gf_neon_lpc_tiles.geometry.is_empty.any(), "Contains empty geometries"


# # 40 seconds on my laptop (NCALM flights are not tiled)
# # TODO: find a NCALM example that's very small
# @network
# def test_ncalm_dem(ncalm_dem_bbox):
#     coincident.io.download.download_ncalm_dem(aoi=ncalm_dem_bbox,
#                                                 dataset_id='WA18_Wall',
#                                                 product='dtm',
#                                                 output_dir='/tmp')
#     assert Path("/tmp/clipped_WALL_GEG_1M_dtm.tif").exists()
#     with rasterio.open("/tmp/clipped_WALL_GEG_1M_dtm.tif") as src:
#         assert src.height > 0, "TIF has no height"
#         assert src.width > 0, "TIF has no width"
#         assert src.count > 0, "TIF has no bands"


@network
def test_ncalm_tile_gdf(ncalm_dem_bbox):
    gf_ncalm_lpc_tiles = coincident.io.download.fetch_ncalm_lpc_tiles(
        aoi=ncalm_dem_bbox, dataset_name="WA18_Wall"
    )
    assert gf_ncalm_lpc_tiles.shape == (1, 3), (
        f"Unexpected shape: {gf_ncalm_lpc_tiles.shape}, expected (1, 3)"
    )
    expected_columns = ["name", "url", "geometry"]
    assert all(col in gf_ncalm_lpc_tiles.columns for col in expected_columns), (
        f"Missing columns. Found: {gf_ncalm_lpc_tiles.columns}"
    )
    for col in expected_columns:
        assert gf_ncalm_lpc_tiles[col].notna().all(), (
            f"Column {col} contains null values"
        )
    assert not gf_ncalm_lpc_tiles.geometry.is_empty.any(), "Contains empty geometries"


# 30 seconds on my laptop (time depends on total NOAA flight area)
# TODO: find a NOAA example that's very small
# @network
# def test_noaa_dem(noaa_dem_bbox):
#    coincident.io.download.download_noaa_dem(aoi=noaa_dem_bbox,
#                                            dataset_id=10150,
#                                            output_dir="/tmp")
#    assert Path("/tmp/clipped_2022_220000e_3305000n_dem.tif").exists()
#    with rasterio.open("/tmp/clipped_2022_225000e_3305000n_dem.tif") as src:
#        assert src.height > 0, "TIF has no height"
#        assert src.width > 0, "TIF has no width"
#        assert src.count > 0, "TIF has no bands"
#    assert Path("/tmp/clipped_2022_225000e_3305000n_dem.tif").exists()
#    with rasterio.open("/tmp/clipped_2022_225000e_3305000n_dem.tif") as src:
#        assert src.height > 0, "TIF has no height"
#        assert src.width > 0, "TIF has no width"
#        assert src.count > 0, "TIF has no bands"


@network
def test_noaa_tile_gdf(noaa_dem_bbox):
    gf_noaa_lpc_tiles = coincident.io.download.fetch_noaa_lpc_tiles(
        aoi=noaa_dem_bbox, dataset_id="10149"
    )
    assert gf_noaa_lpc_tiles.shape == (3, 3), (
        f"Unexpected shape: {gf_noaa_lpc_tiles.shape}, expected (1, 3)"
    )
    expected_columns = ["name", "url", "geometry"]
    assert all(col in gf_noaa_lpc_tiles.columns for col in expected_columns), (
        f"Missing columns. Found: {gf_noaa_lpc_tiles.columns}"
    )
    for col in expected_columns:
        assert gf_noaa_lpc_tiles[col].notna().all(), (
            f"Column {col} contains null values"
        )
    assert not gf_noaa_lpc_tiles.geometry.is_empty.any(), "Contains empty geometries"
