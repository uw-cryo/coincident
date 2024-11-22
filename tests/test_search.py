from __future__ import annotations

import typing

import geopandas as gpd
import pytest
from geopandas.testing import assert_geodataframe_equal

import coincident

# Decorate tests requiring internet (slow & flaky)
network = pytest.mark.network

try:
    import maxar_platform.discovery  # noqa: F401

    not_authenticated = False
except:  # noqa: E722
    not_authenticated = True
maxar_authenticated = pytest.mark.skipif(
    not_authenticated, reason="tests for linux only"
)


@pytest.fixture
def aoi():
    # 11 vertices, 1,361km^2
    aoi_url = "https://raw.githubusercontent.com/SlideRuleEarth/sliderule-python/main/data/grandmesa.geojson"
    return gpd.read_file(aoi_url)


@pytest.fixture
def large_aoi():
    # 260 vertices, large area 269,590 km^2
    aoi_url = "https://raw.githubusercontent.com/unitedstates/districts/refs/heads/gh-pages/states/CO/shape.geojson"
    return gpd.read_file(aoi_url)


@typing.no_type_check
def test_no_dataset_specified():
    with pytest.raises(
        TypeError, match="missing 1 required positional argument: 'dataset'"
    ):
        coincident.search.search(intersects="-120, 40, -121, 41")


def test_unknown_dataset_specified():
    with pytest.raises(ValueError, match="is not a supported dataset"):
        coincident.search.search(dataset="typo", intersects="-120, 40, -121, 41")


def test_polygon_invalid_type():
    with pytest.raises(
        ValueError, match="intersects value must be a GeoDataFrame or GeoSeries"
    ):
        coincident.search.search(dataset="3dep", intersects="-120, 40, -121, 41")


def test_to_geopandas_empty_search_result():
    with pytest.raises(ValueError, match="ItemCollection is empty"):
        coincident.search.stac.to_geopandas([])


def test_unconstrained_search_warns():
    with pytest.warns(match="Neither `bbox` nor `intersects` provided"):
        coincident.search.search(dataset="tdx")


@network
def test_cascading_search(aoi):
    aoi["datetime"] = gpd.pd.to_datetime("2019-06-12", utc=True)
    pad = 30
    secondary_datasets = [("icesat-2", pad), ("gedi", pad)]
    results = coincident.search.cascading_search(aoi, secondary_datasets)

    expected_min = aoi.datetime.iloc[0] - gpd.pd.Timedelta(days=pad)
    expected_max = aoi.datetime.iloc[0] + gpd.pd.Timedelta(days=pad)
    actual_min = results[0].datetime.min()
    actual_max = results[0].datetime.max()

    assert isinstance(results, list)
    assert len(results) == 2
    assert len(results[0]) == 5
    assert len(results[1]) == 4
    assert actual_min >= expected_min
    assert actual_max <= expected_max


# TODO: add more assertions / tests for this section
@network
@maxar_authenticated
@pytest.mark.filterwarnings("ignore:Server does not conform")
def test_maxar_search(aoi):
    gf = coincident.search.search(
        dataset="maxar",
        intersects=aoi,
        datetime="2023",
        filter="eo:cloud_cover < 20",
    )
    assert len(gf) == 10
    assert isinstance(gf.datetime.iloc[0], gpd.pd.Timestamp)
    # NOTE: add more comprehensive check of column names & schema?
    assert "id" in gf.columns
    assert "assets" in gf.columns


@network
@maxar_authenticated
def test_maxar_large_aoi(large_aoi):
    gf = coincident.search.search(
        dataset="maxar",
        intersects=large_aoi,
        datetime="2023",
        # NOTE: this limits total *mono* results, which may or may not have stereo
        max_items=10,
    )
    assert len(gf) <= 10


# NASA
# =======
@network
def test_icesat2_search(aoi):
    gf = coincident.search.search(
        dataset="icesat-2",
        intersects=aoi,
        datetime="2023",
    )
    assert len(gf) == 25
    assert isinstance(gf.start_datetime.iloc[0], gpd.pd.Timestamp)


@network
def test_gedi_search(aoi):
    gf = coincident.search.search(
        dataset="gedi",
        intersects=aoi,
        datetime="2022",
    )
    assert len(gf) == 33


@network
def test_tdx_search(aoi):
    gf = coincident.search.search(
        dataset="tdx", intersects=aoi, datetime=["2009", "2020"]
    )
    assert len(gf) == 48
    assert gf["sar:product_type"].unique() == "SSC"


# MS PLANETARY COMPUTER
# =======
@network
def test_cop30_search(aoi):
    gf = coincident.search.search(dataset="cop30", intersects=aoi)
    assert len(gf) == 4


@network
def test_worldcover_search(aoi):
    gf = coincident.search.search(dataset="worldcover", intersects=aoi, datetime="2020")
    assert len(gf) == 4


@network
def test_round_trip_parquet(aoi):
    outpath = "/tmp/search_results.parquet"
    A = coincident.search.search(dataset="cop30", intersects=aoi)
    A.to_parquet(outpath)
    B = gpd.read_parquet(outpath)
    assert_geodataframe_equal(A, B)


# USGS
# =======
@network
def test_wesm_search(aoi):
    gf = coincident.search.search(
        dataset="3dep",
        intersects=aoi,
    )
    assert len(gf) == 5


# NOTE ~10s on wifi
@network
def test_get_swath_polygon():
    gf = coincident.search.wesm.get_swath_polygons("CO_CameronPkFire_1_2021")
    assert isinstance(gf, gpd.GeoDataFrame)
    assert len(gf) == 51
    assert "start_datetime" in gf.columns
    assert isinstance(gf.datetime.iloc[0], gpd.pd.Timestamp)


@network
def test_swath_polygon_not_found():
    with pytest.raises(ValueError, match="No swath polygons found for workunit="):
        coincident.search.wesm.get_swath_polygons("AL_SWCentral_1_B22")
