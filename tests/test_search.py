from __future__ import annotations

import typing
from importlib.metadata import version

import geopandas as gpd
import pytest
from geopandas.testing import assert_geodataframe_equal
from shapely.geometry import Polygon

import coincident

# Decorate tests requiring internet (slow & flaky)
network = pytest.mark.network

try:
    import maxar_platform.discovery  # noqa: F401

    not_authenticated = False
except:  # noqa: E722
    not_authenticated = True
maxar_authenticated = pytest.mark.skipif(
    not_authenticated, reason="Not authenticated with Maxar API"
)

expected_stac_version = "1.1.0" if version("pystac") >= "1.12" else "1.0.0"


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


@network
def test_3d_poly_search():
    # define a 3D Polygon for search
    bbox = [-105, 40, -104, 41]
    gf_polygon_z = gpd.GeoDataFrame(
        geometry=[
            Polygon(
                [
                    [bbox[0], bbox[1], 0],
                    [bbox[2], bbox[1], 0],
                    [bbox[2], bbox[3], 0],
                    [bbox[0], bbox[3], 0],
                ]
            )
        ]
    )
    gf_is2_stac = coincident.search.search(
        dataset="icesat-2",
        intersects=gf_polygon_z,
        datetime=["2020-01-01", "2020-02-01"],
    )
    assert isinstance(gf_is2_stac, gpd.GeoDataFrame)
    assert len(gf_is2_stac) == 4


# MAXAR
# =======
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
    assert gf.shape == (10, 69)
    assert gf.iloc[0].stac_version == expected_stac_version
    assert {"browse", "cloud-cover", "sample-point-set"}.issubset(
        gf.iloc[0].assets.keys()
    )
    assert isinstance(gf.datetime.iloc[0], gpd.pd.Timestamp)
    # NOTE: add more here depending on what is needed in lidar_tools
    assert {"id", "dayofyear", "stereo_pair_identifiers", "platform"}.issubset(
        set(gf.columns)
    )


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


@network
@maxar_authenticated
def test_maxar_specific_ids():
    gf = coincident.search.search(
        dataset="maxar", ids=["102001008EC5AC00", "102001008BE9BB00"]
    )
    assert set(gf.id) == {"102001008BE9BB00", "102001008EC5AC00"}


# NASA
# =======
expected_nasa_columns = {
    "assets",
    "bbox",
    "collection",
    "geometry",
    "id",
    "links",
    "stac_extensions",
    "stac_version",
    "type",
    "datetime",
    "end_datetime",
    "start_datetime",
    "dayofyear",
}


@network
def test_icesat2_search(aoi):
    gf = coincident.search.search(
        dataset="icesat-2",
        intersects=aoi,
        datetime="2023",
    )
    actual_columns = set(gf.columns)
    data_assets = list(
        filter(lambda x: x["roles"] == "data", gf.iloc[0].assets.values())
    )
    assert gf.iloc[0].stac_version == expected_stac_version
    assert gf.shape == (25, 13)
    assert actual_columns == expected_nasa_columns
    assert "roles" in gf.iloc[0].assets["browse"]
    assert len(data_assets) == 1
    assert gf.iloc[0].collection.startswith("ATL03")
    assert isinstance(gf.start_datetime.iloc[0], gpd.pd.Timestamp)


@network
def test_gedi_search(aoi):
    gf = coincident.search.search(
        dataset="gedi",
        intersects=aoi,
        datetime="2022",
    )
    actual_columns = set(gf.columns)
    data_assets = list(
        filter(lambda x: x["roles"] == "data", gf.iloc[0].assets.values())
    )
    assert gf.iloc[0].stac_version == expected_stac_version
    assert gf.shape == (33, 13)
    assert actual_columns == expected_nasa_columns
    assert "roles" in gf.iloc[0].assets["browse"]
    assert len(data_assets) == 1
    assert gf.iloc[0].collection.startswith("GEDI02_A")
    assert isinstance(gf.start_datetime.iloc[0], gpd.pd.Timestamp)


@network
def test_gliht_search():
    from shapely.geometry import box

    aoi = gpd.GeoDataFrame(
        geometry=[box(*[-71.907258, 41.097413, -71.088571, 42.018798])], crs="EPSG:4326"
    )
    gf = coincident.search.search(dataset="gliht", intersects=aoi, datetime=["2012"])
    actual_columns = set(gf.columns)
    data_assets = list(
        filter(lambda x: x["roles"] == "data", gf.iloc[0].assets.values())
    )
    assert gf.iloc[0].stac_version == expected_stac_version
    assert gf.shape == (10, 13)
    assert actual_columns == expected_nasa_columns
    assert "roles" in gf.iloc[0].assets["browse"]
    assert len(data_assets) == 3
    assert gf.iloc[0].collection.startswith("GLDSMT_001")
    assert isinstance(gf.start_datetime.iloc[0], gpd.pd.Timestamp)


# NASA/CSDA
# =======
@network
def test_tdx_search(aoi):
    gf = coincident.search.search(
        dataset="tdx", intersects=aoi, datetime=["2009", "2020"]
    )
    expected_columns = {
        "assets",
        "bbox",
        "collection",
        "geometry",
        "id",
        "links",
        "stac_extensions",
        "stac_version",
        "type",
        "constellation",
        "datetime",
        "end_datetime",
        "platform",
        "sar:center_frequency",
        "sar:frequency_band",
        "sar:instrument_mode",
        "sar:looks_azimuth",
        "sar:looks_range",
        "sar:polarizations",
        "sar:product_type",
        "sar:resolution_azimuth",
        "sar:resolution_range",
        "start_datetime",
        "dayofyear",
    }
    actual_columns = set(gf.columns)
    # NOTE: I think technically assets should have 'roles' not 'role' key...
    data_assets = list(
        filter(lambda x: x["role"] == "data", gf.iloc[0].assets.values())
    )
    assert gf.iloc[0].stac_version == expected_stac_version
    assert gf.shape == (48, 24)
    assert actual_columns == expected_columns
    assert gf["sar:product_type"].unique() == "SSC"
    assert len(data_assets) >= 1


# MS PLANETARY COMPUTER
# =======
@network
def test_cop30_search(aoi):
    gf = coincident.search.search(dataset="cop30", intersects=aoi)
    expected_columns = {
        "assets",
        "bbox",
        "collection",
        "geometry",
        "id",
        "links",
        "stac_extensions",
        "stac_version",
        "type",
        "datetime",
        "gsd",
        "platform",
        "proj:code",
        "proj:shape",
        "proj:transform",
        "dayofyear",
    }
    actual_columns = set(gf.columns)
    assert gf.iloc[0].stac_version == expected_stac_version
    assert gf.shape == (4, 16)
    assert actual_columns == expected_columns
    assert "data" in gf.iloc[0].assets
    assert gf.iloc[0]["proj:code"] == "EPSG:4326"


@network
def test_worldcover_search(aoi):
    gf = coincident.search.search(dataset="worldcover", intersects=aoi, datetime="2020")
    # NOTE: could test these separately...
    stac_columns = {
        "assets",
        "bbox",
        "collection",
        "geometry",
        "id",
        "links",
        "stac_extensions",
        "stac_version",
        "type",
        "datetime",
    }
    additional_columns = {
        "created",
        "description",
        "end_datetime",
        "esa_worldcover:product_tile",
        "esa_worldcover:product_version",
        "grid:code",
        "instruments",
        "mission",
        "platform",
        "proj:code",
        "start_datetime",
    }
    coincident_columns = {"dayofyear"}
    expected_columns = stac_columns | additional_columns | coincident_columns
    actual_columns = set(gf.columns)
    assert gf.iloc[0].stac_version == expected_stac_version
    assert gf.shape == (4, 22)
    assert actual_columns == expected_columns
    assert "map" in gf.iloc[0].assets


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
    expected_columns = {
        "workunit",
        "project",
        "start_datetime",
        "end_datetime",
        "dayofyear",
        "duration",
        "p_method",
    }
    actual_columns = set(gf.columns)
    assert gf.shape == (5, 37)
    assert expected_columns.issubset(actual_columns)


# NOTE can take ~10s on wifi... maybe only run this on demand...
@network
def test_get_swath_polygon():
    gf = coincident.search.wesm.get_swath_polygons("CO_CameronPkFire_1_2021")
    expected_columns = {
        "start_datetime",
        "end_datetime",
        "collection",
        "datetime",
        "dayofyear",
        "duration",
    }
    actual_columns = set(gf.columns)
    assert isinstance(gf, gpd.GeoDataFrame)
    assert gf.shape == (51, 16)
    assert expected_columns.issubset(actual_columns)
    assert isinstance(gf.datetime.iloc[0], gpd.pd.Timestamp)


@network
def test_swath_polygon_not_found():
    with pytest.raises(ValueError, match="No swath polygons found for workunit="):
        coincident.search.wesm.get_swath_polygons("AL_SWCentral_1_B22")


# opentopo (NCALM and NOAA)
# =======
@network
def test_noaa_search(bathy_aoi):
    gf = coincident.search.search(
        dataset="noaa", intersects=bathy_aoi, datetime=["2019-01-01", "2023-12-31"]
    )

    expected_columns = {
        "id",
        "name",
        "title",
        "start_datetime",
        "end_datetime",
        "geometry",
    }
    assert gf.shape == (2, 6)
    assert set(gf.columns) == expected_columns
    assert all(isinstance(geom, Polygon) for geom in gf["geometry"])


@network
def test_ncalm_search(large_aoi):
    gf = coincident.search.search(
        dataset="ncalm", intersects=large_aoi, datetime=["2019-01-01", "2023-12-31"]
    )
    expected_columns = {
        "id",
        "name",
        "title",
        "start_datetime",
        "end_datetime",
        "geometry",
    }
    assert gf.shape == (6, 6)
    assert set(gf.columns) == expected_columns
    assert all(isinstance(geom, Polygon) for geom in gf["geometry"])


# NEON
# TODO: add a test for provisional NEON datasets
# =======
@network
@pytest.mark.filterwarnings("ignore:Geometry is in a geographic CRS:UserWarning")
def test_neon_search():
    intersects = gpd.read_file(
        "https://raw.githubusercontent.com/unitedstates/districts/refs/heads/gh-pages/states/MA/shape.geojson"
    )
    gf = coincident.search.search(
        dataset="neon", intersects=intersects, datetime=["2019"]
    )

    expected_columns = {
        "id",
        "title",
        "start_datetime",
        "end_datetime",
        "product_url",
        "geometry",
    }
    assert gf.shape == (2, 6)
    assert set(gf.columns) == expected_columns
    assert all(isinstance(geom, Polygon) for geom in gf["geometry"])
