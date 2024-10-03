from __future__ import annotations

import geopandas as gpd
import pytest

import coincident as m

# Decorate tests requiring internet
network = pytest.mark.network


@pytest.fixture
def aoi():
    aoi_url = "https://raw.githubusercontent.com/SlideRuleEarth/sliderule-python/main/data/grandmesa.geojson"
    return gpd.read_file(aoi_url)


def test_no_dataset_specified():
    with pytest.raises(
        TypeError, match="missing 1 required positional argument: 'dataset'"
    ):
        m.search.search(intersects="-120, 40, -121, 41")  # type: ignore[call-arg]


def test_unknown_dataset_specified():
    with pytest.raises(ValueError, match="is not a supported dataset"):
        m.search.search(dataset="typo", intersects="-120, 40, -121, 41")


def test_polygon_invalid_type():
    with pytest.raises(
        ValueError, match="intersects value must be a GeoDataFrame or GeoSeries"
    ):
        m.search.search(dataset="3dep", intersects="-120, 40, -121, 41")


# For now raise an error to not complicate things
def test_polygon_out_of_bounds():
    feature_coll = {
        "type": "FeatureCollection",
        "features": [
            {
                "id": "0",
                "properties": {"col1": "name1"},
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": (1.0, 2.0)},
            }
        ],
    }
    aoi = gpd.GeoDataFrame.from_features(feature_coll, crs="EPSG:4326")
    # with pytest.warns(UserWarning, match="Requested search polygon not within"):
    with pytest.raises(ValueError, match="Requested search polygon not within"):
        m.search.search(dataset="3dep", intersects=aoi)


# TODO: add more assertions / tests for this section
@network
@pytest.mark.filterwarnings("ignore:Server does not conform")
def test_maxar_search(aoi):
    gf = m.search.search(
        dataset="maxar",
        intersects=aoi,
        datetime="2023",
        filter="eo:cloud_cover < 20",
    )
    assert len(gf) == 10
    assert isinstance(gf.datetime.iloc[0], gpd.pd.Timestamp)
    assert "stac_id" in gf.columns
    assert "browse" in gf.columns


# NASA
# =======
@network
def test_icesat2_search(aoi):
    gf = m.search.search(
        dataset="icesat-2",
        intersects=aoi,
        datetime="2023",
    )
    # expected_cols = ['geometry', 'datetime', 'start_datetime', 'end_datetime', 'stac_id', 'data', 'browse'] # Also GEDI
    assert len(gf) == 25
    assert isinstance(gf.start_datetime.iloc[0], gpd.pd.Timestamp)


@network
def test_gedi_search(aoi):
    gf = m.search.search(
        dataset="gedi",
        intersects=aoi,
        datetime="2022",
    )
    assert len(gf) == 33


@network
def test_tdx_search(aoi):
    gf = m.search.search(dataset="tdx", intersects=aoi, datetime=["2009", "2020"])
    assert len(gf) == 48


# MS PLANETARY COMPUTER
# =======
@network
def test_cop30_search(aoi):
    gf = m.search.search(dataset="cop30", intersects=aoi)
    # expected_cols = ['geometry','gsd','datetime','platform','proj:epsg','proj:shape','proj:transform','stac_id','data']
    assert len(gf) == 4


@network
def test_worldcover_search(aoi):
    gf = m.search.search(dataset="worldcover", intersects=aoi, datetime="2020")
    assert len(gf) == 4


# USGS
# =======
@network
def test_wesm_search(aoi):
    gf = m.search.search(
        dataset="3dep",
        intersects=aoi,
    )
    assert len(gf) == 3
