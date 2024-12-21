from __future__ import annotations

import geopandas as gpd
import pytest

import coincident

# import sliderule
# sliderule.init(verbose=True)

network = pytest.mark.network


def test_gdf_to_sliderule_polygon():
    gf = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(
            [-106.46953038438801, -106.47209604357627],
            [38.781089777057986, 38.78292421778072],
            crs="EPSG:4326",
        )
    )
    poly = coincident.io.sliderule._gdf_to_sliderule_polygon(gf)
    assert isinstance(poly, list)
    assert len(poly) == 5


@network
@pytest.mark.skip(reason="https://github.com/uw-cryo/coincident/issues/35")
class TestSlideRule:
    def test_subset_gedi02a(self, tinyaoi):
        gf_gedi = coincident.search.search(
            dataset="gedi", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.subset_gedi02a(
            gf_gedi, aoi=tinyaoi, include_worldcover=True
        )
        assert isinstance(data, gpd.GeoDataFrame)
        assert "elevation_lm" in data.columns
        assert "elevation_hr" in data.columns
        assert data.iloc[0].name.toordinal() == 738192
        assert data.iloc[0]["worldcover.value"] == "Tree cover"
        assert data.shape == (262, 11)

    def test_subset_atl06(self, tinyaoi):
        gf_is2 = coincident.search.search(
            dataset="icesat-2", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.subset_atl06(gf_is2, aoi=tinyaoi)
        assert isinstance(data, gpd.GeoDataFrame)
        assert "h_li" in data.columns
        assert 1034 in data.rgt.unique()
        assert data.iloc[0].name.toordinal() == 738304
        assert data.shape == (413, 21)

    def test_process_atl06sr(self, tinyaoi):
        gf_is2 = coincident.search.search(
            dataset="icesat-2", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.process_atl06sr(
            gf_is2, aoi=tinyaoi, include_3dep=True
        )
        assert isinstance(data, gpd.GeoDataFrame)
        assert "n_fit_photons" in data.columns
        assert "3dep.value" in data.columns
        assert data.shape == (367, 20)
        assert data.iloc[0].name.toordinal() == 738304
        assert int(data.iloc[0].h_mean) == 3178

    @pytest.mark.filterwarnings(
        "ignore:Geometry is in a geographic CRS:UserWarning:geopandas"
    )
    def test_sample_3dep(self):
        gf = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy(
                [-106.46953038438801, -106.47209604357627],
                [38.781089777057986, 38.78292421778072],
                crs="EPSG:4326",
            )
        )
        samples = coincident.io.sliderule.sample_3dep(gf)
        assert {"file", "geometry", "time", "value"} == set(samples.columns)
        assert samples.time.dtype == "<M8[ns]"
        assert str(samples.time.dt.date.iloc[0]) == "2019-09-04"
        # NOTE: is this on edge of 2 tiles?! sometimes comes back as {3110, 3178}...
        # assert set(samples.value.astype('i2')) == {3111, 3178}
