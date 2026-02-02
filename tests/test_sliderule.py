from __future__ import annotations

import os
import pprint

import geopandas as gpd
import pytest
import sliderule

import coincident

network = pytest.mark.network


@network
# Always dump sliderule version info to stdout for debugging
def test_sliderule_versions(capsys):
    with capsys.disabled():
        if os.environ.get("SLIDERULE_GITHUB_TOKEN"):
            message = "SLIDERULE_GITHUB_TOKEN is set, so using Private Cluster instead of public cluster"
            print("\n-------\n" + message + "\n-------\n")
        version_dict = sliderule.get_version()
        print("\nSliderule Version Info:\n")
        pprint.pprint(version_dict)
    assert True


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


@pytest.fixture(scope="class")
def expected_gedi_columns():
    return {
        "beam",
        "elevation_hr",
        "elevation_lm",
        "flags",
        "geometry",
        "orbit",
        "sensitivity",
        "shot_number",
        "solar_elevation",
        "srcid",
        "track",
    }


@pytest.fixture(scope="class")
def expected_is2_columns():
    return {
        "atl06_quality_summary",
        "bsnow_conf",
        "bsnow_h",
        "cycle",
        "dh_fit_dx",
        "geometry",
        "gt",
        "h_li",
        "h_li_sigma",
        "h_robust_sprd",
        "n_fit_photons",
        "r_eff",
        "region",
        "rgt",
        "seg_azimuth",
        "segment_id",
        "sigma_geo_h",
        "spot",
        "srcid",
        "tide_ocean",
        "w_surface_window_final",
        "x_atc",
        "y_atc",
    }


@network
@pytest.mark.usefixtures("initialize_sliderule")
class TestSlideRule:
    def test_subset_gedi02a(self, tinyaoi, expected_gedi_columns):
        gf_gedi = coincident.search.search(
            dataset="gedi", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.subset_gedi02a(
            gf_gedi, aoi=tinyaoi, include_worldcover=True
        )

        assert "request" in data.attrs["meta"]
        assert isinstance(data, gpd.GeoDataFrame)
        assert set(data.columns) == expected_gedi_columns
        # assert data.iloc[0].name.to_pydatetime().toordinal() == 738192
        # assert data.iloc[0]["worldcover.value"] == "Tree cover"
        assert data.shape == (262, 14)

    def test_include_raster_samples(self, tinyaoi, expected_is2_columns):
        gf_is2 = coincident.search.search(
            dataset="icesat-2", intersects=tinyaoi, datetime="2022"
        )

        data = coincident.io.sliderule.subset_atl06(
            gf_is2, aoi=tinyaoi, include_worldcover=True
        )
        expected_is2_columns.update(
            {"worldcover.time_ns", "worldcover.value", "worldcover.fileid"}
        )

        assert isinstance(data, gpd.GeoDataFrame)
        assert len(data.attrs["meta"]["request"]["resources"]) == 6
        assert len(set(data.attrs["meta"]["srctbl"].values())) == 6
        assert set(data.columns) == expected_is2_columns
        assert {1034, 615} == set(data.rgt)
        # NOTE: not sure order is deterministic...
        assert data.iloc[0].name.toordinal() == 738395
        assert data.shape == (413, 23)

    def test_subset_atl06(self, tinyaoi, expected_is2_columns):
        gf_is2 = coincident.search.search(
            dataset="icesat-2", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.subset_atl06(gf_is2, aoi=tinyaoi)

        assert isinstance(data, gpd.GeoDataFrame)
        assert len(data.attrs["meta"]["request"]["resources"]) == 6
        assert len(set(data.attrs["meta"]["srctbl"].values())) == 6
        assert set(data.columns) == expected_is2_columns
        assert {1034, 615} == set(data.rgt)
        # NOTE: not sure order is deterministic...
        assert data.iloc[0].name.toordinal() == 738395
        assert data.shape == (413, 23)

    @pytest.mark.xfail(reason="https://github.com/uw-cryo/coincident/issues/126")
    def test_process_atl06sr(self, tinyaoi, expected_is2_columns):
        gf_is2 = coincident.search.search(
            dataset="icesat-2", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.process_atl06sr(gf_is2, aoi=tinyaoi)

        assert isinstance(data, gpd.GeoDataFrame)
        assert set(data.columns) == expected_is2_columns
        assert data.shape == (367, 20)
        assert data.iloc[0].name.toordinal() == 738304
        assert int(data.iloc[0].h_mean) == 3178

    @pytest.mark.filterwarnings(
        "ignore:Geometry is in a geographic CRS:UserWarning:geopandas"
    )
    @pytest.mark.xfail(reason="USGS TNM is often down", strict=False)
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
