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
        if os.environ.get("PS_GITHUB_TOKEN"):
            message = (
                "PS_GITHUB_TOKEN is set, so using UW Cluster instead of public cluster"
            )
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


@network
@pytest.mark.usefixtures("initialize_sliderule")
class TestSlideRule:
    def test_subset_gedi02a(self, tinyaoi):
        gf_gedi = coincident.search.search(
            dataset="gedi", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.subset_gedi02a(
            gf_gedi, aoi=tinyaoi, include_worldcover=True
        )

        expected_columns = {
            "orbit",
            "beam",
            "elevation_lm",
            "solar_elevation",
            "elevation_hr",
            "sensitivity",
            "track",
            "flags",
            "geometry",
            "worldcover.flags",
            "worldcover.file_id",
            "worldcover.value",
            "worldcover.time",
        }

        assert isinstance(data, gpd.GeoDataFrame)
        assert set(data.columns) == expected_columns
        assert data.iloc[0].name.toordinal() == 738192
        assert data.iloc[0]["worldcover.value"] == "Tree cover"
        assert data.shape == (262, 13)

    def test_subset_atl06(self, tinyaoi):
        gf_is2 = coincident.search.search(
            dataset="icesat-2", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.subset_atl06(gf_is2, aoi=tinyaoi)

        expected_columns = {
            "bsnow_h",
            "w_surface_window_final",
            "h_robust_sprd",
            "r_eff",
            "n_fit_photons",
            "sigma_geo_h",
            "segment_id",
            "gt",
            "cycle",
            "y_atc",
            "atl06_quality_summary",
            "spot",
            "dh_fit_dx",
            "h_li_sigma",
            "seg_azimuth",
            "h_li",
            "bsnow_conf",
            "rgt",
            "x_atc",
            "tide_ocean",
            "geometry",
        }

        assert isinstance(data, gpd.GeoDataFrame)
        assert set(data.columns) == expected_columns
        assert 1034 in data.rgt.unique()
        assert data.iloc[0].name.toordinal() == 738304
        assert data.shape == (413, 21)

    @pytest.mark.xfail(reason="https://github.com/uw-cryo/coincident/issues/126")
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
