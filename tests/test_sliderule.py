from __future__ import annotations

import os
import pprint

import geopandas as gpd
import numpy as np
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
def expected_atl06SR_columns():
    return {
        "cycle",
        "dh_fit_dx",
        "geometry",
        "gt",
        "h_mean",
        "h_sigma",
        "n_fit_photons",
        "pflags",
        "photon_start",
        "region",
        "rgt",
        "rms_misfit",
        "spot",
        "srcid",
        "w_surface_window_final",
        "x_atc",
        "y_atc",
    }


@pytest.fixture(scope="class")
def expected_atl06_columns():
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
            gf_gedi,
            aoi=tinyaoi,
        )

        assert "request" in data.attrs["meta"]
        assert isinstance(data, gpd.GeoDataFrame)
        assert data.crs == "EPSG:7912"
        assert set(data.columns) == expected_gedi_columns
        assert set(data["track"]) == {11156, 2618, 9477}
        assert set(data["beam"]) == {5, 3, 11, 8, 6}
        assert data.shape == (262, 11)
        assert data.index.min().toordinal() == 738192
        assert int(data["elevation_lm"].min()) == 3041

    # TODO: add test for small area where we span multiple 3DEP UTM zones
    # Can select one or the other via "substr" https://slideruleearth.io/web/rtd/user_guide/raster_sampling.html
    def test_subset_gedi02a_with_raster_sampling(self, tinyaoi):
        gf_gedi = coincident.search.search(
            dataset="gedi", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.subset_gedi02a(
            gf_gedi, aoi=tinyaoi, include_worldcover=True, include_3dep=True
        )

        expected_columns = {
            "3dep.value",
            "3dep.fileid",
            "3dep.time_ns",
            "worldcover.value",
            "worldcover.fileid",
            "worldcover.time_ns",
        }

        threedep_tiffs = {
            value
            for key in data.attrs["meta"]["srcdata"]
            for value in data.attrs["meta"]["srcdata"][key]["samples.3dep"].values()
        }
        # For this small AOI we only sample 1 1m tile
        expected_tif = "/vsis3/prd-tnm/StagedProducts/Elevation/1m/Projects/CO_WestCentral_2019_A19/TIFF/USGS_1M_13_x37y430_CO_WestCentral_2019_A19.tif"

        assert "request" in data.attrs["meta"]
        assert isinstance(data, gpd.GeoDataFrame)
        assert expected_tif in threedep_tiffs
        assert type(data["3dep.value"].iloc[0]) is np.ndarray
        assert data["3dep.time_ns"].iloc[0].dtype == "<M8[ns]"
        assert data["3dep.value"].min().astype(int) == 3041
        assert expected_columns.issubset(set(data.columns))
        assert set(data["worldcover.value"]) == {
            "Grassland",
            "Moss and lichen",
            "Tree cover",
        }
        assert data.shape == (262, 17)

    def test_subset_atl06(self, tinyaoi, expected_atl06_columns):
        gf_is2 = coincident.search.search(
            dataset="icesat-2", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.subset_atl06(gf_is2, aoi=tinyaoi)

        assert isinstance(data, gpd.GeoDataFrame)
        assert data.crs == "EPSG:9989"
        assert len(data.attrs["meta"]["request"]["resources"]) == 6
        assert len(set(data.attrs["meta"]["srctbl"].values())) == 6
        assert set(data.columns) == expected_atl06_columns
        assert {1034, 615} == set(data.rgt)
        assert data.index.min().toordinal() == 738304
        assert data.shape == (405, 23)
        # TODO: investigate Alert <-7>: Failure on resource ATL06_20220131125843_06151402_007_01.h5 beam gt3l: H5Coro::Future read failure on /gt3l/land_ice_segments/latitude
        # assert data.shape == (413, 23)
        assert int(data.h_li.max()) == 3569

    def test_process_atl06sr(self, tinyaoi, expected_atl06SR_columns):
        gf_is2 = coincident.search.search(
            dataset="icesat-2", intersects=tinyaoi, datetime="2022"
        )
        data = coincident.io.sliderule.process_atl06sr(gf_is2, aoi=tinyaoi)

        assert isinstance(data, gpd.GeoDataFrame)
        assert data.crs == "EPSG:9989"
        assert set(data.columns) == expected_atl06SR_columns
        # TODO: difference when updating to sliderule v5/icesat-2 v7 ?...
        # assert data.shape == (367, 20)
        assert data.shape == (494, 17)
        assert data.index.min().toordinal() == 738304
        assert {1034, 615} == set(data.rgt)
        assert int(data.h_mean.max()) == 3568

    @pytest.mark.filterwarnings(
        "ignore:Geometry is in a geographic CRS:UserWarning:geopandas"
    )
    @pytest.mark.xfail(reason="USGS TNM is often down", strict=False)
    def test_sample_raster_3dep_tnm(self):
        gf = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy(
                [-106.46953038438801, -106.47209604357627],
                [38.781089777057986, 38.78292421778072],
                crs="EPSG:9989",
            )
        )
        samples = coincident.io.sliderule.sample_raster(gf)

        assert isinstance(samples, gpd.GeoDataFrame)
        # NOTE: results in same CRS as input
        assert samples.crs == "EPSG:9989"
        assert {"file", "geometry", "time", "value"} == set(samples.columns)
        assert samples.time.dtype == "<M8[ns]"
        assert str(samples.time.dt.date.iloc[0]) == "2019-09-04"
        assert set(samples.value.astype("i2")) == {3110, 3178}
