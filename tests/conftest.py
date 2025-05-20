# ruff: noqa: ARG001
from __future__ import annotations

import geopandas as gpd
import pytest
import xarray as xr
from shapely.geometry import box


@pytest.fixture(scope="package")
def aoi():
    # 11 vertices, 1,361km^2
    aoi_url = "https://raw.githubusercontent.com/SlideRuleEarth/sliderule-python/main/data/grandmesa.geojson"
    return gpd.read_file(aoi_url)


@pytest.fixture
def large_aoi(scope="package"):
    # 260 vertices, large area 269,590 km^2
    aoi_url = "https://raw.githubusercontent.com/unitedstates/districts/refs/heads/gh-pages/states/CO/shape.geojson"
    return gpd.read_file(aoi_url)


# TODO: add a small geojson AOI for testing bathy sites
@pytest.fixture
def bathy_aoi(scope="package"):
    # 631 vertices, large area ~ 27,000 km^2
    aoi_url = "https://raw.githubusercontent.com/unitedstates/districts/refs/heads/gh-pages/states/MA/shape.geojson"
    return gpd.read_file(aoi_url)


@pytest.fixture(scope="package")
def tinyaoi():
    # 4 vertices
    return gpd.read_file("tests/data/tiny.geojson")


@pytest.fixture(scope="package")
def dem_tiny():
    # from CO PCD site
    ds_cop30 = xr.open_dataset("tests/data/dem_tiny.tif")
    return ds_cop30.rename({"band_data": "elevation"})


@pytest.fixture(scope="package")
def points_tiny():
    # random points as is2/gedi facade overlapping and same CRS as dem_tiny
    # has cols h_li, elevation_lm, and elevation_hr
    # ^ random elevation values
    return gpd.read_file("tests/data/points_tiny.gpkg")


@pytest.fixture(scope="package")
def dem_tiny_utm():
    ds_cop30 = xr.open_dataset("tests/data/dem_tiny_utm.tif")
    return ds_cop30.rename({"band_data": "elevation"})


@pytest.fixture(scope="package")
def points_tiny_utm():
    return gpd.read_file("tests/data/points_tiny_utm.gpkg")


@pytest.fixture(scope="package")
def usgs_neon_dem_bbox():
    bbox_geometry = box(
        -102.51096149997586, 39.673538237665866, -102.51058942897973, 39.67378984946736
    )
    return gpd.GeoDataFrame(geometry=[bbox_geometry], crs="EPSG:4326")


@pytest.fixture(scope="package")
def ncalm_dem_bbox():
    bbox_geometry = box(-121.39282517, 46.50320809, -121.39262917, 46.50334309)
    return gpd.GeoDataFrame(geometry=[bbox_geometry], crs="EPSG:4326")


@pytest.fixture(scope="package")
def noaa_dem_bbox():
    bbox_geometry = box(-122.50657366, 42.07633832, -122.50637366, 42.07653832)
    return gpd.GeoDataFrame(geometry=[bbox_geometry], crs="EPSG:4326")
