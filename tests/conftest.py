# ruff: noqa: ARG001
from __future__ import annotations

import geopandas as gpd
import pytest

# import os
# if not os.environ.get('MAXAR_API_KEY'):
#    os.environ['MAXAR_API_KEY'] = 'fake-test-key'


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
