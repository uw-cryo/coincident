from __future__ import annotations

import geopandas as gpd
import pytest

import coincident.overlaps


@pytest.fixture
def geodataframe():
    """1x1 degree ll corner at (0,0)"""
    fc = {
        "type": "FeatureCollection",
        "features": [
            {
                "id": "0",
                "properties": {
                    "start_datetime": gpd.pd.Timestamp("2021-01-01"),
                    "end_datetime": gpd.pd.Timestamp("2021-01-06"),
                },
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [0.5, -0.5],
                            [0.5, 0.5],
                            [-0.5, 0.5],
                            [-0.5, -0.5],
                            [0.5, -0.5],
                        ]
                    ],
                },
            }
        ],
    }
    return gpd.GeoDataFrame.from_features(fc, crs="EPSG:4326")


def test_geographic_area(geodataframe):
    area = (
        (coincident.overlaps.geographic_area(geodataframe) * 1e-6)
        .to_numpy()
        .astype(int)
    )
    assert area == 12309


def test_subset_by_minimum_area(geodataframe):
    big_enough = coincident.overlaps.subset_by_minimum_area(geodataframe)
    too_small = coincident.overlaps.subset_by_minimum_area(geodataframe, 20100)
    assert len(big_enough) == 1
    assert len(too_small) == 0


def test_subset_by_temporal_overlap(geodataframe):
    start = gpd.pd.Timestamp("2020-12-01")
    end = gpd.pd.Timestamp("2020-12-31")
    nonoverlapping = coincident.overlaps.subset_by_temporal_overlap(
        geodataframe, start_datetime=start, end_datetime=end, temporal_buffer=0
    )
    overlapping = coincident.overlaps.subset_by_temporal_overlap(
        geodataframe, start_datetime=start, end_datetime=end, temporal_buffer=5
    )
    assert len(nonoverlapping) == 0
    assert len(overlapping) == 1
