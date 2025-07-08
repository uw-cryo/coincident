from __future__ import annotations

import pytest
import xarray as xr

import coincident
from coincident.io.xarray import to_dataset

network = pytest.mark.network


@network
def test_to_dataset_with_cop30_lazy(aoi):
    """Test `to_dataset` functionality with COP30 dataset."""
    gf_cop30 = coincident.search.search(
        dataset="cop30",
        intersects=aoi,
    )
    ds = to_dataset(
        gf_cop30,
        aoi=aoi,
        resolution=0.1,  # ~1km
    )
    assert isinstance(ds, xr.Dataset)
    assert "data" in ds.data_vars


@network
def test_to_dataset_with_worldcover_lazy(aoi):
    """Test `to_dataset` functionality with WorldCover dataset."""
    gf_wc = coincident.search.search(
        dataset="worldcover",
        intersects=aoi,
        datetime=["2020"],
    )
    ds = to_dataset(
        gf_wc,
        bands=["map"],
        aoi=aoi,
        resolution=0.1,  # ~1km
    )
    assert isinstance(ds, xr.Dataset)
    assert "map" in ds.data_vars


@network
@pytest.mark.filterwarnings(
    "ignore:Supplying chunks as dimension-order tuples is deprecated.*:DeprecationWarning"
)
def test_load_usgs_dem(usgs_neon_dem_bbox):
    da_usgs_dem = coincident.io.xarray.load_usgs_dem(
        usgs_neon_dem_bbox, "CO_CentralEasternPlains_2020_D20"
    )
    assert da_usgs_dem.name == "elevation"
    assert da_usgs_dem.shape == (28, 32)
    assert da_usgs_dem.dtype == "float32"
    assert "x" in da_usgs_dem.coords
    assert "y" in da_usgs_dem.coords
    assert "spatial_ref" in da_usgs_dem.coords
    assert da_usgs_dem.rio.resolution() == (1.0, -1.0)
    assert da_usgs_dem.notnull().any()


# marks needed for chunks='auto' in rioxarray.open_rasterio() calls
@network
@pytest.mark.filterwarnings(
    "ignore:Supplying chunks as dimension-order tuples is deprecated.*:DeprecationWarning"
)
@pytest.mark.filterwarnings("ignore:TIFFReadDirectoryCheckOrder.*:RuntimeWarning")
def test_load_neon_dem(usgs_neon_dem_bbox):
    da_neon_dem = coincident.io.xarray.load_neon_dem(
        usgs_neon_dem_bbox, datetime_str="2020-06-30", site_id="ARIK", product="dsm"
    )
    assert da_neon_dem.name == "elevation"
    assert da_neon_dem.shape == (28, 32)
    assert da_neon_dem.dtype == "float32"
    assert "x" in da_neon_dem.coords
    assert "y" in da_neon_dem.coords
    assert "spatial_ref" in da_neon_dem.coords
    assert da_neon_dem.rio.resolution() == (1.0, -1.0)
    assert da_neon_dem.notnull().any(), "DataArray contains no valid data"


# NOTE: this test takes 26 seconds on my local laptop to run
# this is because NCALM DEM products are not tiled or COG'd, so no matter how small your AOI is,
# the entire DEM will be loaded into memory. Hoping to find a better use case for this test
@network
def test_load_ncalm_dem(ncalm_dem_bbox):
    da_ncalm_dtm = coincident.io.xarray.load_ncalm_dem(
        ncalm_dem_bbox, product="dtm", dataset_id="WA18_Wall"
    )
    assert da_ncalm_dtm.name == "elevation"
    assert da_ncalm_dtm.shape == (15, 15)
    assert da_ncalm_dtm.dtype == "float32"
    assert "x" in da_ncalm_dtm.coords
    assert "y" in da_ncalm_dtm.coords
    assert "spatial_ref" in da_ncalm_dtm.coords
    assert da_ncalm_dtm.rio.resolution() == (1.0000000000083153, -1.0)
    assert da_ncalm_dtm.notnull().any(), "DataArray contains no valid data"


# NOTE: this test takes 15 seconds on my local laptop to run for the same reason as the NCALM test
# there aren't many smaller NOAA sites but there's always room for improvement
@network
@pytest.mark.filterwarnings(
    "ignore:Supplying chunks as dimension-order tuples is deprecated.*:DeprecationWarning"
)
def test_load_noaa_dem(noaa_dem_bbox):
    da_noaa_dem = coincident.io.xarray.load_noaa_dem(noaa_dem_bbox, 8869)
    assert da_noaa_dem.name == "elevation"
    assert da_noaa_dem.shape == (25, 19)
    assert da_noaa_dem.dtype == "float32"
    assert "x" in da_noaa_dem.coords
    assert "y" in da_noaa_dem.coords
    assert "spatial_ref" in da_noaa_dem.coords
    assert da_noaa_dem.rio.resolution() == (3.0, -3.0)
    assert da_noaa_dem.notnull().any(), "DataArray contains no valid data"


@network
def test_load_gliht_raster():
    import os

    if not os.getenv("EARTHDATA_USER") or not os.getenv("EARTHDATA_PASS"):
        pytest.skip("EARTHDATA_USER and EARTHDATA_PASS environment variables required")
    import geopandas as gpd
    from shapely.geometry import box

    mini_aoi = gpd.GeoDataFrame(
        geometry=[box(*[-76.55240202, 38.8885878, -76.54340202, 38.8975878])],
        crs="EPSG:4326",
    )
    da_chm = coincident.io.xarray.load_gliht_raster(
        aoi=mini_aoi,
        dataset_id="GLMETRICS_SERC_CalTarps_31July2017_am_l0s0",
        product="chm",
    )
    assert da_chm.shape == (1010, 637)
    assert da_chm.dtype == "float32"
    assert "x" in da_chm.coords
    assert "y" in da_chm.coords
    assert "spatial_ref" in da_chm.coords
    assert da_chm.rio.resolution() == (1.0, -1.0)
    assert da_chm.rio.crs.to_epsg() == 32618
    assert da_chm.notnull().any(), "DataArray contains no valid data"
