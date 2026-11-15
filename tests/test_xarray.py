from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import xarray as xr
from shapely.geometry import box

import coincident
from coincident._utils import vantor_api_key
from coincident.io.xarray import to_dataset

network = pytest.mark.network


def test_mosaic_max_adjacent_tiles():
    # Tiles with different x coords; pytest turns warnings into errors, so this
    # also checks the xarray join FutureWarning isn't emitted
    y = [1.0, 0.0]
    left = xr.DataArray(np.full((2, 3), 1.0), coords={"y": y, "x": [0.0, 1.0, 2.0]})
    right = xr.DataArray(np.full((2, 3), 2.0), coords={"y": y, "x": [2.0, 3.0, 4.0]})

    merged = coincident.io.xarray._mosaic_max([left, right])

    assert merged.x.values.tolist() == [0.0, 1.0, 2.0, 3.0, 4.0]
    # overlapping column x=2 takes the max; no NaNs in the mosaic
    np.testing.assert_array_equal(merged.sel(x=2.0).values, [2.0, 2.0])
    assert not merged.isnull().any()


vantor_authenticated = pytest.mark.skipif(
    not vantor_api_key(), reason="VANTOR_API_KEY not set"
)


@network
@vantor_authenticated
def test_open_vantor_browse():
    """Browse images require the API key to be passed as an HTTP header."""
    gf = coincident.search.search(dataset="vantor", ids=["102001008EC5AC00"])
    item = coincident.search.to_pystac_items(gf)[0]
    da = coincident.io.xarray.open_vantor_browse(item, overview_level=2)
    assert isinstance(da, xr.DataArray)
    assert da.dtype == "uint8"


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
@pytest.mark.xfail(reason="USGS TNM is often down", strict=False)
def test_load_usgs_dem(usgs_neon_dem_bbox):
    da_usgs_dem = coincident.io.xarray.load_usgs_dem(
        usgs_neon_dem_bbox, "CO_CentralEasternPlains_2020_D20"
    )
    assert da_usgs_dem.name == "elevation"
    assert da_usgs_dem.shape == (30, 34)
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
@pytest.mark.skip(reason="temporarily skip b/c very slow :(")
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
def test_load_gliht():
    mini_aoi = gpd.GeoDataFrame(
        geometry=[box(*[-76.55240202, 38.8885878, -76.54340202, 38.8975878])],
        crs="EPSG:4326",
    )
    gf = coincident.search.search(
        dataset="gliht", ids=["GLCHMT_SERC_CalTarps_31July2017_am_l0s0_CHM"]
    )

    da_chm = coincident.io.xarray.load_gliht(
        item=gf.iloc[0],
        aoi=mini_aoi,
    )

    assert da_chm.shape == (1, 1013, 637)
    assert da_chm.dtype == "float32"
    assert {"x", "y", "time", "spatial_ref"} == set(da_chm.coords)
    assert da_chm.rio.resolution() == (1.0, -1.0)
    assert da_chm.rio.crs.to_epsg() == 32618
