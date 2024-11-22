"""
Subset and Sample using Sliderule
https://slideruleearth.io
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import geopandas as gpd

try:
    from sliderule import earthdata, gedi, icesat2, raster, toregion
except ImportError:
    warnings.warn(
        "'sliderule' package not found. Install for GEDI & ICESat2 functionality: https://slideruleearth.io/web/rtd/getting_started/Install.html",
        stacklevel=2,
    )
from coincident._utils import depends_on_optional


def _gdf_to_sliderule_polygon(gf: gpd.GeoDataFrame) -> list[dict[str, float]]:
    # Ignore type necessary I think b/c sliderule doesn't have type hints?
    return toregion(gf[["geometry"]].dissolve())["poly"]  # type: ignore[no-any-return]


def _gdf_to_sliderule_params(gf: gpd.GeoDataFrame) -> dict[str, Any]:
    # Sliderule expects single 'geometry' column
    # Dissolve creates a single multipolygon from all geometries
    sliderule_aoi = _gdf_to_sliderule_polygon(gf)
    params = {}
    params["poly"] = sliderule_aoi
    params["t0"] = gf.start_datetime.min().tz_localize(None).isoformat()
    params["t1"] = gf.end_datetime.max().tz_localize(None).isoformat()
    return params


def _granule_from_assets(assets: gpd.GeoDataFrame) -> str:
    # NOTE: change to while loop in case tons of assets?
    for _k, v in assets.items():
        if v.get("roles") == "data":
            granule = Path(v.get("href")).name

    return granule


@depends_on_optional("sliderule")
def subset_gedi02a(
    gf: gpd.GeoDataFrame,
    aoi: gpd.GeoDataFrame = None,
    # NOTE: investigate B006 best-practice
    sliderule_params: dict[str, Any] = {},  # noqa: B006
) -> gpd.GeoDataFrame:
    """
    Subsets GEDI L2A data based on a given GeoDataFrame and optional area of interest (AOI).

    Parameters
    ----------
    gf : gpd.GeoDataFrame
        A GeoDataFrame of results from coincident.search.search(dataset='gedi')

    aoi : gpd.GeoDataFrame, optional
        A GeoDataFrame with a POLYGON to subset. If not provided the union of geometries in gf is used.

    sliderule_params : dict[Any], optional
        A dictionary of additional parameters to be passed to the sliderule `gedi.gedi02ap`.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the subsetted GEDI L2A data.

    Notes
    -----
    The function sets `degrade_filter=True` and `l2_quality_filter=True` by default.
    """
    params = _gdf_to_sliderule_params(gf)

    granule_names = gf.assets.apply(_granule_from_assets).to_list()

    if aoi is not None:
        params["poly"] = _gdf_to_sliderule_polygon(aoi)

    params.update(
        {
            "degrade_filter": True,
            "l2_quality_filter": True,
        }
    )

    # User-provided parameters take precedence
    params.update(sliderule_params)

    return gedi.gedi02ap(params, resources=granule_names)


@depends_on_optional("sliderule")
def subset_atl06(
    gf: gpd.GeoDataFrame,
    aoi: gpd.GeoDataFrame = None,
    dropna: bool = True,
    sliderule_params: dict[str, Any] = {},  # noqa: B006
) -> gpd.GeoDataFrame:
    params = _gdf_to_sliderule_params(gf)

    # Note necessary but avoids another CMR search
    granule_names = gf.assets.apply(_granule_from_assets).to_list()

    if aoi is not None:
        params["poly"] = _gdf_to_sliderule_polygon(aoi)

    # User-provided parameters take precedence
    params.update(sliderule_params)

    data = icesat2.atl06sp(params, resources=granule_names)

    # NOTE: add server-side filtering for NaNs in sliderule.atl06sp?
    if dropna:
        return data.dropna(subset=["h_li"])
    return data


@depends_on_optional("sliderule")
def sample_worldcover(gf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    points = [[x, y] for x, y in zip(gf.geometry.x, gf.geometry.y, strict=False)]
    return raster.sample("esa-worldcover-10meter", points)


@depends_on_optional("sliderule")
def sample_3dep(gf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Sample 3DEP 1-meter DEM with POINT geometries in GeoDataFrame.

    Points should be (longitude, latitude) not UTM

    Parameters
    ----------
    gf : gpd.GeoDataFrame
        A GeoDataFrame containing POINT geometries

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame with sampled elevation data from the 3DEP 1-meter DEM.
    """
    points = [[x, y] for x, y in zip(gf.geometry.x, gf.geometry.y, strict=False)]

    poly = _gdf_to_sliderule_polygon(gf)
    geojson = earthdata.tnm(
        short_name="Digital Elevation Model (DEM) 1 meter", polygon=poly
    )

    # NOTE: make sliderule aware of 3dep polygons to not sample NaN regions for 1m grid?
    data_3dep = raster.sample("usgs3dep-1meter-dem", points, {"catalog": geojson})

    return data_3dep[data_3dep.value.notna()]
