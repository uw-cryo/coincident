"""
Wrapper for various data provider search APIs

Results returned as GeoDataFrame with consistent STAC column headings and best dtypes
"""

from __future__ import annotations

from typing import Any

import geopandas as _gpd

# Initialize a fake pystac_client search for datetime parsing
from pystac_client.item_search import ItemSearch as _ItemSearch
from shapely.geometry import box as _box

from coincident.datasets import _alias_to_Dataset
from coincident.datasets.general import Dataset as _Dataset
from coincident.search import _stac, _wesm

_pystac_client = _ItemSearch("no_url")


def search(
    dataset: str | _Dataset,
    intersects: _gpd.GeoDataFrame | _gpd.GeoSeries,
    datetime: str | list[str] | None = None,
    # NOTE: change to explicitly stac_kwargs? not sure of proper kwargs type (Optional[dict[str, Any]]?)
    **kwargs: Any,
) -> _gpd.GeoDataFrame:
    """
    Perform a search for geospatial data using STAC or non-STAC APIs for a single dataset.
    Parameters
    ----------
    dataset : str or _Dataset
        The dataset to search. Can be a string alias or a _Dataset object.
    intersects : _gpd.GeoDataFrame
        A GeoDataFrame containing a single row, or GeoSeries containing a geometry to restrict the search area.
    datetime : str or None, optional
        The datetime range for the search in ISO 8601 format. If None, no datetime filter is applied.
    **kwargs : Any
        Additional keyword arguments to pass to the search functions.
        see https://pystac-client.readthedocs.io/en/latest/api.html#item-search
    Returns
    -------
    _gpd.GeoDataFrame
        A GeoDataFrame containing the search results.
    Raises
    ------
    ValueError
        If the provided dataset alias is not supported.
    """
    # Validation of inputs
    if isinstance(dataset, str):
        try:
            dataset = _alias_to_Dataset[dataset]
        except KeyError as e:
            message = (
                f"{dataset} is not a supported dataset: {_alias_to_Dataset.keys()}"
            )
            raise ValueError(message) from e

    _validate_temporal_bounds(dataset, datetime)
    _validate_spatial_bounds(intersects)

    # NOTE: not very robust, explode() demotes MultiPolygons to single Polygon (seems many GeoJSONs have this)
    # ANd 'exterior' not available for Multipolygons, just
    shapely_geometry = intersects.geometry.explode().iloc[0]
    if (
        not shapely_geometry.exterior.is_ccw
    ):  # Apparently NASA CMR enforces polygon CCW order
        shapely_geometry = shapely_geometry.reverse()
    aoi = _pystac_client._format_intersects(shapely_geometry)  # to JSON geometry

    # STAC API Searches
    if dataset.has_stac_api:
        stac_api_kwargs = {**kwargs, **dataset.stac_kwargs}  # NOTE: override order?
        stac_api_kwargs["datetime"] = datetime
        stac_api_kwargs["collections"] = dataset.collections
        stac_api_kwargs["intersects"] = aoi

        if dataset.alias == "maxar":
            # NOTE: not sure how to avoid incompatible type "str | None"; expected "str" for Dataset.attrs
            client = _stac.configure_maxar_client(dataset.area_based_calc)  # type: ignore[attr-defined]
            results = _stac.search(client, **stac_api_kwargs)
            gf = _stac.to_geopandas(results)
            # Client-side reduce to only acquisitions having stereo pairs
            gf = gf.loc[gf.stereo_pair_identifiers.str[0].dropna().index]

        elif dataset.alias in ["icesat-2", "gedi"]:
            client = _stac.configure_nasa_client(dataset.search)  # type: ignore[arg-type]
            results = _stac.search(client, **stac_api_kwargs)
            gf = _stac.to_geopandas(results)

        elif dataset.alias in ["cop30"]:
            client = _stac.configure_mspc_client(dataset.search)  # type: ignore[arg-type]
            results = _stac.search(client, **stac_api_kwargs)
            gf = _stac.to_geopandas(results)

    # Non-STAC Searches
    elif dataset.alias == "3dep":
        gf = _wesm.search_gpkg(
            dataset.search,  # type: ignore[arg-type]
            mask=intersects,
            **kwargs,
        )
        # NOTE: Apply datetime interval filter after reading GPKG?

    return gf


def _validate_temporal_bounds(
    dataset: _Dataset,
    datetime: str | list[str] | None = None,
) -> None:
    """
    Validates the temporal bounds of a dataset against a given datetime range.
    Parameters
    ----------
    dataset : _Dataset
        The dataset object containing start and end attributes.
    datetime : str, list, or None, optional
        The datetime range to validate against the dataset's temporal bounds.
        It can be a single string, a list of strings, or None. The datetime
        range should be in a format supported by pystac_client.
    Raises
    ------
    ValueError
        If the requested search start date is after the dataset's end date or
        if the requested search end date is before the dataset's start date.
    """
    if datetime is not None:
        # Use pystac_client to parse all supported datetime inputs
        start, end = _pystac_client._format_datetime(datetime).split("/")

        if start != ".." and dataset.end is not None:
            search_start = _gpd.pd.to_datetime(start)
            dataset_end = _gpd.pd.to_datetime(dataset.end, utc=True)
            if search_start > dataset_end:
                message = f"Requested search start {search_start} > {dataset.alias} end date: {dataset_end}"
                raise ValueError(message)

        if end != ".." and dataset.start is not None:
            search_end = _gpd.pd.to_datetime(end)
            dataset_start = _gpd.pd.to_datetime(dataset.start, utc=True)
            if search_end < dataset_start:
                message = f"Requested search end {search_end} < {dataset.alias} start date: {dataset_start}"
                raise ValueError(message)


def _validate_spatial_bounds(
    intersects: _gpd.GeoSeries,
) -> None:
    """
    Validate that the specified area of interest (AOI) within CONUS or Alaska.
    Parameters
    ----------
    dataset : _Dataset
        The dataset to validate.
    intersects : _gpd.GeoSeries
        A GeoSeries containing the area of interest (AOI) geometry.
    Raises
    ------
    ValueError
        If the AOI is not within the bounds of the Continental U.S. or Alaska.
    """
    if not isinstance(intersects, _gpd.GeoDataFrame | _gpd.GeoSeries):
        message = f"intersects value must be a GeoDataFrame or GeoSeries, not {type(intersects)}"
        raise ValueError(message)
    if len(intersects) > 1:
        message = "GeoDataFrame contains multiple geometries, search requires a single geometry"
        raise ValueError(message)

    # NOTE: use geopandas to first convert to projected CRS?
    aoi = intersects.to_crs("EPSG:4326").geometry.iloc[0]  # shapely geometry
    CONUS = _box(*(-124.84, 24.39, -66.88, 49.38))
    AK = _box(*(-179.99, 51.21, -129.63, 71.35))
    if not aoi.within(CONUS) and not aoi.within(AK):
        message = "Requested search polygon not within Continental U.S. or Alaska"
        raise ValueError(message)
        # warnings.warn(f'Requested search polygon not within Continental U.S. or Alaska')
