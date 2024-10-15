from __future__ import annotations

from typing import Any

import geopandas as gpd

# Used to access formatters
from pystac_client.item_search import ItemSearch as _ItemSearch
from shapely.geometry import box

from coincident.datasets import _alias_to_Dataset
from coincident.datasets.general import Dataset
from coincident.search import stac, wesm

_pystac_client = _ItemSearch("no_url")


def search(
    dataset: str | Dataset,
    intersects: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    datetime: str | list[str] | None = None,
    # NOTE: change to explicitly stac_kwargs? not sure of proper kwargs type (Optional[dict[str, Any]]?)
    **kwargs: Any,
) -> gpd.GeoDataFrame:
    """
    Perform a search for geospatial data using STAC or non-STAC APIs for a single dataset.

    Parameters
    ----------
    dataset : str or Dataset
        The dataset to search. Can be a string alias or a Dataset object.
    intersects : gpd.GeoDataFrame
        A GeoDataFrame containing a single row, or GeoSeries containing a geometry to restrict the search area.
    datetime : str or None, optional
        The datetime range for the search in ISO 8601 format. If None, no datetime filter is applied.
    **kwargs : Any
        Additional keyword arguments to pass to the search functions.
        see https://pystac-client.readthedocs.io/en/latest/api.html#item-search

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the search results.

    Raises
    ------
    ValueError
        If the provided dataset alias is not supported.
    """
    # Validate Dataset
    if isinstance(dataset, str):
        try:
            dataset = _alias_to_Dataset[dataset]
        except KeyError as e:
            message = (
                f"{dataset} is not a supported dataset: {_alias_to_Dataset.keys()}"
            )
            raise ValueError(message) from e

    # Validate Datetimes
    _validate_temporal_bounds(dataset, datetime)
    if datetime is not None:
        start, end = _pystac_client._format_datetime(datetime).split("/")
        search_start = (
            gpd.pd.to_datetime(start).tz_localize(None) if start != ".." else None
        )
        search_end = gpd.pd.to_datetime(end).tz_localize(None) if end != ".." else None
    else:
        search_start = None  # or  gpd.pd.Timestamp(dataset.start) ?
        search_end = None  # or gpd.pd.Timestamp.today()?

    # Validate Intersects
    if intersects is not None:
        _validate_spatial_bounds(intersects)

        # NOTE: not very robust, explode() demotes MultiPolygons to single Polygon (seems many GeoJSONs have this)
        # ANd 'exterior' not available for Multipolygons, just
        shapely_geometry = intersects.geometry.explode().iloc[0]
        if (
            not shapely_geometry.exterior.is_ccw
        ):  # Apparently NASA CMR enforces polygon CCW order
            shapely_geometry = shapely_geometry.reverse()
        aoi = _pystac_client._format_intersects(shapely_geometry)  # to JSON geometry
    else:
        aoi = None

    # STAC API Searches
    # NOTE: refactor so that search config lives alongside dataset? easier to add, would need conditionals below
    if dataset.has_stac_api:
        stac_api_kwargs = {**kwargs, **dataset.stac_kwargs}  # NOTE: override order?
        stac_api_kwargs["datetime"] = datetime
        stac_api_kwargs["collections"] = dataset.collections
        stac_api_kwargs["intersects"] = aoi

        if dataset.provider == "maxar":
            # NOTE: not sure how to avoid incompatible type "str | None"; expected "str" for Dataset.attrs
            client = stac.configure_maxar_client(dataset.area_based_calc)  # type: ignore[attr-defined]
            results = stac.search(client, **stac_api_kwargs)
            gf = stac.to_geopandas(results)
            # Client-side reduce to only acquisitions having stereo pairs
            gf = gf.loc[gf.stereo_pair_identifiers.str[0].dropna().index]

        else:
            if dataset.provider == "microsoft":
                client = stac.configure_mspc_client(dataset.search)  # type: ignore[arg-type]

            # NOTE: NASA-CMR-STAC seems to require GET, but maxar prefers POST (errors for large polygons)?!
            elif dataset.provider == "nasa":
                stac_api_kwargs["method"] = "GET"
                client = stac.configure_stac_client(dataset.search)  # type: ignore[arg-type]
            # Generic STAC endpoint w/o additional config
            else:
                client = stac.configure_stac_client(dataset.search)  # type: ignore[arg-type]
            results = stac.search(client, **stac_api_kwargs)
            gf = stac.to_geopandas(results)

    # Non-STAC Searches
    elif dataset.alias == "3dep":
        gf = wesm.search_convex_hulls(
            intersects=intersects,
            search_start=search_start,
            search_end=search_end,
            **kwargs,
        )

    return gf


def _validate_temporal_bounds(
    dataset: Dataset,
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
            search_start = gpd.pd.to_datetime(start)
            dataset_end = gpd.pd.to_datetime(dataset.end, utc=True)
            if search_start > dataset_end:
                message = f"Requested search start {search_start} > {dataset.alias} end date: {dataset_end}"
                raise ValueError(message)

        if end != ".." and dataset.start is not None:
            search_end = gpd.pd.to_datetime(end)
            dataset_start = gpd.pd.to_datetime(dataset.start, utc=True)
            if search_end < dataset_start:
                message = f"Requested search end {search_end} < {dataset.alias} start date: {dataset_start}"
                raise ValueError(message)


def _validate_spatial_bounds(
    intersects: gpd.GeoSeries,
) -> None:
    """
    Validate that the specified area of interest (AOI) wis of type GeoDataFrame or GeoSeries
    ----------
    dataset : _Dataset
        The dataset to validate.
    intersects : gpd.GeoSeries
        A GeoSeries containing the area of interest (AOI) geometry.
    Raises
    ------
    ValueError
        If the AOI is not type GeoDataFrame or GeoSeries
    """
    if not isinstance(intersects, gpd.GeoDataFrame | gpd.GeoSeries):
        message = f"intersects value must be a GeoDataFrame or GeoSeries, not {type(intersects)}"
        raise ValueError(message)
    if len(intersects) > 1:
        message = "GeoDataFrame contains multiple geometries, search requires a single geometry"
        raise ValueError(message)

    # NOTE: use geopandas to first convert to projected CRS?
    #aoi = intersects.to_crs("EPSG:4326").geometry.iloc[0]  # shapely geometry
    #CONUS = box(*(-124.84, 24.39, -66.88, 49.38))
    #AK = box(*(-179.99, 51.21, -129.63, 71.35))
    #if not aoi.within(CONUS) and not aoi.within(AK):
    #    message = "Requested search polygon not within Continental U.S. or Alaska"
    #    raise ValueError(message)
        # warnings.warn(f'Requested search polygon not within Continental U.S. or Alaska')
