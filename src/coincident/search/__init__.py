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
    datetime: str | None = None,
    # NOTE: change to explicitly stac_kwargs? not sure of proper kwargs type (Optional[dict[str, Any]]?)
    **kwargs: Any,
) -> _gpd.GeoDataFrame:
    """
    Perform a search for geospatial data using STAC or non-STAC APIs.
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

    if datetime is not None:
        # Use pystac_client to parse all supported datetime inputs
        start, end = _pystac_client._format_datetime(datetime).split(
            "/"
        )  # to tz-aware UTC datetimes
        search_start = _gpd.pd.to_datetime(start)
        search_end = _gpd.pd.to_datetime(end)
    else:
        search_start = search_end = None

    if not isinstance(intersects, _gpd.GeoDataFrame | _gpd.GeoSeries):
        message = f"intersects value must be a GeoDataFrame or GeoSeries, not {type(intersects)}"
        raise ValueError(message)
    if len(intersects) > 1:
        message = "GeoDataFrame contains multiple geometries, search requires a single geometry"
        raise ValueError(message)

    _validate_search_params(dataset, intersects, search_start, search_end)

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
        # Apply datetime interval filter after reading GPKG
        if datetime is not None:
            results_intervals = _gpd.pd.IntervalIndex.from_arrays(
                gf["start_datetime"], gf["end_datetime"], closed="both"
            )
            search_interval = _gpd.pd.Interval(search_start, search_end, closed="both")
            keep = results_intervals.overlaps(search_interval)
            gf = gf.loc[keep, :]

    return gf


def _validate_search_params(
    dataset: _Dataset,
    intersects: _gpd.GeoSeries,
    search_start: _gpd.pd.Timestamp | None = None,
    search_end: _gpd.pd.Timestamp | None = None,
) -> None:
    """
    Validate search parameters including date range and spatial bounds.
    This function checks if the requested search polygon and date range fall within the
    Continental United States (CONUS) and the dataset's date range.
    Parameters
    ----------
    dataset : _Dataset
        The dataset object containing start and end dates.
    intersects : _gpd.GeoDataFrame
        A GeoDataFrame representing the search polygon.
    search_start : _gpd.pd.Timestamp, optional
        The start date of the search range. Default is None.
    search_end : _gpd.pd.Timestamp, optional
        The end date of the search range. Default is None.
    Raises
    ------
    ValueError
        If the search end date is earlier than the dataset start date.
        If the search start date is later than the dataset end date.
        If the intersects parameter is not a GeoDataFrame.
        If the search polygon is not within the Continental U.S. or Alaska.
    """

    # NOTE: currently assumes single dataset searched at a time
    dataset_start = _gpd.pd.to_datetime(dataset.start, utc=True)
    dataset_end = _gpd.pd.to_datetime(dataset.end, utc=True)
    if (search_end is not None) and (search_end < dataset_start):
        message = f"Requested search end {search_end} < {dataset.alias} start date: {dataset_start}"
        raise ValueError(message)
    if (search_start is not None) and (search_start > dataset_end):
        message = f"Requested search start {search_start} > {dataset.alias} end date: {dataset_end}"
        raise ValueError(message)

    # Check spatial bounds
    # NOTE: use geopandas to first convert to projected CRS?
    aoi = intersects.to_crs("EPSG:4326").geometry.iloc[0]  # shapely geometry
    CONUS = _box(*(-124.84, 24.39, -66.88, 49.38))
    AK = _box(*(-179.99, 51.21, -129.63, 71.35))
    if not aoi.within(CONUS) and not aoi.within(AK):
        message = "Requested search polygon not within Continental U.S. or Alaska"
        raise ValueError(message)
        # warnings.warn(f'Requested search polygon not within Continental U.S. or Alaska')
