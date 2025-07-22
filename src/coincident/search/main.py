from __future__ import annotations

import warnings
from typing import Any

import geopandas as gpd

# Used to access formatters
from pystac_client.item_search import ItemSearch as _ItemSearch

from coincident.datasets import _alias_to_Dataset
from coincident.datasets.general import Dataset
from coincident.overlaps import subset_by_minimum_area
from coincident.search import neon_api, opentopo_api, stac, wesm

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
            msg_unsupported = (
                f"{dataset} is not a supported dataset: {_alias_to_Dataset.keys()}"
            )
            raise ValueError(msg_unsupported) from e

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
        # NOTE: force_2d as some STAC searches fail with 3D polygons
        # https://github.com/uw-cryo/coincident/issues/101#issuecomment-3104277451
        shapely_geometry = intersects.geometry.force_2d().explode().iloc[0]

        if not shapely_geometry.exterior.is_ccw:
            shapely_geometry = (
                shapely_geometry.reverse()
            )  # Apparently NASA CMR enforces polygon CCW order

        aoi = _pystac_client._format_intersects(shapely_geometry)  # to JSON geometry

    else:
        if "bbox" not in kwargs and "ids" not in kwargs:
            msg_unconstrained = (
                "Neither `bbox` nor `intersects` provided... search will be global"
            )
            warnings.warn(msg_unconstrained, stacklevel=2)
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
            item_collection = stac.search(client, **stac_api_kwargs)
            gf = stac.to_geopandas(item_collection)
            # Client-side reduce to only acquisitions having stereo pairs
            gf = gf.loc[gf.stereo_pair_identifiers.str[0].dropna().index]

        else:
            # NOTE: NASA-CMR-STAC seems to require GET, but default is POST (maxar API w/ GET throws error for large polygons)
            # NOTE: In July 2025, ICESat-2 ATL03 endpoint moved from SIDC_ECS to NSIDC_CPRD, this works with POST
            # if (dataset.provider == "nasa") & (dataset.alias != "icesat-2"):
            if dataset.provider == "nasa":
                stac_api_kwargs["method"] = "GET"
                client = stac.configure_stac_client(dataset.search)  # type: ignore[arg-type]
            # Generic STAC endpoint w/o additional config
            else:
                client = stac.configure_stac_client(dataset.search)  # type: ignore[arg-type]
            item_collection = stac.search(client, **stac_api_kwargs)

            # Per-dataset munging
            # https://github.com/uw-cryo/coincident/issues/8#issuecomment-2449810481
            if dataset.alias == "tdx":
                # Drop columns with messy schema
                dropcols = [
                    "sceneInfo",
                    "missionInfo",
                    "previewInfo",
                    "imageDataInfo",
                    "generationInfo",
                    "acquisitionInfo",
                    "productVariantInfo",
                ]
                for item in item_collection:
                    for col in dropcols:
                        item.properties.pop(col)

            gf = stac.to_geopandas(item_collection)

    # Non-STAC Searches
    elif dataset.alias == "3dep":
        gf = wesm.search_bboxes(
            intersects=intersects,
            search_start=search_start,
            search_end=search_end,
            **kwargs,
        )

    elif dataset.alias == "neon":
        gf = neon_api.search_bboxes(
            intersects=intersects,
            search_start=search_start,
            search_end=search_end,
        )

    elif dataset.provider == "opentopography":
        # TODO: best way to address lpc vs dem argument in opentopo_api.search_ncalm_noaa()?
        gf = opentopo_api.search_opentopo(
            intersects=intersects,
            search_start=search_start,
            search_end=search_end,
            dataset=dataset.alias,
        )

    # Keep track of dataset alias in geodataframe metadata
    # NOTE: attrs not always retrained and not saved to file
    # https://github.com/geopandas/geopandas/issues/3320
    gf.attrs["dataset"] = dataset.alias

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
    Validate that the specified area of interest (AOI) is type GeoDataFrame or GeoSeries
    Parameters
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


def cascading_search(
    primary_dataset: gpd.GeoDataFrame,
    secondary_datasets: list[tuple[str, int]] = [  # noqa: B006
        ("maxar", 14),
        ("icesat", 40),
        ("gedi", 40),
    ],
    min_overlap_area: float = 20,
) -> list[gpd.GeoDataFrame]:
    """
    Perform an cascading search to find overlapping datasets acquired within specific time ranges.

    Secondary datasets are searched based only on spatial overlap areas with previous datasets. In other words, the overlapping area is progressively reduced.

    Temporal buffer is applied as either (datetime-buffer <= acquisition <= datetime+buffer)
     or (start_datetime-buffer <= acquisition <= end_datetime+buffer)

    Parameters
    ----------
    primary_dataset : gpd.GeoDataFrame
        The primary dataset having 'datetime' or 'start_dateteime' and 'end_datetime' columns.

    secondary_datasets : list of tuple, optional
        Each tuple contains the name of the secondary dataset and temporal buffer in days.

    min_overlap_area : int, optional
        The minimum overlap area in km^2. Default is 20.

    Returns
    -------
    list of gpd.GeoDataFrame
        A list of GeoDataFrames containing the search results for each secondary dataset.
    """
    # Do searches on simple geometry, but intersect results with original geometry
    search_geometry = primary_dataset.simplify(0.01)  # or convex_hull?
    detailed_geometry = primary_dataset[["geometry"]]

    if "end_datetime" in primary_dataset.columns:
        start = primary_dataset.start_datetime.iloc[0]
        end = primary_dataset.end_datetime.iloc[0]
    else:
        start = end = primary_dataset.datetime.iloc[0]

    results = []
    for dataset, temporal_buffer in secondary_datasets:
        pad = gpd.pd.Timedelta(days=temporal_buffer)
        date_range = [start - pad, end + pad]

        # Search secondary dataset
        gfs = search(
            dataset=dataset,
            intersects=search_geometry,
            datetime=date_range,
        )

        if dataset == "maxar":
            gfs["stereo_pair_id"] = gfs.stereo_pair_identifiers.str[0]
            gfs = gfs.dissolve(by="stereo_pair_id", as_index=False)

        # Keep track of original footprints
        gfs["original_geometry"] = gfs["geometry"]

        gf_i = gfs.overlay(detailed_geometry, how="intersection")
        gf_i = subset_by_minimum_area(gf_i, min_overlap_area)
        results.append(gf_i)

        # We've refined our search polygon again, so update search GeoDataFrame
        detailed_geometry = gf_i.dissolve()[["geometry"]]
        search_geometry = detailed_geometry.simplify(0.01)

    return results
