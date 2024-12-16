"""
Searching for NCALM and NOAA data with an input search polygon is made easy via the OpenTopography API
https://portal.opentopography.org/apidocs/
There is also a GitHub repo containing flight extent polygons for NOAA
https://github.com/OpenTopography/Data_Catalog_Spatial_Boundaries/tree/main/NOAA_Coastal_Lidar
"""

from __future__ import annotations

from datetime import datetime

import geopandas as gpd
import requests
from shapely.geometry import shape


def search_ncalm_noaa(
    aoi: gpd.GeoDataFrame | gpd.GeoSeries,
    start: str,
    end: str,
    dataset: str | None,  # Default dataset is noaa for now
) -> gpd.GeoDataFrame:
    """
    Perform a search for geospatial data using OpenTopography API.
    This function dynamically adjusts the API URL based on the dataset.

    Parameters
    ----------
    aoi : gpd.GeoDataFrame | gpd.GeoSeries
        A GeoDataFrame or GeoSeries containing a geometry to restrict the search area.
    start : str
        The start datetime in ISO 8601 format (e.g., "2018-01-01").
    end : str
        The end datetime in ISO 8601 format (e.g., "2024-12-31").
    dataset : str
        The dataset type (either "noaa" or "ncalm").
    **kwargs : Any
        Additional keyword arguments to pass to the OpenTopography API.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the search results.

    Raises
    ------
    ValueError
        If no valid `aoi` is provided or the API request fails.
    """
    # convex hull for search
    # convex_hull works better than simplify for more-complex geometry (ie. Louisiana)
    # https://raw.githubusercontent.com/unitedstates/districts/refs/heads/gh-pages/states/LA/shape.geojson
    search_poly = aoi.to_crs(4326).union_all()
    search_poly_chull = search_poly.convex_hull
    coords = ",".join([f"{x},{y}" for x, y in search_poly_chull.exterior.coords])

    # alter the API URL based on the dataset
    if dataset == "noaa":
        url = f"https://portal.opentopography.org/API/otCatalog?productFormat=PointCloud&polygon={coords}&detail=true&outputFormat=json&include_federated=true"
    elif dataset == "ncalm":
        url = f"https://portal.opentopography.org/API/otCatalog?productFormat=PointCloud&polygon={coords}&detail=true&outputFormat=json&include_federated=false"
    else:
        msg = f"Unsupported dataset: {dataset}"
        raise ValueError(msg)

    response = requests.get(url)
    if response.status_code != 200:
        msg = f"Error querying OpenTopography API: {response.status_code}"
        raise ValueError(msg)

    catalog = response.json()
    if dataset == "noaa":
        catalog = [
            ds
            for ds in catalog["Datasets"]
            if "NOAA" in ds["Dataset"]["identifier"]["propertyID"]
        ]
    else:
        catalog = catalog["Datasets"]

    # temporal filtering
    start_dt, end_dt = datetime.fromisoformat(start), datetime.fromisoformat(end)
    filtered_catalog = [
        entry
        for entry in catalog
        if any(
            start_dt <= datetime.fromisoformat(date) <= end_dt
            for date in (entry["Dataset"]["temporalCoverage"].split(" / "))
        )
    ]

    # convert the filtered catalog into a GeoDataFrame
    gf_lidar = gpd.GeoDataFrame(
        [
            {
                "id": entry["Dataset"]["identifier"]["value"],
                "title": entry["Dataset"]["name"],
                "start_datetime": datetime.fromisoformat(
                    entry["Dataset"]["temporalCoverage"].split(" / ")[0]
                ),
                "end_datetime": datetime.fromisoformat(
                    entry["Dataset"]["temporalCoverage"].split(" / ")[-1]
                ),
                "geometry": shape(
                    entry["Dataset"]["spatialCoverage"]["geo"]["geojson"]["features"][
                        0
                    ]["geometry"]
                ),
            }
            for entry in filtered_catalog
        ],
        crs="EPSG:4326",
    )

    # filter the GeoDataFrame to only include datasets that intersect with the search polygon
    return gf_lidar[gf_lidar.intersects(search_poly)].reset_index(drop=True)


# TODO: add meaningful error messages when search fails
def search_opentopo(
    dataset: str | None,
    intersects: gpd.GeoDataFrame | gpd.GeoSeries,
    datetime: list[str] | str,
) -> gpd.GeoDataFrame:
    """
    Integrate the OpenTopography dataset search into the main search function.

    Parameters
    ----------
    dataset : str
        Dataset alias ('noaa' or 'ncalm').
    intersects : gpd.GeoDataFrame | gpd.GeoSeries
        The geometry to restrict the search.
    datetime : list[str] | str
        The datetime range for the search in ISO 8601 format.
    **kwargs : Any
        Additional parameters for the search.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the search results from the OpenTopography API.
    """
    # extract start and end dates
    start, end = datetime if isinstance(datetime, list) else [datetime, datetime]

    # search using the OpenTopography API
    return search_ncalm_noaa(intersects, start=start, end=end, dataset=dataset)
