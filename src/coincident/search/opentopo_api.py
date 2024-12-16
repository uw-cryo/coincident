"""
Searching for NCALM and NOAA data with an input search polygon is made easy via the OpenTopography API
https://portal.opentopography.org/apidocs/
There is also a GitHub repo containing flight extent polygons for NOAA
https://github.com/OpenTopography/Data_Catalog_Spatial_Boundaries/tree/main/NOAA_Coastal_Lidar
"""

from __future__ import annotations

from datetime import datetime

import geopandas as gpd
import pandas as pd
import requests
from pandas import Timestamp
from shapely.geometry import shape


def search_ncalm_noaa(
    aoi: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: Timestamp,
    search_end: Timestamp,
    dataset: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Perform a search for geospatial data using OpenTopography API.
    This function dynamically adjusts the API URL based on the dataset.

    Parameters
    ----------
    aoi : gpd.GeoDataFrame | gpd.GeoSeries
        A GeoDataFrame or GeoSeries containing a geometry to restrict the search area.
    search_start : Timestamp, optional
        The start datetime for the search, by default None.
    search_end : Timestamp, optional
        The end datetime for the search, by default None.
    dataset : str
        The dataset type (either "noaa" or "ncalm").

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the search results.

    Raises
    ------
    ValueError
        If no valid `aoi` is provided or the API request fails.
    """
    # If `aoi` is None, set the entire world as the search area
    if aoi is None:
        search_poly = gpd.GeoSeries(
            [
                shape(
                    {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [-180, -90],
                                [180, -90],
                                [180, 90],
                                [-180, 90],
                                [-180, -90],
                            ]
                        ],
                    }
                )
            ]
        )
    else:
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
            and ds["Dataset"]["temporalCoverage"] != ""
        ]
    else:
        catalog = [
            ds for ds in catalog["Datasets"] if ds["Dataset"]["temporalCoverage"] != ""
        ]

    # temporal filtering
    start_dt, end_dt = search_start.to_pydatetime(), search_end.to_pydatetime()
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
    search_start: Timestamp | None = None,
    search_end: Timestamp | None = None,
) -> gpd.GeoDataFrame:
    """
    Integrate the OpenTopography dataset search into the main search function.

    Parameters
    ----------
    dataset : str
        Dataset alias ('noaa' or 'ncalm').
    intersects : gpd.GeoDataFrame | gpd.GeoSeries
        The geometry to restrict the search.
    search_start : Timestamp, optional
        The start time for the temporal search, by default None.
    search_end : Timestamp, optional
        The end time for the temporal search, by default None.
    **kwargs : Any
        Additional parameters for the search.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the search results from the OpenTopography API.
    """
    if search_start is None and search_end is None:
        search_start = pd.Timestamp(
            "1996-01-01"
        )  # Starting from January 1, 1996 (first NOAAA dataset)
        search_end = pd.Timestamp.today()  # Default to today's date

    # search using the OpenTopography API
    return search_ncalm_noaa(
        intersects, search_start=search_start, search_end=search_end, dataset=dataset
    )