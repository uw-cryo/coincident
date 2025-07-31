"""
Functions for searching, with a spatiotemporal filter, for:
- NGS survey mark control monuments
- Permanent GNSS control position(s)
- Potentially USGS 3DEP GCPs in the future...
"""

from __future__ import annotations

import logging
import time

import geopandas as gpd
import pandas as pd
import requests  # type: ignore[import-untyped]

log = logging.getLogger(__name__)

# module-level constants for metadata variable mapping
_DATUM_TO_EPSG_MAP = {
    "NAD 83(2011)": "EPSG:6318",
    "NAD 83(NSRS2007)": "EPSG:4759",
    "NAD 83(HARN)": "EPSG:4152",
    "NAD 83(FBN)": "EPSG:4152",
    "NAD 83(1991)": "EPSG:4152",
    "NAD 83(1992)": "EPSG:4134",
    "NAD 83(1993)": "EPSG:4135",
    "NAD 83(1986)": "EPSG:4269",
    "NAD 27": "EPSG:4267",
    "NAD 83(PA11)": "EPSG:6323",
    "NAD 83(MA11)": "EPSG:6322",
    "NAD 83": "EPSG:4269",
    # "NAD 83(CORS)": "EPSG:6318",
}

_CONTROL_TYPE_MAPPING = {
    "h1": "Horizontal Control Adjusted",
    "h2": "Horizontal Control Other",
    "v1": "Vertical Control Adjusted",
    "v2": "Vertical Control Other",
    "v1h1": "Combined Control (v1h1)",
    "v1h2": "Combined Control (v1h2)",
    "v2h1": "Combined Control (v2h1)",
    "v2h2": "Combined Control (v2h2)",
    "": "Not Control",
}


def _parse_ngs_recovery_date(date_str: str | None) -> pd.Timestamp:
    """Parses an NGS 'lastRecovered' field (YYYY or YYYYMMDD) into a datetime object."""
    if not isinstance(date_str, str) or not date_str.strip():
        return pd.NaT
    date_str = date_str.strip()
    try:
        if len(date_str) == 8:
            return pd.to_datetime(date_str, format="%Y%m%d")
        if len(date_str) == 4:
            return pd.to_datetime(f"{date_str}-01-01")
        return pd.NaT
    except (ValueError, TypeError):
        return pd.NaT


def _determine_ngs_control_type(row: pd.Series) -> str:
    """Determines the control type and returns a descriptive string."""
    h_part = ""
    v_part = ""
    pos_source = str(row.get("posSource", "")).strip()
    vert_source = str(row.get("vertSource", "")).strip()

    if pos_source:
        h_part = "h1" if pos_source in ("ADJUSTED", "SCALED") else "h2"
    if vert_source:
        v_part = "v1" if vert_source in ("ADJUSTED", "SCALED") else "v2"

    code = f"{v_part}{h_part}"
    return _CONTROL_TYPE_MAPPING.get(code, "Not Control")


def _query_ngs_stations_by_bounds(
    minlat: float, maxlat: float, minlon: float, maxlon: float
) -> gpd.GeoDataFrame:
    """
    Fetches all NGS survey stations within a specified bounding box.

    This function bypasses the 500-item API limit by recursively subdividing
    the bounding box and groups stations by their source datum for rigorous
    CRS transformation.

    Args:
        minlat: Minimum latitude in decimal degrees.
        maxlat: Maximum latitude in decimal degrees.
        minlon: Minimum longitude in decimal degrees.
        maxlon: Maximum longitude in decimal degrees.

    Returns:
        A GeoDataFrame with station data in ITRF2020 (EPSG:9989), or an empty
        GeoDataFrame if an error occurs.
    """
    seen_pids = set()
    all_stations = []

    log.info("Starting NGS station fetch. This may take time for large areas...")

    def _fetch_ngs_recursive(
        b_minlat: float, b_maxlat: float, b_minlon: float, b_maxlon: float
    ) -> None:
        """
        Nested function to perform the recursive fetch.
        API return limit is 500, so here we recursively split areas over that limit
        into quadrants and search the quadrants
        """
        base_url = "https://geodesy.noaa.gov/api/nde/bounds"
        params = {
            "minlat": b_minlat,
            "maxlat": b_maxlat,
            "minlon": b_minlon,
            "maxlon": b_maxlon,
        }
        try:
            time.sleep(0.1)
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            data = response.json()
            if len(data) == 500:
                msg = f"Limit reached for box ({b_minlat:.4f}, {b_minlon:.4f} to {b_maxlat:.4f}, {b_maxlon:.4f}). Subdividing..."
                log.info(msg)
                mid_lat, mid_lon = (b_minlat + b_maxlat) / 2, (b_minlon + b_maxlon) / 2
                _fetch_ngs_recursive(b_minlat, mid_lat, b_minlon, mid_lon)
                _fetch_ngs_recursive(b_minlat, mid_lat, mid_lon, b_maxlon)
                _fetch_ngs_recursive(mid_lat, b_maxlat, b_minlon, mid_lon)
                _fetch_ngs_recursive(mid_lat, b_maxlat, mid_lon, b_maxlon)
            else:
                for station in data:
                    if (pid := station.get("pid")) and (pid not in seen_pids):
                        seen_pids.add(pid)
                        all_stations.append(station)
        except requests.exceptions.RequestException as e:
            msg = f"API request failed for box ({b_minlat:.4f}, {b_minlon:.4f}): {e}"
            log.error(msg)

    _fetch_ngs_recursive(minlat, maxlat, minlon, maxlon)

    if not all_stations:
        log.warning("No stations found in the specified bounding box.")
        return gpd.GeoDataFrame()

    msg = f"\nTotal unique stations found: {len(all_stations)}"
    log.info(msg)
    df = pd.DataFrame(all_stations)

    df["controlType"] = df.apply(_determine_ngs_control_type, axis=1)
    df["recoveryDate"] = (
        df["lastRecovered"].apply(_parse_ngs_recovery_date).dt.tz_localize(None)
    )
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lon"] = pd.to_numeric(df["lon"], errors="coerce")
    df = df.dropna(subset=["lat", "lon", "posDatum"])
    df = df[df["posDatum"].str.strip() != ""]

    gdf_list = []
    log.info("Grouping stations by datum for rigorous CRS transformation...")
    for datum_name, group in df.groupby("posDatum"):
        if epsg_code := _DATUM_TO_EPSG_MAP.get(datum_name.strip()):
            msg = f"  - Processing {len(group)} stations in {datum_name} ({epsg_code})."
            log.info(msg)
            temp_gdf = gpd.GeoDataFrame(
                group,
                geometry=gpd.points_from_xy(group.lon, group.lat),
                crs=epsg_code,
            )
            gdf_list.append(temp_gdf.to_crs("EPSG:9989"))  # Transform to ITRF2020
        else:
            msg = f"Warning: No EPSG code for datum '{datum_name}'. Skipping {len(group)} stations."
            log.warning(msg)

    if not gdf_list:
        log.warning("No stations could be processed due to missing datum information.")
        return gpd.GeoDataFrame()

    log.info("Concatenating all processed groups into a final GeoDataFrame...")
    gdf = pd.concat(gdf_list, ignore_index=True)
    log.info("Transformation complete.")

    final_cols = [
        "pid",
        "corsId",
        "name",
        "recoveryDate",
        "condition",
        "controlType",
        "posDatum",
        "vertDatum",
        "monumentType",
        "geometry",
    ]
    existing_cols = [col for col in final_cols if col in gdf.columns]
    return gdf[existing_cols]


def search_ngs_stations(
    intersects: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: pd.Timestamp | None = None,
) -> gpd.GeoDataFrame:
    """
    Searches for NGS stations and applies precise spatial and temporal filters.

    This function is designed to be called from a our main search. It uses
    the total bounds of the input geometry to query the NGS API, then performs
    a spatial intersection and a temporal filter on the results.
    Point geometries are defined in their NGS-reported horizontal CRS and then
    transformed to ITRF 2020.

    Args:
        intersects: A GeoDataFrame or GeoSeries containing the search geometry.
        search_start: The start of the temporal search window. Results will have a
                      recovery date on or after this time.

    Returns:
        A GeoDataFrame containing the filtered NGS station results.
    """
    # Get the simple bounding box from the complex geometry for the API call
    bounds = intersects.total_bounds
    minlon, minlat, maxlon, maxlat = bounds

    # Fetch all stations within the simple bounding box
    gdf = _query_ngs_stations_by_bounds(minlat, maxlat, minlon, maxlon)

    if gdf.empty:
        return gdf

    # filter using the actual geometry
    search_geom = intersects.union_all()
    gdf = gdf[gdf.intersects(search_geom)]

    # temporal filter, here we just use the recoveryDate  as being later than the
    # minimum search datetime so it's difficult to get an end datetime for NGS stations
    # with the API we're using
    if search_start is not None:
        # Ensure search_start is timezone-naive for comparison
        search_start_naive = search_start.tz_localize(None)
        gdf = gdf[gdf["recoveryDate"] >= search_start_naive]

    return gdf.reset_index(drop=True).rename(
        columns={"recoveryDate": "latest_recoveryDate"}
    )


def get_ngs_opus_solutions(input_gdf: gpd.GeoDataFrame, date: str) -> gpd.GeoDataFrame:
    """
    Finds the OPUS solution for each station in a GeoDataFrame closest to a target date.

    NOTE: Most stations returned from coincident.search.search(dataset='ngs-stations') will
    NOT have OPUS solutions

    Args:
        input_gdf: GeoDataFrame containing NGS stations with a 'pid' column.
        date: The desired date for the solution in "YYYY-MM-DD" format.

    Returns:
        A DataFrame with detailed OPUS solutions for each station.
    """
    base_url = "https://geodesy.noaa.gov/api/opus/pid"
    target_date = pd.to_datetime(date).tz_localize(None)
    results = []

    msg = f"Fetching OPUS solutions closest to {date}..."
    log.info(msg)
    for station in input_gdf.itertuples():
        pid = station.pid
        params = {"pid": pid}
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            solutions = response.json()

            if not isinstance(solutions, list) or not solutions:
                msg = f"  - No OPUS solutions found for PID: {pid}"
                log.warning(msg)
                continue

            best_solution = min(
                solutions,
                key=lambda s: abs(
                    pd.to_datetime(s["observed"]).tz_localize(None) - target_date
                ),
            )
            solution_date = pd.to_datetime(best_solution["observed"]).tz_localize(None)

            if solution_date.date() != target_date.date():
                msg = f"  - No solution on {date} for PID {pid}. Using closest: {solution_date.strftime('%Y-%m-%d')}."
                log.warning(msg)

            results.append(
                {
                    "pid": pid,
                    "corsId": station.corsId if "corsId" in station._fields else "",
                    "name": best_solution.get("designation", station.name),
                    "solutionDate": solution_date,
                    "horizontalDatum": best_solution.get("refFrame"),
                    "verticalDatum": best_solution.get("datum"),
                    "ellipsoidHeight_m": pd.to_numeric(best_solution.get("ellHt")),
                    "orthometricHeight_m": pd.to_numeric(best_solution.get("orthoHt")),
                    "accuracy_lat_p2p_m": pd.to_numeric(best_solution.get("latP2p")),
                    "accuracy_lon_p2p_m": pd.to_numeric(best_solution.get("lonP2p")),
                    "accuracy_ellHt_p2p_m": pd.to_numeric(
                        best_solution.get("ellHtP2p")
                    ),
                    "accuracy_orthoHt_p2p_m": pd.to_numeric(
                        best_solution.get("orthoHtP2p")
                    ),
                    "lat": pd.to_numeric(best_solution.get("lat")),
                    "lon": pd.to_numeric(best_solution.get("lon")),
                }
            )
        except requests.exceptions.RequestException as e:
            msg = f"  - API request failed for PID {pid}: {e}"
            log.error(msg)
        except (ValueError, KeyError) as e:
            msg = f"  - Failed to parse data for PID {pid}: {e}"
            log.warning(msg)

    if not results:
        log.warning("No valid OPUS solutions could be compiled.")
        return gpd.GeoDataFrame()

    df = pd.DataFrame(results)

    gdf_list = []
    for datum_name, group in df.groupby("horizontalDatum"):
        # Normalize datum name from OPUS (e.g., 'NAD_83') to match our map ('NAD 83')
        normalized_datum = datum_name.replace("_", " ")
        if epsg_code := _DATUM_TO_EPSG_MAP.get(normalized_datum):
            temp_gdf = gpd.GeoDataFrame(
                group,
                geometry=gpd.points_from_xy(group.lon, group.lat),
                crs=epsg_code,
            )
            gdf_list.append(temp_gdf.to_crs("EPSG:9989"))  # Transform to ITRF2020
        else:
            msg = f"Warning: No EPSG code for OPUS datum '{datum_name}'. Assuming NAD 83(2011)."
            log.warning(msg)
            temp_gdf = gpd.GeoDataFrame(
                group,
                geometry=gpd.points_from_xy(group.lon, group.lat),
                crs="EPSG:6318",
            )
            gdf_list.append(temp_gdf.to_crs("EPSG:9989"))

    log.info("Processing complete.")
    return pd.concat(gdf_list, ignore_index=True).drop(columns=["lat", "lon"])


def search_gage_stations(
    intersects: gpd.GeoDataFrame | gpd.GeoSeries,
    search_start: pd.Timestamp | None = None,
    search_end: pd.Timestamp | None = None,
) -> gpd.GeoDataFrame:
    """
    Searches for GAGE GNSS/GPS stations using spatial and temporal filters.
    https://www.unavco.org/data/web-services/documentation/documentation.html#!/GNSS47GPS/getGpsSiteMetadata

    This function queries the GAGE (formerly UNAVCO) web services for stations
    within the bounding box of the input geometry. It then performs precise
    spatial and temporal filtering on the results.

    The final coordinates are transformed to ITRF2020 (EPSG:9989) from WGS84

    Args:
        intersects: A GeoDataFrame or GeoSeries defining the spatial search area.
        search_start: The start of the temporal search window. Stations active
                      on or after this time will be included.
        search_end: The end of the temporal search window. Stations active
                    on or before this time will be included.

    Returns:
        A GeoDataFrame with filtered GAGE station data in ITRF2020 (EPSG:9989),
        or an empty GeoDataFrame if no matching stations are found.
    """
    if intersects.empty:
        log.warning("Input geometry is empty. Returning empty GeoDataFrame.")
        return gpd.GeoDataFrame()

    # simple bounding box for the initial API query
    minlon, minlat, maxlon, maxlat = intersects.total_bounds
    url = (
        "https://web-services.unavco.org/gps/metadata/sites/v1"
        f"?minlatitude={minlat}&maxlatitude={maxlat}"
        f"&minlongitude={minlon}&maxlongitude={maxlon}"
    )

    # fetch data from the GAGE API
    try:
        response = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
        response.raise_for_status()  # Raise an exception for bad status codes
        data = response.json()
    except (
        requests.exceptions.RequestException,
        requests.exceptions.JSONDecodeError,
    ) as e:
        msg_api_error = f"Failed to fetch or parse data from GAGE API: {e}"
        log.error(msg_api_error)
        return gpd.GeoDataFrame()

    if not data:
        msg_empty_search = (
            "No stations found in the specified bounding box from GAGE API."
        )
        log.info(msg_empty_search)
        return gpd.GeoDataFrame()

    df = pd.DataFrame(data)

    # Convert date columns, coercing errors to NaT
    date_cols = ["session_start_time", "session_stop_time", "latest_data_time"]
    for col in date_cols:
        df[col] = pd.to_datetime(df[col], errors="coerce")

    # drop rows where critical information is missing
    df = df.dropna(subset=["id", "latitude", "longitude", "latest_data_time"])
    if df.empty:
        log.info("No valid stations after initial data cleaning.")
        return gpd.GeoDataFrame()

    # deduplicate
    df = df.sort_values("latest_data_time", ascending=False).drop_duplicates(
        subset="id", keep="first"
    )

    # Temporal Filtering
    if search_start:
        # ensure search_start is timezone-aware for comparison
        search_start_aware = (
            search_start.tz_localize("UTC")
            if search_start.tzinfo is None
            else search_start
        )
        df = df[df["session_stop_time"] >= search_start_aware]

    if search_end:
        # ensure search_end is timezone-aware
        search_end_aware = (
            search_end.tz_localize("UTC") if search_end.tzinfo is None else search_end
        )
        df = df[df["session_start_time"] <= search_end_aware]

    if df.empty:
        log.info("No stations found within the specified time window.")
        return gpd.GeoDataFrame()

    # initial GDF in WGS84 for spatial filtering
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.longitude, df.latitude),
        crs="EPSG:4326",  # GAGE API uses WGS 84
    )
    search_geom = intersects.to_crs(gdf.crs).union_all()
    gdf = gdf[gdf.intersects(search_geom)]

    if gdf.empty:
        log.info("No stations intersect the precise search geometry.")
        return gpd.GeoDataFrame()

    # Transform to ITRF2020
    gdf_itrf = gdf.to_crs("EPSG:9989")

    return gdf_itrf.reset_index(drop=True)
