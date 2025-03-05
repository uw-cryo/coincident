##########################################################
#
# inspired by sliderule at
# https://github.com/SlideRuleEarth/sliderule/blob/396fb6f85ea46b97cad43e367107811ae9664f33/datasets/usgs3dep/utils/build_3dep_DEM_geojson.py
# This is just for testing purposes at branch als_dem_fetch that will eventually be deleted
#
##########################################################
from __future__ import annotations

import geopandas as gpd
import requests  # type: ignore[import-untyped]
import rioxarray as rxr
import xarray as xr


# TODO: look to omit this given the main search.search functionality
def aoi_to_polygon_string(aoi: gpd.GeoDataFrame) -> str:
    """
    Helper function for load_usgs_dem()

    Convert an AOI from a GeoDataFrame to a polygon string
    format required by the TNM API.

    Parameters:
      aoi (gpd.GeoDataFrame): GeoDataFrame containing one or more polygon geometries.

    Returns:
      str: A string in the format "x1 y1, x2 y2, ..." representing the polygon coordinates.
    """
    geometry_union = aoi.unary_union
    if geometry_union.geom_type == "Polygon":
        # For a single polygon, use its exterior coordinates
        coords = list(geometry_union.exterior.coords)
    elif geometry_union.geom_type == "MultiPolygon":
        # For multiple polygons, use the convex hull's exterior
        coords = list(geometry_union.convex_hull.exterior.coords)
    else:
        msg_geom_error = "Unsupported geometry type for AOI"
        raise ValueError(msg_geom_error)

    # Format each coordinate pair into the required string format for the TNM API
    return ", ".join([f"{x} {y}" for x, y in coords])


def query_tnm_api(
    polygon_str: str, tnmdataset: str = "Digital Elevation Model (DEM) 1 meter"
) -> list[dict]:
    """
    Helper function for load_usgs_dem()

    Query the TNM API for USGS DEM products that intersect the given polygon.

    Parameters:
      polygon_str (str): Polygon coordinates string in the required API format.
      tnmdataset (str): The dataset name for the TNM API (default is "Digital Elevation Model (DEM) 1 meter").

    Returns:
      list: A list of JSON items returned from the TNM API.
    """
    items = []
    offset = 0
    url = "https://tnmaccess.nationalmap.gov/api/v1/products"

    # paginate through all available products
    # TODO: probably don't want to use a "While True:" loop...
    while True:
        params = {
            "datasets": tnmdataset,
            "polygon": polygon_str,
            "prodFormats": "GeoTIFF",
            "outputFormat": "JSON",
            "max": 100,
            "offset": offset,
        }
        # print(f"Querying TNM API with offset={offset}...")
        response = requests.get(url, params=params)
        if response.status_code != 200:
            msg_tnm_api = (
                f"TNM API request failed with status code {response.status_code}"
            )
            raise RuntimeError(msg_tnm_api)

        data = response.json()
        # Extend the items list with the current batch of results
        items.extend(data.get("items", []))
        total_items = data.get("total", 0)
        offset += len(data.get("items", []))

        # Break the loop if all items have been retrieved
        if offset >= total_items:
            break
    return items


# TODO: fix workunit vs project identifier issue
def filter_items_by_project(api_items: list, project_identifier: str) -> list[dict]:
    """
    Helper function for load_usgs_dem()

    Filter TNM API items to only include those corresponding to a specific project in case
    there are multiple USGS site tiles for a given AOI.

    Parameters:
      api_items (list): List of items (dict) from the TNM API.
      project_identifier (str): The project identifier to filter by (e.g., 'WI_12County_B22').

    Returns:
      list: A filtered list of API items matching the specified project.

    Note: I would ideally like to use the 'workunit' identifier for this but there are some
    inconsistencies in the tile naming conventions where the workunit does not always match the
    expected tile (e.g. CO_CentralEasternPlains_1_2020 maps to tiles for CO_CentralEasternPlains_3_2020)
    """
    filtered_items = []
    for item in api_items:
        vendor_meta_url = item.get("vendorMetaUrl", "")
        if vendor_meta_url:
            try:
                # Extract the project identifier from the vendorMetaUrl.
                remainder = vendor_meta_url.split("metadata/")[-1]
                parts = remainder.split("/")
                if len(parts) >= 2:
                    candidate_project = parts[0]
                    if candidate_project == project_identifier:
                        filtered_items.append(item)
            except Exception as e:
                msg_vendor_url = f"Warning: Unable to parse vendorMetaUrl '{vendor_meta_url}'. Error: {e}"
                print(msg_vendor_url)  # noqa: T201
    return filtered_items


def load_usgs_dem(
    aoi: gpd.GeoDataFrame,
    project: str,
    tnmdataset: str = "Digital Elevation Model (DEM) 1 meter",
    res: int = 1,
    clip: bool = True,
) -> xr.DataArray:
    """
    Load and merge USGS 1-meter DEM tiles based on an AOI by querying the TNM API.

    Steps:
      1. Reproject the AOI to EPSG:4326 (required by the TNM API).
      2. Convert the AOI into a polygon string for the API query.
      3. Query the TNM API for matching DEM products.
      4. Filter the API items by the given project identifier.
      5. Open each GeoTIFF as a lazy xarray.DataArray using rioxarray.
      6. Optionally coarsen the DEM data to a lower resolution by aggregating over
         res x res blocks (default res=1 means no coarsening).
      7. Merge the DEM tiles by taking the maximum value along the new 'band' dimension
         (preferring non-NaN values).
      8. Optionally clip the merged DEM mosaic to the original AOI geometry.

    Parameters:
      aoi (gpd.GeoDataFrame): The input Area of Interest as a GeoDataFrame.
      project (str): A project identifier string that is part of the vendor metadata.
      tnmdataset (str): The dataset identifier for the TNM API query (default "Digital Elevation Model (DEM) 1 meter").
      res (int): The coarsening factor. A value > 1 aggregates 1x1 m pixels into res x res m pixels (default 1).
      clip (bool): Whether to clip the final mosaic to the AOI geometry (default True).

    Returns:
      xr.DataArray: An xarray.DataArray containing the merged (and optionally clipped) DEM mosaic,
                    with a variable name "elevation".
    """
    # 1: Reproject AOI to EPSG:4326 for the TNM API query.
    aoi_in_4326: gpd.GeoDataFrame = aoi.to_crs(epsg=4326)

    # 2: Convert AOI geometry to the polygon string format required by the API.
    aoi_polygon_string: str = aoi_to_polygon_string(aoi_in_4326)
    # msg_poly = "Polygon string for TNM API:", aoi_polygon_string
    # print(msg_poly)

    # 3: Query the TNM API for DEM products intersecting the AOI.
    tnm_api_items: list[dict] = query_tnm_api(aoi_polygon_string, tnmdataset=tnmdataset)
    if not tnm_api_items:
        msg_no_hits = "TNM API returned no products for the given AOI and dataset."
        raise ValueError(msg_no_hits)

    # 4: Filter the returned items to include only those matching the specified project.
    filtered_api_items: list[dict] = filter_items_by_project(tnm_api_items, project)
    if not filtered_api_items:
        msg_no_products = (
            f"No products found for project '{project}' intersecting the AOI."
        )
        raise ValueError(msg_no_products)
    # msg_products = f"Found {len(filtered_api_items)} products for project '{project}'."
    # print(msg_products)

    # 5: Extract TIFF URLs from the filtered API items.
    tiff_urls: list = [
        item["urls"]["TIFF"]
        for item in filtered_api_items
        if "urls" in item and "TIFF" in item["urls"]
    ]

    # 6: Open each GeoTIFF file lazily using rioxarray.
    dem_dataarrays: list = []
    for url in tiff_urls:
        try:
            dem_tile: xr.DataArray = rxr.open_rasterio(url, masked=True, chunks="auto")
            # If a coarser resolution is requested, coarsen
            if res > 1:
                dem_tile = dem_tile.coarsen(x=res, y=res, boundary="trim").mean()
            dem_dataarrays.append(dem_tile)
        except Exception as error:
            msg_url_open = f"Warning: Failed to open {url}. Error: {error}"
            print(msg_url_open)  # noqa: T201

    if not dem_dataarrays:
        msg_runtime = "None of the TIFF files could be opened."
        raise RuntimeError(msg_runtime)

    # 7: Assume all DEM tiles share the same CRS (usually a UTM). Use the CRS from the first tile.
    tile_crs = dem_dataarrays[0].rio.crs

    # Merge the DEM tiles along a new "band" dimension and select non-NaN values using max().
    # This creates a single DEM mosaic.
    merged_dem: xr.DataArray = (
        xr.concat(dem_dataarrays, dim="band").max(dim="band").rename("elevation")
    )

    # 8: If clipping is enabled, clip the mosaic to the AOI geometry.
    if clip:
        # Reproject AOI to the CRS of the DEM tiles if needed.
        aoi_tile_crs: gpd.GeoDataFrame = aoi.to_crs(tile_crs)
        merged_dem = merged_dem.rio.clip(aoi_tile_crs.geometry, tile_crs, drop=True)

    return merged_dem
