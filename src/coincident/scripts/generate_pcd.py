from __future__ import annotations

import logging
from pathlib import Path

import geopandas as gpd
import requests  # type: ignore[import-untyped]

OWNER = "uw-cryo"
REPO = "coincident"
# RELEASE_TAG = "v0.4"
RELEASE_TAG = "latest"
# create a dir to store downloaded and generated assets?
OUTPUT_DIR = Path("data/pcd_assets")


def fetch_and_convert_PCD_assets() -> None:
    """
    Fetches, downloads, and converts all PCD .geojson and .parquet assets
    from the specified GitHub release.
    """
    msg_release_fetch = f" Fetching release information for tag '{RELEASE_TAG}'..."
    logging.debug(msg_release_fetch)
    if RELEASE_TAG == "latest":
        api_url = f"https://api.github.com/repos/{OWNER}/{REPO}/releases/latest"
    else:
        api_url = (
            f"https://api.github.com/repos/{OWNER}/{REPO}/releases/tags/{RELEASE_TAG}"
        )

    response = requests.get(api_url)
    response.raise_for_status()  # raise an exception for HTTP errors
    assets = response.json().get("assets", [])

    if not assets:
        msg_no_asset = "No assets found in the release."
        logging.debug(msg_no_asset)
        return

    # create the output directory if it doesn't exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    msg_print_out_dir = f" Output directory is: {OUTPUT_DIR}"
    logging.debug(msg_print_out_dir)

    pcd_assets = [
        asset
        for asset in assets
        if "PCD" in asset["name"]
        and (asset["name"].endswith(".geojson") or asset["name"].endswith(".parquet"))
    ]

    for asset in pcd_assets:
        asset_name = asset["name"]
        download_url = asset["browser_download_url"]

        temp_path = OUTPUT_DIR / asset_name

        msg_download_asset = f"\nDownloading {asset_name}..."
        logging.debug(msg_download_asset)
        with requests.get(download_url, stream=True) as r:
            r.raise_for_status()
            with Path(temp_path).open("wb") as f:
                for chunk in r.iter_content():
                    f.write(chunk)

        msg_asset_err = f"  - Reading {asset_name} into GeoDataFrame..."
        logging.debug(msg_asset_err)
        try:
            if temp_path.suffix == ".parquet":
                gdf = gpd.read_parquet(temp_path)
            else:  # .geojson
                gdf = gpd.read_file(temp_path)
        except Exception as e:
            msg_read_error = f"  -  Could not read file {asset_name}: {e}"
            logging.debug(msg_read_error)
            temp_path.unlink()  # Clean up failed download
            continue

        output_basename = temp_path.stem

        gpkg_path = OUTPUT_DIR / f"{output_basename}.geojson"
        msg_save_gj = f"  -   Saving to {gpkg_path.name}..."
        logging.debug(msg_save_gj)
        gdf.to_file(gpkg_path, driver="GeoJSON")

        geoparquet_path = OUTPUT_DIR / f"{output_basename}.geoparquet"
        msg_save_parquet = f"  -  Saving to {geoparquet_path.name}..."
        logging.debug(msg_save_parquet)
        gdf.to_parquet(geoparquet_path)

        # Clean up the original downloaded file
        temp_path.unlink()

    msg_files_downloaded = "\n All assets processed successfully."
    logging.debug(msg_files_downloaded)


if __name__ == "__main__":
    fetch_and_convert_PCD_assets()
