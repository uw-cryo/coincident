# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository.

## Project Overview

`coincident` is a Python package for searching and analyzing STV (Space
Technology Vehicle) Precursor Coincident Datasets across multiple geospatial
data archives (lidar, stereo imagery, altimetry, SAR). It's developed at the UW
Cryosphere Lab with NASA funding.

## Commands

This project uses `pixi` for environment management. All commands below run
inside the `dev` environment.

```bash
# Run tests (excludes network tests)
pixi run test

# Run all tests including network/API calls
pixi run networktest

# Run a single test
pixi run pytest tests/test_search.py::test_function_name

# Lint
pixi run lint

# Pre-commit hooks
pixi run precommit
```

Tests are marked with `network` for any test requiring internet access. The
default `test` task excludes these via `-m 'not network'`.

Required secrets for network tests: `VANTOR_API_KEY`, `EARTHDATA_TOKEN` (set as
env vars).

## Architecture

Source code lives in `src/coincident/` with five main modules:

### `datasets/`

Dataset metadata definitions. Each dataset (e.g., `ICESat2`, `ThreeDEP`,
`Vantor`) is a class with metadata: alias, STAC collections, search endpoint,
temporal bounds, spatial reference, data type
(`lidar`/`stereo`/`altimeter`/`sar`), and provider info. Files map to providers:
`nasa.py`, `usgs.py`, `vantor.py`, `planetary_computer.py`, `csda.py`,
`opentopo.py`, `neon.py`.

### `search/`

Search APIs that return GeoDataFrames with consistent STAC-compatible columns.
The main entry points are `search()` and `cascading_search()` in `main.py`.
Provider-specific modules handle different APIs: `stac.py` (pystac-client),
`wesm.py` (USGS TNM), `neon_api.py`, `opentopo_api.py`, `vantor.py`. Results are
normalized via `to_geopandas()`.

### `io/`

Reading and downloading data. Organized by return type/library:

- `download.py` — download STAC items via stac-asset (handles Vantor API keys,
  concurrent downloads)
- `xarray.py` — return rasters as xarray DataArrays
- `gdal.py` — GDAL-based reading
- `sliderule.py` — server-side lidar processing via the SlideRule service
- `proj.py` — CRS/projection utilities

### `overlaps/`

Spatial/temporal filtering functions that operate on GeoDataFrames:
`subset_by_temporal_overlap()`, `subset_by_minimum_area()`,
`subset_by_maximum_duration()`, `geographic_area()`.

### `plot/`

Visualization utilities. `matplotlib.py` has plotting functions for DEMs,
elevation distributions, and terrain comparisons. `utils.py` has helpers: aspect
calculation, haversine distance, DEM sampling, map tile fetching.

## Key Design Patterns

- All search functions return `geopandas.GeoDataFrame` with consistent STAC
  columns, enabling uniform downstream processing regardless of data source.
- Dataset classes are data descriptors, not service clients — they hold metadata
  about where/how to search, while `search/` modules handle the actual API
  calls.
- The `rustac` library (Rust-based Arrow STAC parsing) is used for
  performance-critical STAC item parsing.
- Type hints are enforced strictly via mypy for all `coincident.*` modules
  (`strict = true` in pyproject.toml).

## Tooling

- **Packaging:** hatchling (PEP 517), version from `src/coincident/__init__.py`
- **Linting:** ruff (isort + pylint + flake8 plugins), pylint, mypy
- **CI:** GitHub Actions — pre-commit gate, then matrix test (Python
  3.12/3.13/3.14 × ubuntu/macos) using `uv`
- **Docs:** Sphinx + MyST-NB (Jupyter notebooks as examples), hosted on
  ReadTheDocs

## Development

ALWAYS run `pixi run test` after making changes to code. RUN
`pixi run networktest` if you modified any code that interacts with external
APIs or requires internet access.
