# Introduction

`coincident` simplifies access to a curated set of datasets of relevance to NASA
STV studies. It is designed to simplify working with disparate metadata for
areal and satellite remote sensing datasets. `coincident` relies heavily on
[GeoPandas](https://geopandas.org/en/stable/index.html) in that metadata records
are always returned as GeoDataFrame objects, and most methods are written to
operate either on entire dataframes or single rows within a dataframe.

## Dataset aliases

```python
import coincident

coincident.datasets.aliases
```

## Unified search function

the `coincident` package provides a [search()](#coincident.search.search) method
that has the same syntax regardless of which dataset you are searching. Behind
the scenes, polygons intersecting your area of interest are efficiently located
and returned as a geodataframe.

```python
aoi = gpd.read_file(
    "https://raw.githubusercontent.com/unitedstates/districts/refs/heads/gh-pages/states/CO/shape.geojson"
)
gf = coincident.search(
    dataset="3dep",
    intersects=aoi,
    datetime=["2018", "2024"],
)
gf.explore(column="workunit", popup=True)
```

## Convenience functions

`coincident` also provides a number of convenience functions, some of which only
pertain to specific datasets. For example, loading raster imagery via
[Xarray](https://docs.xarray.dev/en/stable) or creating visualizations of browse
imagery. Refer to [the API Docs](../api) for a listing of functions.
