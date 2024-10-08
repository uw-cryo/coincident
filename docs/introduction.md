# Introduction

`coincident` simplifies access to a curated set of datasets of relevance to NASA
STV studies

```python
import coincident

coincident.datasets.aliases
```

the `coincident` package provides a `search()` method that has the same syntax
regardless of which dataset you are searching. Behind the scenes, polygons
intersecting your area of interest are efficiently located and returned as a
geodataframe.

```python
aoi = gpd.read_file(
    "https://raw.githubusercontent.com/unitedstates/districts/refs/heads/gh-pages/states/CO/shape.geojson"
)
gf = coincident.search.search(
    dataset="3dep",
    intersects=aoi,
    datetime=["2018", "2024"],
)
gf.explore(column="workunit", popup=True)
```
