# Input/Output

This section covers suggestions on writing outputs from using `coincident`, And
how to use those outputs with other software libraries. For saving data we
recommend using
[Cloud-Optimized Geospatial Formats](https://guide.cloudnativegeo.org).

## Saving search results

Coincident results are returned as a GeoDataFrame. For small numbers of results
(<1000) saving to a GeoJSON file is a convenient option
(`gf.to_file('results.geojson')`).

However, for larger searches saving to
[GeoParquet](https://geopandas.org/en/stable/docs/user_guide/io.html#apache-parquet-and-feather-file-formats)
is preferable (`gf.to_parquet('results.parquet')`). GeoParquet has the advantage
of saving nested datatypes that are sometimes present in STAC metadata that are
not supported in GeoJSON (for example lists). This allows re-opening results
later (`gf = gpd.read_parquet('results.parquet')`) or converting to STAC JSON
format.

## Saving data

### Point data

Querying raster values or retrieving altimeter point measurements can result in
very large tables. Saving to GeoParquet as described above is preferred.

### Raster data

Coincident returns raster data as Xarray objects. For single 2D images, we
recommend saving Cloud-optimized geotiffs via the
[rioxarray extension](https://corteva.github.io/rioxarray/stable/examples/convert_to_raster.html)
(`da.rio.to_raster('mysubset.tif', driver='COG')`)

For larger multi-dimensional or multi-image datasets saving to
[Zarr](https://docs.xarray.dev/en/latest/user-guide/io.html#zarr) format is
recommended (`ds.to_zarr('mydataset.zarr')`).

## QGIS

To work with GeoParquet in QGIS using the most recent released version is
recommended (GDAL, which QGIS depends on has been steadily adding support and
performance enhancements for GeoParquet). There are different way to install
QGIS, and we recommend the following approach to install the latest release with
pixi:

```bash
pixi global install qgis
pixi global add libgdal-arrow-parquet -e qgis
```

Then from a terminal launch QGIS simply by typing `qgis`
