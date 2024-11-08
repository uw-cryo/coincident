# Supported datasets

Below we provide a short table summarizing datasets that are searchable with
`coincident`. Note that many of these datasets (or subsets) are available from
different providers, the _provider_ column identifies the source of the data
used by this library.

| Dataset        | Alias      | Type      | Start      | End        | Extent        | Source                                                                      |
| -------------- | ---------- | --------- | ---------- | ---------- | ------------- | --------------------------------------------------------------------------- |
| TanDEM-X       | tdx        | SAR       | 2007-07-01 |            | global        | [NASA CSDAP](https://csdap.earthdata.nasa.gov/stac/collections/airbus)      |
| Maxar Stereo   | maxar      | VHR       | 2007-07-01 |            | global        | [Maxar](https://developers.maxar.com/docs/discovery/)                       |
| Coperincus DEM | cop30      | SAR       | 2021-04-22 |            | global        | [Microsoft](https://planetarycomputer.microsoft.com/dataset/cop-dem-glo-30) |
| ICESat-2 ATL06 | atl06      | Altimeter | 2018-10-13 |            | global        | [NASA](https://nsidc.org/data/atl03)                                        |
| GEDI L2A       | gedi       | Altimeter | 2019-04-04 | 2023-03-17 | mid-latitudes | [NASA](https://lpdaac.usgs.gov/products/gedi02_av002/)                      |
| 3DEP LiDAR     | 3dep       | LiDAR     | 2000-12-01 |            | CONUS         | [USGS](https://www.usgs.gov/3d-elevation-program)                           |
| ESA WorldCover | worldcover | LULC      | 2020-01-01 | 2021-12-31 | global        | [Microsoft](https://planetarycomputer.microsoft.com/dataset/esa-worldcover) |

## Other data sources

If you are interested in working with additional data, feel free to open an
[issue](https://github.com/uw-cryo/coincident/issues).
