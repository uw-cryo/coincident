"""
PROJ functions for accurate 3D Reprojections
"""

from __future__ import annotations

import pandas as pd
import pyproj
from pyproj.crs import CompoundCRS, ProjectedCRS, VerticalCRS
from pyproj.crs.coordinate_operation import UTMConversion
from pyproj.transformer import TransformerGroup


# mypy: disable-error-code="no-any-return"
def _timestamp_to_decimal_year(ts: pd.Timestamp) -> float:
    """Convert pandas Timestamp to decimal year, needed for PROJ epochs"""
    year_start = pd.Timestamp(year=ts.year, month=1, day=1, tz=ts.tz)
    next_year_start = pd.Timestamp(year=ts.year + 1, month=1, day=1, tz=ts.tz)

    year_length = (next_year_start - year_start).total_seconds()
    elapsed = (ts - year_start).total_seconds()

    return ts.year + elapsed / year_length


def construct_custom_utm_crs(
    epsg: str, a_srs: str | None = None, geoid: str | None = None
) -> pyproj.CRS:
    """Custom 3D CRS for UTM + NAVD88 height for NAD83(2011), northern hemisphere

    Parameters
    ----------
    epsg : str
        EPSG code for the horizontal CRS (default 'EPSG:26913' NAD83
    geoid : str | None
        Geoid model (e.g. EGM96, EGM08, GEOID09, GEOID12A, GEOID12B, GEOID18)

    Examples
    -----
    usgs 3dep 1m: NAD83(2011) / UTM zone 13N + NAVD88 height with GEOID18
        construct_custom_utm_crs('EPSG:26913', a_srs='EPSG:6318', geoid='GEOID18')

    copdem: WGS84 (G1150) / UTM zone 13N + EGM2008 height
    """
    crs_default = pyproj.CRS(epsg)
    if a_srs is not None:  # noqa: SIM108
        crs_ellipsoid = pyproj.CRS(a_srs)
    else:
        crs_ellipsoid = pyproj.CRS("EPSG:6318")  # specific NAD83(2011) ellipsoid

    proj_crs = ProjectedCRS(
        name=f"{crs_ellipsoid.name} / UTM zone {crs_default.utm_zone}",
        geodetic_crs=crs_ellipsoid,
        conversion=UTMConversion(crs_default.utm_zone.strip("N")),
    )

    # (EPSG:5703)
    if geoid:
        if geoid.startswith("GEOID"):
            vert_crs = VerticalCRS(
                name="NAVD88 height",
                datum="North American Vertical Datum 1988",
                geoid_model=geoid,
            )
        elif geoid == "EGM96":
            vert_crs = VerticalCRS.from_epsg(5773)
        elif geoid == "EGM08":
            vert_crs = VerticalCRS.from_epsg(3855)

        tmp_crs = CompoundCRS(
            name=f"{crs_ellipsoid.name} / UTM zone {crs_default.utm_zone} + {vert_crs.name}",
            components=[proj_crs, vert_crs],
        )
    else:
        tmp_crs = proj_crs

    # Match area of use from default UTM Zone
    crs_dict = crs_default.to_json_dict()
    tmp = tmp_crs.to_json_dict()
    tmp["area"] = crs_dict["area"]
    tmp["bbox"] = crs_dict["bbox"]

    return pyproj.CRS.from_json_dict(tmp)


# Need AOI or dataset
def get_proj_transform(srcSRS: str, dstSRS: str) -> pyproj.Transformer:
    """
    Get the PROJ pipeline string for a given dataset and destination SRS.

    Parameters
    ----------
    srcSRS : str
        Source spatial reference system. e.g. 'EPSG:7912' or custom pyproj.CRS
    dstSRS : str
        Destination spatial reference system. e.g. 'EPSG:7912' or custom pyproj.CRS
    """
    # Get 3D CRS of dataset from metadata
    # crs_from = coincident.datasets._alias_to_Dataset[dataset].spatial_ref

    # We use TransformerGroup, b/c Transformer sometimes returns None if multiple transformations are possible
    trans_group = TransformerGroup(
        crs_from=srcSRS, crs_to=dstSRS, allow_ballpark=False, always_xy=True
    )

    # based on transform accuracy and PROJ heuristics
    # print(trans_group.transformers[0].to_proj4(pretty=True))
    return trans_group.transformers[0]
