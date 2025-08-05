"""
Access the metadata on-demand for the Precursor Coincident Dataset (PCD),
which consists of 13 high-quality 3D measurement sites across the contiguous
United States (CONUS).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TypedDict

import geopandas as gpd
import pandas as pd
from shapely import wkt

from coincident.search import search, wesm


# types for mypy:
class _3depFilter(TypedDict):
    workunit: str


class _NeonFilter(TypedDict):
    site: str


class _MaxarFilter(TypedDict):
    stereo_ids: list[str]


class _SatelliteFilter(TypedDict):
    granule_ids: list[str]


_PcdFilters = TypedDict(
    "_PcdFilters",
    {
        "3dep": _3depFilter,
        "neon": _NeonFilter,
        "ncalm": dict[str, str],
        "maxar": _MaxarFilter,
        "icesat-2": _SatelliteFilter,
        "gedi": _SatelliteFilter,
    },
    total=False,  # makes all keys in this TypedDict optional
)


class PcdSiteParams(TypedDict):
    """Defines the top-level structure for each site in PCD_SITES"""

    provider: str
    datetime: list[str]
    overlap_geometry: str
    filters: _PcdFilters


PCD_SITES: dict[str, PcdSiteParams] = {
    "CA_SanFrancisco_1_B23": {
        "provider": "USGS",
        "datetime": ["2023-01-17", "2023-07-18"],
        "overlap_geometry": "POLYGON ((-122.49869579540385 37.7225321435048, -122.51421754728509 37.78265684820127, -122.47378332685508 37.80994995226257, -122.39815883627537 37.80791177118956, -122.4616240341796 37.74745501281027, -122.49869579540385 37.7225321435048))",
        "filters": {
            "3dep": {"workunit": "CA_SanFrancisco_1_B23"},
            "maxar": {"stereo_ids": ["10200100DE097E00", "10200100DC70E900"]},
            "icesat-2": {
                "granule_ids": [
                    "ATL03_20230117092727_04251806_006_02.h5",
                    "ATL03_20230418050714_04251906_006_02.h5",
                    "ATL03_20230522151905_09511902_006_02.h5",
                    "ATL03_20230718004632_04252006_006_02.h5",
                ]
            },
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2023028212552_O23387_02_T09397_02_003_02_V002"
                ]
            },
        },
    },
    "AZ_PimaCo_2_2021": {
        "provider": "USGS",
        "datetime": ["2021-09-18", "2021-11-29"],
        "overlap_geometry": "POLYGON ((-111.0813927224041 31.894609808770994, -111.01953572428549 31.94798227506963, -111.01317668582358 31.947098123594788, -111.01244165184261 31.886079708871957, -111.07866693825639 31.82683230205409, -111.0813927224041 31.894609808770994))",
        "filters": {
            "3dep": {"workunit": "AZ_PimaCo_2_2021"},
            "maxar": {"stereo_ids": ["104001006D788E00", "104001006E942500"]},
            "icesat-2": {
                "granule_ids": [
                    "ATL03_20211018062519_03941306_006_01.h5",
                    "ATL03_20211129162037_10421302_006_01.h5",
                ]
            },
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2021261011319_O15667_02_T08065_02_003_02_V002",
                    "GEDI02_A_2021264234230_O15728_02_T10605_02_003_02_V002",
                    "GEDI02_A_2021285225842_O16053_03_T08636_02_003_02_V002",
                    "GEDI02_A_2021332204942_O16780_02_T08218_02_003_02_V002",
                ]
            },
        },
    },
    "NE_Northeast_Phase2_2_2020": {
        "provider": "USGS",
        "datetime": ["2020-11-02", "2020-12-22"],
        "overlap_geometry": "POLYGON ((-98.283830547 40.8319992894, -98.28805470114001 41.185077194343826, -98.26500850845709 41.26193292321868, -98.12052091014777 41.26082060738351, -98.16070277564374 40.89950149849092, -98.283830547 40.8319992894))",
        "filters": {
            "3dep": {"workunit": "NE_Northeast_Phase2_2_2020"},
            "maxar": {
                "stereo_ids": [
                    "10200100A1228A00",
                    "102001009D48BC00",
                    "10200100A0BEE600",
                    "10200100A0A9E600",
                ]
            },
            "icesat-2": {
                "granule_ids": [
                    "ATL03_20201119213030_08660906_006_01.h5",
                    "ATL03_20201220075046_13310902_006_01.h5",
                ]
            },
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2020307133033_O10715_03_T09706_02_003_02_V002",
                    "GEDI02_A_2020311115609_O10776_03_T11282_02_003_02_V002",
                    "GEDI02_A_2020314045702_O10818_02_T06489_02_003_02_V002",
                    "GEDI02_A_2020315102144_O10837_03_T10012_02_003_02_V002",
                    "GEDI02_A_2020318032241_O10879_02_T06642_02_003_02_V002",
                    "GEDI02_A_2020323071439_O10959_03_T07166_02_003_02_V002",
                    "GEDI02_A_2020326001628_O11001_02_T10911_02_003_02_V002",
                    "GEDI02_A_2020327054135_O11020_03_T09859_02_003_02_V002",
                    "GEDI02_A_2020329224322_O11062_02_T09335_02_003_02_V002",
                    "GEDI02_A_2020331040823_O11081_03_T07013_02_003_02_V002",
                    "GEDI02_A_2020333211005_O11123_02_T07912_02_003_02_V002",
                    "GEDI02_A_2020335023508_O11142_03_T08436_02_003_02_V002",
                    "GEDI02_A_2020337193638_O11184_02_T10758_02_003_02_V002",
                    "GEDI02_A_2020341180306_O11245_02_T09488_02_003_02_V002",
                    "GEDI02_A_2020345162931_O11306_02_T11064_02_003_02_V002",
                    "GEDI02_A_2020346215427_O11325_03_T07166_02_003_02_V002",
                    "GEDI02_A_2020349145551_O11367_02_T08218_02_003_02_V002",
                    "GEDI02_A_2020353132211_O11428_02_T06795_02_003_02_V002",
                    "GEDI02_A_2020357114828_O11489_02_T09641_02_003_02_V002",
                ]
            },
        },
    },
    "WI_Brown_2_2020": {
        "provider": "USGS",
        "datetime": ["2020-04-07", "2020-05-09"],
        "overlap_geometry": "POLYGON ((-88.08851055495393 44.65002443734368, -88.07680590851425 44.67237307709082, -88.03832665179429 44.655758940856884, -88.0173511796351 44.62591728077161, -88.03829737314108 44.49921466897285, -88.0855997216971 44.49804855713017, -88.08851055495393 44.65002443734368))",
        "filters": {
            "3dep": {"workunit": "WI_Brown_2_2020"},
            "maxar": {"stereo_ids": ["104001005ACC6E00", "104001005B461800"]},
            "icesat-2": {"granule_ids": ["ATL03_20200503062653_05760706_006_01.h5"]},
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2020098225152_O07481_03_T04610_02_003_01_V002",
                    "GEDI02_A_2020126115921_O07908_03_T04457_02_003_01_V002",
                    "GEDI02_A_2020130102505_O07969_03_T01764_02_003_01_V002",
                ]
            },
        },
    },
    "GA_Central_3_2019": {
        "provider": "USGS",
        "datetime": ["2020-02-02", "2020-04-01"],
        "overlap_geometry": "POLYGON ((-84.12671569644924 31.01941095632088, -84.13329816636715 31.568656984115982, -84.11305907559789 31.746612946981553, -84.05059072866084 31.69169387974824, -84.0310090441035 31.606427307310568, -84.02359554302821 31.222141221316225, -84.0384748132213 31.02889768223571, -84.12671569644924 31.01941095632088))",
        "filters": {
            "3dep": {"workunit": "GA_Central_3_2019"},
            "maxar": {"stereo_ids": ["10300100A26FD700", "10300100A132F900"]},
            "icesat-2": {
                "granule_ids": [
                    "ATL03_20200203102127_05910606_006_01.h5",
                    "ATL03_20200316201642_12390602_006_01.h5",
                ]
            },
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2020040225749_O06582_03_T02744_02_003_01_V002",
                    "GEDI02_A_2020044212416_O06643_03_T05590_02_003_01_V002",
                    "GEDI02_A_2020048195041_O06704_03_T01321_02_003_01_V002",
                    "GEDI02_A_2020052181703_O06765_03_T04320_02_003_01_V002",
                    "GEDI02_A_2020056164321_O06826_03_T04167_02_003_01_V002",
                    "GEDI02_A_2020060150939_O06887_03_T01474_02_003_01_V002",
                    "GEDI02_A_2020064133554_O06948_03_T04320_02_003_01_V002",
                    "GEDI02_A_2020068120206_O07009_03_T00051_02_003_01_V002",
                    "GEDI02_A_2020072102816_O07070_03_T00204_02_003_01_V002",
                    "GEDI02_A_2020076085423_O07131_03_T01780_02_003_01_V002",
                    "GEDI02_A_2020080072040_O07192_03_T03050_02_003_01_V002",
                    "GEDI02_A_2020084054801_O07253_03_T01627_02_003_01_V002",
                    "GEDI02_A_2020088041519_O07314_03_T02897_02_003_01_V002",
                    "GEDI02_A_2020092024236_O07375_03_T04167_02_003_01_V002",
                ]
            },
        },
    },
    "CA_YosemiteNP_2019": {
        "provider": "USGS",
        "datetime": ["2019-10-06", "2019-11-05"],
        "overlap_geometry": "POLYGON ((-119.5492575968371 38.06526583204277, -119.53734167552687 38.158315268657255, -119.5179123158 38.142176044699994, -119.51402521180002 38.1612131398, -119.4894348269374 38.16262019618275, -119.4619867552 38.11630312370001, -119.42909732468549 38.120153556330095, -119.41637399179986 38.06973536669182, -119.5492575968371 38.06526583204277))",
        "filters": {
            "3dep": {"workunit": "CA_YosemiteNP_2019"},
            "maxar": {"stereo_ids": ["102001008E543300", "102001008EE06B00"]},
            "icesat-2": {"granule_ids": ["ATL03_20191008182255_01810506_006_02.h5"]},
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2019279200308_O04626_02_T03812_02_003_01_V002",
                    "GEDI02_A_2019302110619_O04977_02_T00966_02_003_01_V002",
                    "GEDI02_A_2019309145009_O05088_03_T05331_02_003_01_V002",
                ]
            },
        },
    },
    "TX_DesertMountains_B1_2018": {
        "provider": "USGS",
        "datetime": ["2019-08-28", "2019-11-02"],
        "overlap_geometry": "MULTIPOLYGON (((-105.75648307892891 31.48802043744102, -105.74044466388929 31.628299278519673, -105.66532015794972 31.60785307566191, -105.6817150069322 31.462434793354177, -105.75648307892891 31.48802043744102)), ((-105.71457893360392 31.81479372478263, -105.70834640640595 31.866955289316444, -105.6454916290212 31.825676001630708, -105.64514768085348 31.78739286540503, -105.68387842249429 31.787950501271613, -105.71457893360392 31.81479372478263)))",
        "filters": {
            "3dep": {"workunit": "TX_DesertMountains_B1_2018"},
            "maxar": {"stereo_ids": ["1030010097D5CC00", "10300100984A3A00"]},
            "icesat-2": {"granule_ids": ["ATL03_20190924180456_13540406_006_02.h5"]},
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2019249064153_O04152_02_T01470_02_003_01_V002",
                    "GEDI02_A_2019257110927_O04279_03_T04122_02_003_01_V002",
                    "GEDI02_A_2019271214434_O04503_02_T04163_02_003_01_V002",
                    "GEDI02_A_2019284003705_O04691_03_T04428_02_003_01_V002",
                    "GEDI02_A_2019298111055_O04915_02_T04622_02_003_01_V002",
                    "GEDI02_A_2019306154007_O05042_03_T01429_02_003_01_V002",
                ]
            },
        },
    },
    "CO_WestCentral_2019": {
        "provider": "USGS",
        "datetime": ["2019-08-16", "2019-10-08"],
        "overlap_geometry": "MULTIPOLYGON (((-106.8087159664412 38.650074344601826, -106.73016207984253 38.712207247211836, -106.71836013228334 38.708171939366025, -106.72185460200924 38.64390009787436, -106.79085198080406 38.59952250275771, -106.8087159664412 38.650074344601826)), ((-106.80789338847443 38.82538692772677, -106.80703647763674 38.98583952996138, -106.79181281241124 38.997555854547585, -106.75071116255369 38.97004033224472, -106.72765182335152 38.78453016515039, -106.80789338847443 38.82538692772677)))",
        "filters": {
            "3dep": {"workunit": "CO_WestCentral_2019"},
            "maxar": {"stereo_ids": ["102001008BD60800", "1020010088168700"]},
            "icesat-2": {"granule_ids": ["ATL03_20191008052407_01730502_006_02.h5"]},
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2019228215125_O03836_03_T03938_02_003_01_V002",
                    "GEDI02_A_2019247142840_O04126_03_T02362_02_003_02_V002",
                    "GEDI02_A_2019255050055_O04244_02_T03888_02_003_01_V002",
                    "GEDI02_A_2019270053052_O04477_03_T03785_02_003_01_V002",
                ]
            },
        },
    },
    "WY_FEMA_East_B9_2019": {
        "provider": "USGS",
        "datetime": ["2019-07-13", "2019-09-24"],
        "overlap_geometry": "POLYGON ((-109.63272177563503 43.16570451142812, -109.68179261170016 43.510454255175134, -109.66662928596925 43.53324932573008, -109.50051092431917 43.60620231011112, -109.46371826882357 43.57911683711135, -109.48539679398097 43.182155622226006, -109.63272177563503 43.16570451142812))",
        "filters": {
            "3dep": {"workunit": "WY_FEMA_East_B9_2019"},
            "maxar": {"stereo_ids": ["1020010086976A00", "1020010089403600"]},
            "icesat-2": {
                "granule_ids": [
                    "ATL03_20190721091907_03560402_006_02.h5",
                    "ATL03_20190723212639_03940406_006_02.h5",
                    "ATL03_20190819075518_07980402_006_02.h5",
                    "ATL03_20190821200249_08360406_006_02.h5",
                    "ATL03_20190917063123_12400402_006_02.h5",
                    "ATL03_20190919183853_12780406_006_02.h5",
                ]
            },
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2019194055411_O03298_02_T00155_02_003_01_V002",
                    "GEDI02_A_2019200084834_O03393_03_T03035_02_003_02_V002",
                    "GEDI02_A_2019209000122_O03527_02_T04730_02_003_01_V002",
                    "GEDI02_A_2019211043557_O03561_03_T00847_02_003_01_V002",
                    "GEDI02_A_2019219194832_O03695_02_T03965_02_003_02_V002",
                    "GEDI02_A_2019223180803_O03756_02_T03613_02_003_01_V002",
                    "GEDI02_A_2019225224231_O03790_03_T02729_02_003_01_V002",
                    "GEDI02_A_2019242104318_O04046_02_T02343_02_003_02_V002",
                    "GEDI02_A_2019244151907_O04080_03_T02576_02_003_02_V002",
                    "GEDI02_A_2019257045811_O04275_02_T00002_02_003_01_V002",
                    "GEDI02_A_2019261032101_O04336_02_T04883_02_003_01_V002",
                    "GEDI02_A_2019267062026_O04431_03_T05575_02_003_01_V002",
                ]
            },
        },
    },
    "REDB": {
        "provider": "NEON",
        "datetime": ["2021-03-28", "2021-06-28"],
        "overlap_geometry": "POLYGON ((-111.7469617770448 40.78845879512573, -111.80059326210007 40.76843628665632, -111.79312640671719 40.82830749058919, -111.7439004944288 40.828639975426405, -111.7469617770448 40.78845879512573))",
        "filters": {
            "neon": {"site": "REDB"},
            "maxar": {
                "stereo_ids": [
                    "10200100B1301200",
                    "10200100B1172E00",
                    "10200100B1551C00",
                    "10200100B1630600",
                ]
            },
            "icesat-2": {"granule_ids": ["ATL03_20210628115237_00741206_006_01.h5"]},
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2021124140303_O13552_03_T11252_02_003_02_V002",
                    "GEDI02_A_2021087223155_O12984_02_T08647_02_003_02_V002",
                ]
            },
        },
    },
    "BART": {
        "provider": "NEON",
        "datetime": ["2019-06-10", "2019-09-22"],
        "overlap_geometry": "POLYGON ((-71.28063099480259 43.987633244325664, -71.29835822129166 44.02732951660264, -71.20408121118402 44.02583039135407, -71.20107586677035 43.98781654969867, -71.28063099480259 43.987633244325664))",
        "filters": {
            "neon": {"site": "BART"},
            "maxar": {
                "stereo_ids": [
                    "10300100988FA400",
                    "1030010098464200",
                    "1030010093B2F600",
                    "1030010092D6AF00",
                ]
            },
            "icesat-2": {
                "granule_ids": [
                    "ATL03_20190903043904_10250402_006_02.h5",
                    "ATL03_20190807181029_06210406_006_02.h5",
                ]
            },
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2019195081130_O03315_03_T05405_02_003_01_V002",
                    "GEDI02_A_2019210205202_O03556_02_T00995_02_003_01_V002",
                ]
            },
        },
    },
    "WREF": {
        "provider": "NEON",
        "datetime": ["2019-05-22", "2019-08-21"],
        "overlap_geometry": "POLYGON ((-122.00390477808023 45.779209212637085, -122.0344861746451 45.82361234666161, -121.99166117102219 45.845309058661, -121.99238741969438 45.78813979017355, -122.00390477808023 45.779209212637085))",
        "filters": {
            "neon": {"site": "WREF"},
            "maxar": {"stereo_ids": ["10300100965CBA00", "10300100962F3500"]},
            "icesat-2": {
                "granule_ids": [
                    "ATL03_20190529005739_09280306_006_02.h5",
                    "ATL03_20190522125829_08290302_006_02.h5",
                    "ATL03_20190821083817_08290402_006_02.h5",
                ]
            },
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2019206022612_O03482_02_T00370_02_003_01_V002"
                ]
            },
        },
    },
    "OTLAS.092021.32611.1": {
        "provider": "NCALM",
        "datetime": ["2020-01-16", "2020-02-23"],
        "overlap_geometry": "POLYGON ((-115.73124807944741 33.35440577353728, -115.78901631590048 33.401096326313194, -115.79498918906548 33.42186189441148, -115.74678165992853 33.38575032566509, -115.73124807944741 33.35440577353728))",
        "filters": {
            "neon": {"site": "OTLAS.092021.32611.1"},
            "maxar": {"stereo_ids": ["10300100A178D600", "10300100A3456700"]},
            "icesat-2": {
                "granule_ids": [
                    "ATL03_20200126004754_04630602_006_01.h5",
                    "ATL03_20200116132021_03180606_006_01.h5",
                ]
            },
            "gedi": {
                "granule_ids": [
                    "GEDI02_A_2020050133846_O06731_02_T04240_02_003_01_V002",
                    "GEDI02_A_2020046151221_O06670_02_T01394_02_003_01_V002",
                    "GEDI02_A_2020054120506_O06792_02_T01547_02_003_01_V002",
                    "GEDI02_A_2020034195258_O06487_02_T01547_02_003_01_V002",
                    "GEDI02_A_2020026225949_O06365_02_T01394_02_003_01_V002",
                ]
            },
        },
    },
}


# def load_pcd_geometries() -> gpd.GeoDataFrame:
#     """
#     Loads PCD site geometries by downloading them from the latest GitHub release.

#     This function finds the 'pcd_overlap_geometries_2025.parquet' asset from
#     the latest release, downloads it into memory, and loads it into a GeoDataFrame.

#     Returns:
#         gpd.GeoDataFrame: A GeoDataFrame with 'pcd_id' and 'geometry' columns.

#     Raises:
#         FileNotFoundError: If the geometry asset is not found in the latest release.
#     """
#     OWNER = "uw-cryo"
#     REPO = "coincident"
#     GEOMETRY_FILENAME = "pcd_overlap_geometries_2025.parquet"

#     logging.info("Fetching latest release to find geometry file...")
#     api_url = f"https://api.github.com/repos/{OWNER}/{REPO}/releases/latest"

#     try:
#         response = requests.get(api_url)
#         response.raise_for_status()
#         assets = response.json().get("assets", [])
#     except requests.exceptions.RequestException as e:
#         msg_git_release = f"Failed to fetch GitHub release info: {e}"
#         logging.error(msg_git_release)
#         raise

#     asset_url = next(
#         (
#             asset["browser_download_url"]
#             for asset in assets
#             if asset["name"] == GEOMETRY_FILENAME
#         ),
#         None,
#     )

#     if not asset_url:
#         msg_no_asset = (
#             f"Asset '{GEOMETRY_FILENAME}' not found in the latest GitHub release."
#         )
#         raise FileNotFoundError(msg_no_asset)

#     # logging.info(f"Downloading geometry file: {GEOMETRY_FILENAME}")

#     try:
#         response = requests.get(asset_url)
#         response.raise_for_status()

#         # read downloaded content directly from an in-memory buffer
#         buffer = io.BytesIO(response.content)
#         return gpd.read_parquet(buffer)

#     except requests.exceptions.RequestException as e:
#         msg_no_download = f"Failed to download geometry file: {e}"
#         logging.error(msg_no_download)
#         raise


# NOTE: could adjust the below to do a cascading_search but this works easier with how
# I set up the PCD_SITES dict and also the cascading_search automatically dissolves
# maxar stereo where this preserves both images
def read_pcd_site(pcd_id: str) -> dict[str, gpd.GeoDataFrame]:
    """
    Performs on-demand searches to retrieve all data for a given PCD site.

    This function now uses the site's hardcoded 'overlap_geometry' from the
    PCD_SITES fixtures to run live searches.
    """
    if pcd_id not in PCD_SITES:
        msg_invalid_id = f"PCD Site ID '{pcd_id}' not found in fixtures."
        raise ValueError(msg_invalid_id)

    msg_site_id = f"Reading data for PCD Site: {pcd_id}"
    logging.info(msg_site_id)
    site_params = PCD_SITES[pcd_id]
    search_date = list(site_params["datetime"])
    provider = site_params["provider"]
    filters = site_params["filters"]
    wkt_geom = site_params.get("overlap_geometry")

    if not wkt_geom:
        msg_no_geom = f"Overlap geometry not found for PCD Site ID '{pcd_id}'."
        raise ValueError(msg_no_geom)

    # Create gdf from the hardcoded WKT string in our PCD site class
    geom = wkt.loads(wkt_geom)
    gf_overlap_original = gpd.GeoDataFrame(
        {"pcd_id": [pcd_id], "geometry": [geom]}, crs="EPSG:4326"
    )

    # less complex search polygon to improve performance and rid errors
    gf_overlap_search = gf_overlap_original.copy()
    gf_overlap_search.geometry = gf_overlap_search.geometry.envelope

    # Search for lidar
    ms_als_prov = f"Searching for {provider} data..."
    logging.info(ms_als_prov)
    if provider == "USGS":
        gf_als = search(
            dataset="3dep", intersects=gf_overlap_search, datetime=search_date
        )
        gf_als = gf_als[gf_als.workunit == pcd_id]
        if not gf_als.empty:
            gf_als = wesm.load_by_fid(fids=gf_als.index.values)
    elif provider == "NEON":
        gf_als = search(
            dataset="neon", intersects=gf_overlap_search, datetime=search_date
        )
        gf_als = gf_als[gf_als.id == pcd_id].reset_index(drop=True)
    elif provider == "NCALM":
        gf_als = search(
            dataset="ncalm", intersects=gf_overlap_search, datetime=search_date
        )
        gf_als = gf_als[gf_als.id == pcd_id].reset_index(drop=True)
    else:
        gf_als = gpd.GeoDataFrame()

    msg_maxar = "Searching for Maxar data..."
    logging.info(msg_maxar)
    gf_maxar = search(
        dataset="maxar", intersects=gf_overlap_search, datetime=search_date
    )
    gf_maxar = gf_maxar[gf_maxar.id.isin(filters["maxar"]["stereo_ids"])]

    msg_gedi = "Searching for GEDI data..."
    logging.info(msg_gedi)
    gf_gedi = search(dataset="gedi", intersects=gf_overlap_search, datetime=search_date)
    gf_gedi = gf_gedi[gf_gedi.id.isin(filters["gedi"]["granule_ids"])]

    msg_is2 = "Searching for ICESat-2 data..."
    logging.info(msg_is2)
    gf_is2 = search(
        dataset="icesat-2", intersects=gf_overlap_search, datetime=search_date
    )
    gf_is2 = gf_is2[gf_is2.id.isin(filters["icesat-2"]["granule_ids"])]

    return {
        "als": gf_als,
        "maxar": gf_maxar,
        "is2": gf_is2,
        "gedi": gf_gedi,
        "overlap": gf_overlap_original,
    }


def download_pcd_files(
    output_dir: str, geoparquet: bool = True, geojson: bool = True
) -> None:
    """
    Downloads all data for every PCD site, groups it by provider and sensor (als, stereo, gedi, icesat-2, overlap area),
    and saves the results to GeoParquet and/or GeoJSON files.

    Args:
        output_dir (str): The directory where output files will be saved.
        geoparquet (bool): If True, saves files in GeoParquet format. Defaults to True.
        geojson (bool): If True, saves files in GeoJSON format. Defaults to True.
    """
    if not (geoparquet or geojson):
        msg = "No output format selected. Set geoparquet or geojson to True."
        logging.warning(msg)
        return

    # Create output dir if it doesn't exist
    try:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        msg_dir = f"Output directory '{output_dir}' is ready."
        logging.info(msg_dir)
    except OSError as e:
        msg_dir_fail = f"Failed to create output directory: {e}"
        logging.error(msg_dir_fail)
        raise

    # Group sites by provider (USGS, NEON, NCALM)
    sites_by_provider: dict[str, list[str]] = {}
    for pcd_id, params in PCD_SITES.items():
        provider = params.get("provider", "UNKNOWN")
        if provider not in sites_by_provider:
            sites_by_provider[provider] = []
        sites_by_provider[provider].append(pcd_id)

    # iterate over each provider (USGS, NEON, NCALM) group
    for provider, pcd_ids in sites_by_provider.items():
        msg_provider = f"Processing provider: {provider}"
        logging.info(msg_provider)

        # in-memory stores for GeoDataFrames for the current provider
        provider_data: dict[str, list[gpd.GeoDataFrame]] = {
            "als": [],
            "maxar": [],
            "is2": [],
            "gedi": [],
            "overlap": [],
        }

        for pcd_id in pcd_ids:
            try:
                site_data = read_pcd_site(pcd_id)
                for sensor, gdf in site_data.items():
                    if not gdf.empty:
                        # Add the primary PCD ID to each dataframe for consistent tracking
                        gdf["pcd_id"] = pcd_id
                        provider_data[sensor].append(gdf)
            except Exception as e:
                msg_fail = f"Failed to download data for site '{pcd_id}'. Error: {e}"
                logging.error(msg_fail)
                continue  # Move to the next site

        # Write the collected data for the provider to files
        for sensor, gdf_list in provider_data.items():
            if not gdf_list:
                msg_no_data = (
                    f"No data found for sensor '{sensor}' under provider '{provider}'."
                )
                logging.info(msg_no_data)
                continue

            msg_concat = (
                f"Concatenating {len(gdf_list)} GeoDataFrames for {provider}-{sensor}."
            )
            logging.info(msg_concat)
            full_gdf = pd.concat(gdf_list, ignore_index=True)

            # provider-specific ID column name
            id_col_map = {"USGS": "workunit", "NEON": "neon_id", "NCALM": "ncalm_id"}
            id_col_name = id_col_map.get(provider, "pcd_id")

            if id_col_name in full_gdf.columns and id_col_name != "pcd_id":
                # 'workunit' already exists from the wesm
                # so 'pcd_id' column is now redundant, so we drop it.
                full_gdf = full_gdf.drop(columns=["pcd_id"])
            else:
                # The provider-specific column does not exist (e.g., for maxar, gedi, anything not usgs als),
                # rename our 'pcd_id' column to create it.
                full_gdf = full_gdf.rename(columns={"pcd_id": id_col_name})

            # ID column first column
            cols = [id_col_name] + [
                col for col in full_gdf.columns if col != id_col_name
            ]
            full_gdf = full_gdf[cols]

            # quick patch for duplicate 'id' col in NEON and NCALM ALS gdfs because I don't have time to implement a better fix
            if sensor == "als" and "id" in full_gdf.columns and id_col_name != "id":
                full_gdf = full_gdf.drop(columns=["id"])

            base_path = Path(output_dir) / f"PCD_{provider}_{sensor}"
            if geoparquet:
                try:
                    pq_path = f"{base_path}.parquet"
                    msg_write = f"Writing file: {pq_path}"
                    logging.info(msg_write)
                    full_gdf.to_parquet(pq_path)
                except Exception as e:
                    msg_err = f"Failed to write GeoParquet file {pq_path}: {e}"
                    logging.error(msg_err)

            if geojson:
                try:
                    json_path = f"{base_path}.geojson"
                    msg_write = f"Writing file: {json_path}"
                    logging.info(msg_write)
                    full_gdf.to_file(json_path, driver="GeoJSON")
                except Exception as e:
                    msg_err = f"Failed to write GeoJSON file {json_path}: {e}"
                    logging.error(msg_err)

    msg_done = "All providers processed. File download complete."
    logging.info(msg_done)
