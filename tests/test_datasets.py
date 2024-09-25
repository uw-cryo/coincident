from __future__ import annotations

import coincident


def test_valid_data_sources():
    supported = {"maxar", "nasa", "usgs"}
    imported = set(dir(coincident.datasets))
    assert imported & supported == supported


def test_maxar_defaults():
    ds = coincident.datasets.maxar.Stereo()
    mono_collections = ["wv01", "wv02", "wv03-vnir", "ge01"]
    assert set(ds.collections) == set(mono_collections)


def test_threedep_defaults():
    ds = coincident.datasets.usgs.ThreeDEP()
    assert ds.alias == "3dep"
