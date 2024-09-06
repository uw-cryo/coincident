from __future__ import annotations

import importlib.metadata

import coincident as m


def test_version():
    assert importlib.metadata.version("coincident") == m.__version__
