from __future__ import annotations

import coincident as m


def test_version():
    # Just check that conventional __version__ exists
    assert m.__version__
    # NOTE: these do not necessarily match during development I think
    # importlib picks up the version of package at time of install
    # but _version.py / __init__.py is only updated by `pipx run build`
    # https://github.com/python/cpython/issues/94181
    # assert importlib.metadata.version("coincident") == m.__version__
