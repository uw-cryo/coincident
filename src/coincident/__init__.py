"""
Copyright (c) 2024 Scott Henderson. All rights reserved.

coincident: Search and analysis of STV Precursor Coincident Datasets
"""

from __future__ import annotations

from coincident import datasets, overlaps, search
from coincident._version import version as __version__

__all__ = ["__version__", "datasets", "search", "overlaps"]
