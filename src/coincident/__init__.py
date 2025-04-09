"""
coincident: Search and analysis of STV Precursor Coincident Datasets
"""

from __future__ import annotations

from coincident import datasets, io, overlaps, plot, search
from coincident._version import version as __version__

__all__ = ["__version__", "datasets", "io", "overlaps", "plot", "search"]
