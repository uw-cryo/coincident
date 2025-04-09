"""Helper functions for plot testing"""

from __future__ import annotations

import matplotlib.patches
import matplotlib.pyplot as plt
import pytest


@pytest.fixture
def assert_boxplot():
    """Fixture for common boxplot assertions"""

    def _assert_boxplot(axd, expected_boxes):
        """Helper to verify boxplot properties"""
        assert isinstance(axd, plt.Axes), (
            "Return value should be a matplotlib Axes object"
        )
        assert len(axd.get_children()) > 0, "Figure should exist"
        assert (
            len(
                [
                    child
                    for child in axd.get_children()
                    if isinstance(child, matplotlib.patches.PathPatch)
                ]
            )
            == expected_boxes
        ), f"Should have {expected_boxes} boxplot boxes"
        assert axd.get_legend() is not None, "Legend should exist"
        assert len(axd.figure.axes) == 2, "Should have two y-axes"
        assert axd.figure.axes[1].get_ylabel() == "Count", (
            "Second y-axis should be 'Count'"
        )

    return _assert_boxplot
