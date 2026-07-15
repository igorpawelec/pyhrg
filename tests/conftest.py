"""Shared fixtures."""
import numpy as np
import pytest


def _gauss(shape, cy, cx, height, sigma):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return height * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sigma ** 2))


@pytest.fixture
def stand():
    """A small synthetic stand: 12 well-separated trees on a 120x120 CHM."""
    rng = np.random.default_rng(42)
    chm = np.zeros((120, 120))
    tops = []
    for r in range(15, 106, 30):
        for c in range(15, 106, 30):
            h = float(rng.uniform(14, 28))
            chm = np.maximum(chm, _gauss(chm.shape, r, c, h, 6.0))
            tops.append((r, c))
    return chm.astype(np.float32), tops


@pytest.fixture
def two_tops_one_tree():
    """One broad crown carrying two local maxima — classic over-detection.

    The right-hand peak is the taller one.
    """
    shape = (60, 60)
    chm = _gauss(shape, 30, 30, 24.0, 11.0)
    chm += _gauss(shape, 30, 24, 1.5, 3.0)
    chm += _gauss(shape, 30, 36, 2.0, 3.0)
    return chm.astype(np.float32), [(30, 24), (30, 36)]


@pytest.fixture
def flat():
    """A CHM with no trees at all."""
    return np.full((40, 40), 1.0, dtype=np.float32)
