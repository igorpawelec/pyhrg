"""
pyhrg.chm — canopy height model smoothing.

A raw CHM carries pit noise and spurious local maxima: a single laser
return through a gap, or the ragged upper surface of a crown, both create
peaks that tree detection would read as separate trees. Smoothing trades
a little height accuracy for far fewer false tops.

Copyright (C) 2025 Igor Pawelec
Licence: GPLv3 — see LICENSE.
"""

import numpy as np

__all__ = ["smooth_chm", "SMOOTHING_METHODS"]

SMOOTHING_METHODS = ("median", "mean", "gaussian", "maximum")


def _ensure_scipy():
    try:
        from scipy import ndimage
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "scipy is required for CHM smoothing.\n"
            "  conda install -c conda-forge scipy"
        ) from e
    return ndimage


def smooth_chm(chm, ws=3, method="median"):
    """
    Smooth a canopy height model.

    Parameters
    ----------
    chm : ndarray, shape (rows, cols)
        Canopy height model.
    ws : int
        Window size in pixels. Larger windows suppress more false tops but
        also merge genuinely adjacent crowns. Should be smaller than the
        narrowest crown you want to keep.
    method : {'median', 'mean', 'gaussian', 'maximum'}
        - ``'median'`` (default) — removes pits and spikes while keeping
          crown edges sharp. The usual choice.
        - ``'mean'`` — cheaper, but blurs crown boundaries.
        - ``'gaussian'`` — smooth falloff, sigma = ws / 3.
        - ``'maximum'`` — dilates crowns; flattens the apex, so tops are
          detected as plateaus rather than points. Rarely what you want
          before :func:`pyhrg.treetops.detect_tops`.

    Returns
    -------
    smoothed : ndarray, same shape as *chm*

    Notes
    -----
    NaNs are not handled specially: scipy's filters propagate them, so a
    nodata-as-NaN raster will grow its NaN halo by ``ws // 2``. Fill or
    mask nodata before smoothing if that matters.
    """
    ndimage = _ensure_scipy()
    chm = np.asarray(chm)
    if chm.ndim != 2:
        raise ValueError(f"chm must be 2-D, got shape {chm.shape}")
    if ws < 1:
        raise ValueError(f"ws must be >= 1, got {ws}")
    if method not in SMOOTHING_METHODS:
        raise ValueError(
            f"method must be one of {SMOOTHING_METHODS}, got {method!r}"
        )

    if method == "median":
        return ndimage.median_filter(chm, size=ws)
    if method == "mean":
        return ndimage.uniform_filter(chm, size=ws)
    if method == "gaussian":
        return ndimage.gaussian_filter(chm, sigma=ws / 3.0)
    return ndimage.maximum_filter(chm, size=ws)
