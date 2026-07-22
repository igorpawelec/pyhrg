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
    # An even window has no centre pixel, so scipy places it half a pixel
    # off, and the result then depends on which way the raster happens to
    # be oriented: smoothing a scene and smoothing its mirror image gave
    # values up to 8.8 m apart on a 40x55 test scene. Nobody picks ws=4
    # intending a shifted window, so this is refused rather than warned
    # about. 'gaussian' is exempt because ws only scales sigma there and
    # the kernel stays symmetric.
    if method != "gaussian" and ws % 2 == 0:
        raise ValueError(
            f"ws must be odd for method='{method}', got {ws}. An even window "
            f"has no centre pixel, so it sits half a pixel off and the "
            f"result depends on the raster's orientation. Use {ws - 1} or "
            f"{ws + 1}."
        )

    # NaN is skipped, not propagated and not fed to the filters raw.
    #
    # scipy's comparison-based filters are undefined on NaN rather than
    # NaN-aware: the result depends on where in the window the NaN sits.
    # With the nine values 1..9 and one NaN, median_filter returned 9, 5 or
    # 4 depending on the position, and maximum_filter returned nan, 9 or 8.
    # On chm_150_2023.tif that left 1208 pixels (3% of the raster) smoothed
    # to an arbitrary window element, 20 of 253 tree tops standing on them,
    # and errors against a NaN-skipping median of up to 26 m.
    #
    # The fast path is taken whenever there is no NaN, so a clean CHM costs
    # nothing and gives byte-identical results to before.
    if not np.isnan(chm).any():
        if method == "median":
            return ndimage.median_filter(chm, size=ws)
        if method == "mean":
            return ndimage.uniform_filter(chm, size=ws)
        if method == "gaussian":
            return ndimage.gaussian_filter(chm, sigma=ws / 3.0)
        return ndimage.maximum_filter(chm, size=ws)

    return _smooth_with_nan(ndimage, chm, ws, method)


def _smooth_with_nan(ndimage, chm, ws, method):
    """Window statistics over the non-NaN cells; all-NaN windows stay NaN."""
    valid = ~np.isnan(chm)
    filled = np.where(valid, chm, 0.0)
    w = valid.astype(np.float64)

    if method == "maximum":
        # -inf loses every comparison, so it drops out of the max without
        # poisoning it the way NaN does.
        out = ndimage.maximum_filter(np.where(valid, chm, -np.inf), size=ws)
    elif method == "mean":
        total = ndimage.uniform_filter(filled, size=ws)
        count = ndimage.uniform_filter(w, size=ws)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = total / count
    elif method == "gaussian":
        # Normalised convolution: weight the kernel by which cells exist.
        sigma = ws / 3.0
        total = ndimage.gaussian_filter(filled, sigma=sigma)
        count = ndimage.gaussian_filter(w, sigma=sigma)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = total / count
    else:                                    # median
        # No vectorised form; 178x slower than median_filter on a 200x201
        # raster, which is why this path runs only when NaN is present.
        out = ndimage.generic_filter(chm, np.nanmedian, size=ws,
                                     mode="reflect")

    # A window with nothing in it has no statistic. For the gaussian that
    # means the *kernel* reached nothing, and the kernel is wider than ws --
    # 4 sigma either side -- so the test is the denominator rather than a
    # ws-sized window. Using ws there would have marked 7663 cells empty
    # against the 5889 the kernel actually fails to reach, and put rHRG and
    # this package 1774 pixels apart.
    if method == "gaussian":
        empty = count == 0
    else:
        empty = ndimage.maximum_filter(valid.astype(np.uint8), size=ws) == 0
    out = np.asarray(out, dtype=float)
    out[empty] = np.nan
    out[~np.isfinite(out) & ~empty] = np.nan
    return out
