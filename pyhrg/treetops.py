"""
pyhrg.treetops — tree top detection, merging and screening.

Tree tops are found as local maxima of a smoothed CHM, the standard
raster approach to individual tree detection (Popescu & Wynne 2004). The
window size sets the implicit minimum spacing between trees, so it is the
main lever on over- versus under-detection.

Every function here returns tree tops as an ``(n, 2)`` float64 array of
(row, column) coordinates — including the empty case, which is ``(0, 2)``
rather than ``(0,)`` so that callers can index ``[:, 0]`` unconditionally.

References
----------
Popescu, S.C., Wynne, R.H. (2004). Seeing the trees in the forest: using
    lidar and multispectral data fusion with local filtering and variable
    window size for estimating tree height. Photogrammetric Engineering &
    Remote Sensing 70(5), 589-604.

Copyright (C) 2025 Igor Pawelec
Licence: GPLv3 — see LICENSE.
"""

import numpy as np

__all__ = ["detect_tops", "merge_tops", "screen_tops", "as_pixels"]


def _ensure_scipy():
    try:
        from scipy import ndimage
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "scipy is required for tree top detection.\n"
            "  conda install -c conda-forge scipy"
        ) from e
    return ndimage


def _as_tops_array(x):
    """Normalise anything top-like to float64 (n, 2), empty -> (0, 2)."""
    return np.asarray(x, dtype=np.float64).reshape(-1, 2)


def detect_tops(chm, hmin=2.0, ws=3):
    """
    Detect tree tops as local maxima of a (smoothed) CHM.

    A pixel is a candidate when it equals the maximum of its ``ws x ws``
    neighbourhood and exceeds *hmin*. Adjacent candidates — a flat apex
    yields several — are labelled into connected components and each
    component is reduced to its height-weighted centre of mass, giving
    subpixel coordinates.

    Parameters
    ----------
    chm : ndarray, shape (rows, cols)
        Canopy height model. Smooth it first (:func:`pyhrg.chm.smooth_chm`);
        on a raw CHM this will over-detect badly.
    hmin : float
        Minimum height for a pixel to be considered. Sets the floor
        between canopy and understorey/ground.
    ws : int
        Neighbourhood size in pixels. Acts as the minimum spacing between
        detected tops: too small over-detects one tree as many, too large
        merges neighbouring trees. Over-detection is the safer error here,
        since :func:`pyhrg.hrg.HierarchicalRegionGrower.run_all` can merge
        surplus tops back together.

    Returns
    -------
    tops : ndarray, shape (n, 2), dtype float64
        Subpixel (row, column) coordinates. Empty gives shape (0, 2).
    """
    ndimage = _ensure_scipy()
    chm = np.asarray(chm)
    if chm.ndim != 2:
        raise ValueError(f"chm must be 2-D, got shape {chm.shape}")
    if ws < 1:
        raise ValueError(f"ws must be >= 1, got {ws}")

    local_max = ndimage.maximum_filter(chm, size=ws)
    detected = (chm == local_max) & (chm > hmin)
    labels, num = ndimage.label(detected)
    centres = ndimage.center_of_mass(chm, labels, range(1, num + 1))
    return _as_tops_array(centres)


def merge_tops(tops, distance=5.0):
    """
    Merge tree tops that lie within *distance* of each other.

    Tops closer than the threshold are almost always the same tree seen
    through several local maxima. Neighbours within the radius are found
    with a k-d tree, grouped transitively with union-find, and each group
    is replaced by its centroid.

    Note the grouping is transitive: a chain of tops each within
    *distance* of the next collapses into one, even if the endpoints are
    far apart. Keep the threshold well below the typical crown diameter.

    Parameters
    ----------
    tops : array-like, shape (n, 2)
        Tree top coordinates.
    distance : float
        Merge radius in pixels.

    Returns
    -------
    merged : ndarray, shape (m, 2), m <= n
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "scipy is required for merging tree tops.\n"
            "  conda install -c conda-forge scipy"
        ) from e

    tops = _as_tops_array(tops)
    if tops.shape[0] < 2:
        return tops

    kd = cKDTree(tops)
    neighbour_lists = kd.query_ball_point(tops, r=distance)

    n = tops.shape[0]
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i, neighbours in enumerate(neighbour_lists):
        for j in neighbours:
            if j > i:
                union(i, j)

    from collections import defaultdict
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    merged = [tops[idx].mean(axis=0) for idx in groups.values()]
    return _as_tops_array(merged)


def screen_tops(chm, tops, hmin):
    """
    Drop tree tops standing on CHM values below *hmin*.

    Useful as a second pass with a stricter threshold than
    :func:`detect_tops` used — e.g. detect at 2 m to let the growing see
    the understorey, then keep only tops above 10 m.

    Parameters
    ----------
    chm : ndarray, shape (rows, cols)
        The same CHM the tops were detected on.
    tops : array-like, shape (n, 2)
    hmin : float
        Minimum height to keep.

    Returns
    -------
    kept : ndarray, shape (m, 2), m <= n
    """
    chm = np.asarray(chm)
    tops = _as_tops_array(tops)
    if tops.shape[0] == 0:
        return tops

    rows, cols = chm.shape
    rc = np.floor(tops).astype(np.intp)
    inside = ((rc[:, 0] >= 0) & (rc[:, 0] < rows) &
              (rc[:, 1] >= 0) & (rc[:, 1] < cols))
    keep = np.zeros(len(tops), dtype=bool)
    keep[inside] = chm[rc[inside, 0], rc[inside, 1]] >= hmin
    return tops[keep]


def as_pixels(tops):
    """
    Convert subpixel tops to integer (row, col) pixel indices.

    Returns
    -------
    pixels : list of (int, int)
        Suitable for :meth:`pyhrg.hrg.HierarchicalRegionGrower.run_all`.
    """
    tops = _as_tops_array(tops)
    return [(int(r), int(c)) for r, c in np.floor(tops)]
