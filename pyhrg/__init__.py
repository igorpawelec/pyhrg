"""
pyHRG — individual tree crown delineation from canopy height models.

Crowns are delineated by Hierarchical Region Growing: a marker-based
watershed seeded with tree tops, a weighted region adjacency graph over
the resulting regions, and greedy variance-constrained growing that lets
neighbouring trees merge — which is how over-detected tree tops get
corrected.

Quick start::

    from pyhrg import CrownDelineator

    cd = CrownDelineator.from_file("chm.tif")
    cd.smooth(ws=3).detect(hmin=7, ws=5).merge(5.0).screen(10.0)
    crowns = cd.delineate(variance_thresh=2.0)
    cd.to_raster("crowns.tif")

or, in one call::

    from pyhrg import delineate_crowns
    crowns, tops = delineate_crowns("chm.tif", hmin=7, merge_distance=5.0)

Every stage is also a plain function on arrays::

    from pyhrg import smooth_chm, detect_tops, HierarchicalRegionGrower

Copyright (C) 2025 Igor Pawelec
Licence: GPLv3 — see LICENSE.
"""

try:
    from importlib.metadata import version, PackageNotFoundError
    try:
        __version__ = version("pyhrg")
    except PackageNotFoundError:  # not installed, e.g. running from source
        __version__ = "0.2.0"
except ImportError:  # pragma: no cover
    __version__ = "0.2.0"

__author__ = "Igor Pawelec"
__license__ = "GPLv3"

_DEPS_HINT = (
    "conda install -c conda-forge numpy numba scipy scikit-image "
    "rasterio fiona"
)

# Heavy dependencies (numba, skimage) are pulled in on first use rather
# than at import time, so `import pyhrg` stays fast and a broken optional
# dependency surfaces where it is actually needed.
_LAZY = {
    "smooth_chm": ".chm",
    "SMOOTHING_METHODS": ".chm",
    "detect_tops": ".treetops",
    "merge_tops": ".treetops",
    "screen_tops": ".treetops",
    "as_pixels": ".treetops",
    "HierarchicalRegionGrower": ".hrg",
    "resolve_conflicts": ".hrg",
    "CONFLICT_RULES": ".hrg",
    "CrownDelineator": ".delineate",
    "delineate_crowns": ".delineate",
    "read_chm": ".io",
    "save_crowns_raster": ".io",
    "save_segments": ".io",
    "save_tree_tops": ".io",
}


def __getattr__(name):
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from importlib import import_module
    try:
        mod = import_module(module, __name__)
    except (ImportError, OSError) as e:
        raise ImportError(
            f"cannot import {name!r} from pyhrg{module}: {e}\n"
            f"Install the dependencies with:\n  {_DEPS_HINT}"
        ) from e
    value = getattr(mod, name)
    globals()[name] = value  # cache, so this runs once per name
    return value


def __dir__():
    return sorted(list(globals()) + list(_LAZY))


__all__ = sorted(_LAZY)
