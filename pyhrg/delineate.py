"""
pyhrg.delineate — the end-to-end pipeline.

:class:`CrownDelineator` chains smoothing, tree top detection and
hierarchical region growing, holding the intermediate state so the steps
can be inspected and re-run individually. It is a convenience layer: each
stage is an ordinary function in :mod:`pyhrg.chm`, :mod:`pyhrg.treetops`
and :mod:`pyhrg.hrg`, and can be used on its own.

Copyright (C) 2025 Igor Pawelec
Licence: GPLv3 — see LICENSE.
"""

import time

import numpy as np

from .chm import smooth_chm
from .treetops import detect_tops, merge_tops, screen_tops, as_pixels
from .hrg import HierarchicalRegionGrower

__all__ = ["CrownDelineator", "delineate_crowns"]


class CrownDelineator:
    """
    Delineate individual tree crowns from a canopy height model.

    The pipeline is::

        smooth  ->  detect_tops  ->  merge_tops  ->  screen_tops  ->  delineate

    Every step is optional and re-runnable; :meth:`delineate` fills in
    smoothing and detection with defaults if they were skipped.

    Parameters
    ----------
    chm : ndarray, shape (rows, cols)
        Canopy height model. Use :meth:`from_file` to read one from disk.
    transform : affine.Affine, optional
        Geotransform, needed only for the export helpers.
    crs : rasterio.crs.CRS or str, optional
        Coordinate reference system, needed only for the export helpers.
    quiet : bool
        Suppress progress messages.

    Attributes
    ----------
    chm : ndarray
        The input CHM, unmodified.
    smoothed : ndarray or None
        Result of :meth:`smooth`.
    tops : ndarray, shape (n, 2) or None
        Current tree tops, in subpixel (row, col).
    crowns : ndarray or None
        Crown id per pixel, 0 = background.

    Examples
    --------
    >>> cd = CrownDelineator.from_file("chm.tif")          # doctest: +SKIP
    >>> cd.smooth(ws=3).detect(hmin=7, ws=5).merge(5.0)    # doctest: +SKIP
    >>> crowns = cd.delineate(variance_thresh=2.0)         # doctest: +SKIP
    """

    def __init__(self, chm, transform=None, crs=None, quiet=False):
        chm = np.asarray(chm)
        if chm.ndim != 2:
            raise ValueError(f"chm must be 2-D, got shape {chm.shape}")

        self.chm = chm.astype(np.float32, copy=False)
        self.transform = transform
        self.crs = crs
        self.quiet = quiet

        self.smoothed = None
        self.tops = None
        self.crowns = None
        self.grower = None

        self._say(f"CHM {self.chm.shape[0]}x{self.chm.shape[1]}, "
                  f"range {np.nanmin(self.chm):.1f}-{np.nanmax(self.chm):.1f} m")

    # ── construction ─────────────────────────────────────────────────

    @classmethod
    def from_file(cls, path, window=None, band=1, quiet=False):
        """
        Read a CHM raster and return a delineator for it.

        Parameters
        ----------
        path : str or Path
        window : tuple, optional
            ``(col_off, row_off, width, height)`` to work on a sub-region.
        band : int
        quiet : bool

        Returns
        -------
        CrownDelineator
        """
        from .io import read_chm
        chm, transform, crs = read_chm(path, window=window, band=band)
        return cls(chm, transform=transform, crs=crs, quiet=quiet)

    def _say(self, msg):
        if not self.quiet:
            print(f"  {msg}")

    # ── pipeline steps ───────────────────────────────────────────────

    def smooth(self, ws=3, method="median"):
        """Smooth the CHM. See :func:`pyhrg.chm.smooth_chm`. Returns self."""
        self.smoothed = smooth_chm(self.chm, ws=ws, method=method)
        self._say(f"smooth: {method} ws={ws}, range "
                  f"{np.nanmin(self.smoothed):.1f}-{np.nanmax(self.smoothed):.1f} m")
        return self

    def detect(self, hmin=2.0, ws=3):
        """Detect tree tops. See :func:`pyhrg.treetops.detect_tops`. Returns self."""
        if self.smoothed is None:
            self.smooth(ws=3)
        self.tops = detect_tops(self.smoothed, hmin=hmin, ws=ws)
        self._say(f"detect: {len(self.tops)} tops (hmin={hmin}, ws={ws})")
        return self

    def merge(self, distance=5.0):
        """Merge nearby tops. See :func:`pyhrg.treetops.merge_tops`. Returns self."""
        if self.tops is None:
            raise ValueError("call detect() before merge()")
        before = len(self.tops)
        self.tops = merge_tops(self.tops, distance=distance)
        self._say(f"merge: {before} -> {len(self.tops)} (distance={distance})")
        return self

    def screen(self, hmin):
        """Drop short tops. See :func:`pyhrg.treetops.screen_tops`. Returns self."""
        if self.tops is None:
            raise ValueError("call detect() before screen()")
        before = len(self.tops)
        self.tops = screen_tops(self.smoothed, self.tops, hmin=hmin)
        self._say(f"screen: {before} -> {len(self.tops)} (>={hmin} m)")
        return self

    def delineate(self,
                  variance_thresh=2.0,
                  mask_thresh=0.0,
                  morpho_radius=0,
                  alpha=1.0,
                  beta=0.5,
                  gamma=0.1,
                  anneal_lambda=1.0,
                  max_iters=200,
                  conflict_rule="height",
                  protect_seeds=False,
                  retry_rejected=False,
                  n_jobs=1):
        """
        Grow crowns from the tree tops.

        Thin wrapper over
        :meth:`pyhrg.hrg.HierarchicalRegionGrower.run_all` — see that
        method for the full meaning of each parameter.

        Returns
        -------
        crowns : ndarray[int32]
            Crown id per pixel, 0 = background. Ids follow tree-top order;
            merged trees leave gaps, so the crown count can be lower than
            the tree-top count.
        """
        if self.smoothed is None:
            self.smooth(ws=3)
        if self.tops is None:
            self.detect()

        seeds = as_pixels(self.tops)
        self._say(f"delineate: {len(seeds)} seeds, rule={conflict_rule}")

        t0 = time.time()
        self.grower = HierarchicalRegionGrower(self.smoothed)
        self.crowns = self.grower.run_all(
            seeds,
            variance_thresh=variance_thresh,
            mask_thresh=mask_thresh,
            morpho_radius=morpho_radius,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            anneal_lambda=anneal_lambda,
            max_iters=max_iters,
            conflict_rule=conflict_rule,
            protect_seeds=protect_seeds,
            retry_rejected=retry_rejected,
            n_jobs=n_jobs,
        ).astype(np.int32)

        n_crowns = len(np.unique(self.crowns)) - 1
        merged = len(seeds) - n_crowns
        extra = f", {merged} merged" if merged > 0 else ""
        self._say(f"delineate: {n_crowns} crowns{extra} "
                  f"in {time.time() - t0:.2f}s "
                  f"({self.grower.n_contested} contested regions)")
        return self.crowns

    # ── export ───────────────────────────────────────────────────────

    def _need_geo(self, what):
        if self.transform is None:
            raise ValueError(
                f"{what} needs a geotransform; construct with "
                f"CrownDelineator.from_file(...) or pass transform="
            )

    def to_raster(self, path, **kwargs):
        """Write the crown raster as a GeoTIFF."""
        from .io import save_crowns_raster
        if self.crowns is None:
            raise ValueError("call delineate() first")
        self._need_geo("to_raster")
        crs_wkt = self.crs.to_wkt() if hasattr(self.crs, "to_wkt") else self.crs
        save_crowns_raster(self.crowns, str(path), self.transform,
                           crs_wkt, **kwargs)
        self._say(f"wrote {path}")
        return path

    def to_vector(self, out_dir, name="crowns", **kwargs):
        """Write crown polygons (Shapefile by default)."""
        from .io import save_segments
        if self.crowns is None:
            raise ValueError("call delineate() first")
        self._need_geo("to_vector")
        crs_wkt = self.crs.to_wkt() if hasattr(self.crs, "to_wkt") else self.crs
        save_segments(self.crowns, str(out_dir), name, self.transform,
                      crs_wkt, chm_array=self.chm, **kwargs)
        self._say(f"wrote {out_dir}/{name}")
        return out_dir

    def tops_to_vector(self, out_dir, name="treetops", **kwargs):
        """Write tree tops as points."""
        from .io import save_tree_tops
        if self.tops is None:
            raise ValueError("call detect() first")
        self._need_geo("tops_to_vector")
        crs_wkt = self.crs.to_wkt() if hasattr(self.crs, "to_wkt") else self.crs
        save_tree_tops(self.tops, str(out_dir), name, self.transform,
                       crs_wkt, chm=self.chm, **kwargs)
        self._say(f"wrote {out_dir}/{name}")
        return out_dir


def delineate_crowns(chm,
                     smooth_ws=3,
                     smooth_method="median",
                     hmin=2.0,
                     detect_ws=3,
                     merge_distance=None,
                     screen_hmin=None,
                     quiet=True,
                     **delineate_kwargs):
    """
    Run the whole pipeline in one call.

    Parameters
    ----------
    chm : ndarray or str or Path
        A CHM array, or a path to a raster.
    smooth_ws, smooth_method
        Passed to :func:`pyhrg.chm.smooth_chm`.
    hmin, detect_ws
        Passed to :func:`pyhrg.treetops.detect_tops`.
    merge_distance : float, optional
        If given, merge tops within this radius.
    screen_hmin : float, optional
        If given, drop tops below this height.
    quiet : bool
    **delineate_kwargs
        Passed to :meth:`CrownDelineator.delineate`.

    Returns
    -------
    crowns : ndarray[int32]
    tops : ndarray, shape (n, 2)

    Examples
    --------
    >>> crowns, tops = delineate_crowns("chm.tif", hmin=7)   # doctest: +SKIP
    """
    if isinstance(chm, np.ndarray):
        cd = CrownDelineator(chm, quiet=quiet)
    else:
        cd = CrownDelineator.from_file(chm, quiet=quiet)

    cd.smooth(ws=smooth_ws, method=smooth_method)
    cd.detect(hmin=hmin, ws=detect_ws)
    if merge_distance is not None:
        cd.merge(distance=merge_distance)
    if screen_hmin is not None:
        cd.screen(hmin=screen_hmin)
    crowns = cd.delineate(**delineate_kwargs)
    return crowns, cd.tops
