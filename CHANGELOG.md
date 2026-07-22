# Changelog

## [0.5.1] — 2026-07-22

### Fixed
- **0.5.0 did not import at all on Python 3.9.** `max_iters: int | None` is
  PEP 604, which before 3.10 is a runtime expression: `TypeError: unsupported
  operand type(s) for |: 'type' and 'NoneType'`, raised while the class body
  executes, so `import pyhrg.hrg` failed outright. `pyproject.toml` claims
  `requires-python = ">=3.9"`. Now `Optional[int]`.

  It reached a tag because development runs on 3.12 and every local test
  passed. `tests/test_python39_compat.py` now checks the package the way 3.9
  would: no PEP 604 in annotations outside a module with
  `from __future__ import annotations`, and every file `compile()`s — the
  first attempt at this fix put the future import after three `__author__`
  assignments, where `ast.parse` accepts it and Python does not.


## [0.5.0] — 2026-07-22

### Fixed
- **NaN was handed straight to `scipy.ndimage`, which is undefined on it
  rather than NaN-aware.** With the nine values 1..9 and a single NaN,
  `median_filter` returned 9, 5 or 4 depending on *where* in the window the
  NaN sat, and `maximum_filter` returned nan, 9 or 8. A NaN-skipping median
  would have been a constant 5.5. The comparison-based selection is poisoned
  by NaN and picks an arbitrary element.

  On `chm_150_2023.tif` that left **1208 pixels — 3 % of the raster —
  smoothed to an arbitrary window value**, with **20 of 253 tree tops
  standing on them** and errors against a NaN-skipping median of up to
  **26.29 m**. Nothing raised, nothing warned; the output looked like a
  canopy.

  `smooth_chm()` and `detect_tops()` now skip NaN: window statistics are
  taken over the cells that exist, and a window holding nothing but NaN
  stays NaN. The fast path still runs whenever the CHM has no NaN, so a
  clean raster costs nothing and gives byte-identical results to before.

  Found while fixing rHRG, which crashed on the same input. Neither was
  right: R failed loudly with an unhelpful message, Python produced
  plausible garbage. Both now agree, and the shared cross-check gained ten
  scenes with nodata holes so the path is actually covered.


## [0.4.0] — 2026-07-22

### Fixed
- **Exported tree tops sat half a pixel up and to the left.**
  `save_tree_tops()` mapped a top through `transform * (col, row)`, but an
  affine transform treats whole numbers as pixel *corners* while a top from
  `center_of_mass` is an array index, and array indices refer to pixel
  centres. Every exported point was therefore offset by half a pixel — 0.25 m
  on the 0.5 m test rasters, systematic, and enough to matter when comparing
  against field-measured stems. Now matches `rasterio.transform.xy`, verified
  on a whole-pixel and a subpixel top.

  Polygons were never affected: `save_segments()` goes through
  `rasterio.features.shapes`, which handles the convention itself. So points
  and polygons disagreed with each other.
- **`delineate()` under-reported the crown count by one when the scene had no
  background.** It computed `len(np.unique(crowns)) - 1`, subtracting for a
  label 0 that is absent whenever `mask_thresh` sits below the CHM's minimum.
  Two crowns were announced as one. The returned array was always correct;
  only the printed number was wrong.


## [0.3.0] — 2026-07-22

### Changed
- **`max_iters` now defaults to `None`, meaning natural termination.** It was
  200, and the bound bit: on `chm_33_2012.tif` 332 of 492 crowns stopped there
  with candidates still queued, and the crown count read **132 against 63**
  once the cap was lifted. More than a factor of two in the headline number,
  decided by a constant rather than by the canopy.

  The boundaries barely moved — 2.9 % of the partition — because a truncated
  grow blocks merges rather than misplacing pixels. That is what made it hard
  to notice: the segmentation looked right and the tree count did not.

  The cap protected nothing. Growth is bounded anyway, since every iteration
  either accepts a region or records a rejection and there are finitely many
  of both, and natural termination needed at most 484 iterations on that scene
  while running *faster* — 0.67 s against 2.34 s, because twice as many
  surviving crowns cost more in conflict arbitration than the extra merges
  cost in growing.

  `delineate_crowns()` had its own default of 200 and now passes `None` too.
  An explicit `max_iters` still works and warns when it binds; the message
  carries no seed id, so Python's default filter collapses what would
  otherwise have been 332 warnings into one.

## [0.2.0] — 2026-07-22

### Changed
- **An even smoothing or detection window is now refused.** An even window
  has no centre pixel, so `scipy.ndimage` places it half a pixel off and the
  result depends on which way the raster happens to be oriented. Measured:
  smoothing a 40x55 scene and its mirror image differed by up to 8.8 m at
  `ws=4`; `detect_tops` on `chm_150_2023.tif` found 397 tops and 400 on the
  mirror at `ws=4`, and 188 against 206 at `ws=6` — a 9 per cent difference
  in the tree count from orientation alone.

  `ws` was validated only as `>= 1`, and nothing in the output shows that it
  happened, so this raises rather than warns. `method="gaussian"` is exempt:
  `ws` only scales sigma there and the kernel stays symmetric.

  The same guard went into rHRG, which reproduced the asymmetry faithfully.


All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.0] — 2026-07-15

First release. pyHRG supersedes `pycrown_simplified`, which is archived.

### Added
- `CrownDelineator` — chainable pipeline: smooth → detect → merge → screen → delineate.
- `delineate_crowns()` — the whole pipeline in one call.
- Array-first API: `smooth_chm`, `detect_tops`, `merge_tops`, `screen_tops`,
  `as_pixels`, `HierarchicalRegionGrower`, `resolve_conflicts`.
- `conflict_rule` on `run_all()` — `'height'`, `'distance'` or `'similarity'`
  decide which crown wins canopy claimed by two. Ties break to the lower
  crown id, so the segmentation is reproducible whatever order the tree
  tops arrive in.
- `protect_seeds` — guarantee one crown per tree top, disabling merging.
- `retry_rejected` — reconsider regions rejected earlier in a grow.
- `n_contested` — how many regions were claimed by more than one crown.
- `read_chm()` with windowed reads, for rasters that do not fit in memory.
- Command line: `pyhrg -i chm.tif -o crowns.tif`, and `python -m pyhrg`.
- Test suite (65 tests) and CI.

### Changed
- **Split from PyCrown's structure.** Hierarchical region growing is now the
  only delineation method (see *Removed*). Smoothing, tree top detection and
  the algorithm each live in their own module, and none of them touch the disk
  — file access is confined to `pyhrg.io`.
- **Contested regions are arbitrated, not overwritten.** Previously the crown
  with the higher id silently won, so which of two merging trees survived
  depended on tree-top ordering.
- Statistics merging is attributed correctly to the pairwise formula of
  Chan, Golub & LeVeque (1983). The code called it Welford's method, which is
  a different algorithm — it adds one sample at a time, whereas a crown
  absorbs a whole region at once.
- Parallel growing shares the graph per worker instead of per task. The
  CSR graph was previously pickled once per seed: 122 MB of IPC for 1156
  trees, versus 0.4 MB now. Results are identical to sequential.
- `tree_detection` (now `detect_tops`) returns `(n, 2)` float64 rather than a
  list of tuples, matching every other tree-top function.

### Fixed
- **Growing crashed on any CHM without trees.** An empty detection produced a
  `(0,)` array instead of `(0, 2)`, so delineation raised
  `IndexError: too many indices for array`. Tiles covering clearcuts, gaps or
  water hit this.
- **Growing was unusable on in-memory arrays.** The algorithm demanded a file
  path, so `PyCrown(chm_array=...)` raised
  `TypeError: invalid path or file: None` on `hierarchical_crown_delineation`.
- Tree top correction returned early with fewer than two tops without storing
  its result, leaving the object's state inconsistent with every other path.
- Tree tops outside the raster are rejected explicitly; negative coordinates
  previously wrapped around silently.
- Dead code removed: an unused inverse-statistics helper, an unused loop
  variable in the merge formula, and a branch that could never be taken.

### Packaging
- Licence declared as an SPDX expression with `license-files` (PEP 639)
  rather than a TOML table, which setuptools deprecated and will stop
  accepting on 2027-02-18. The redundant licence classifier is gone.
  Requires setuptools >= 77 to build.
- Packages resolved with an explicit `find` directive rather than a
  hard-coded list, so `tests/` and `examples/` are never shipped and any
  future subpackage is picked up instead of being silently dropped.
  Reported against `pycrown_simplified` by Bas van Driel.

### Removed
- **Dalponte & Coomes delineation.** It is well covered by
  [PyCrown](https://github.com/manaakiwhenua/pycrown),
  [lidR](https://github.com/r-lidar/lidR) and
  [itcSegment](https://cran.r-project.org/package=itcSegment); duplicating it
  here served no one. pyHRG does one thing.
- Point cloud I/O. pyHRG takes a CHM raster and returns crowns.
