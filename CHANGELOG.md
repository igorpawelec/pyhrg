# Changelog

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
