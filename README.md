# PyHRG

<img src="https://raw.githubusercontent.com/igorpawelec/pyhrg/main/www/logopy.png" alt="pyHRG logo" align="right" width="200"/>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org)

**Individual tree crown delineation from canopy height models, by Hierarchical Region Growing.**

Pure Python + Numba. No compiled extensions, no external binaries.

> **R users:** an R implementation of the same algorithm lives in [rHRG](https://github.com/igorpawelec/rhrg). The two are separate packages by design — installation, tooling and idioms differ too much to share a repository — but they implement the same method and are validated against each other — exactly on the shared synthetic suite, and to within 0.25 % of watershed pixels on real canopy height models, where the two break plateau ties differently.

## Background

Delineating individual crowns from a canopy height model runs into one persistent problem: **tree tops are over-detected**. A single broad crown has a ragged upper surface, so local-maximum detection finds several peaks on it. Lower the sensitivity and you start losing real trees instead.

pyHRG treats that as the central problem rather than a preprocessing nuisance. Surplus tree tops are allowed, and the growing merges them back:

1. **Watershed.** Tree tops seed a marker-based watershed on the inverted CHM. This yields *exactly one region per tree top* — so the regions **are** the detected trees.
2. **Region adjacency graph.** Neighbouring regions get a weighted edge,
   *w(a,b) = α·|Δμ| + β·|Δσ| + γ/(border+1)*, stored as CSR arrays. Region statistics come from the same single pixel pass.
3. **Growing.** Each tree grows greedily, absorbing neighbours while the combined height variance stays under a threshold. When two trees absorb each other, they are one crown — this is how over-detection is corrected.
4. **Arbitration.** Growing runs independently per seed, so two crowns can claim the same region. Those claims are settled on the data (taller tree, nearest tree, or best height match) rather than on iteration order, which makes the segmentation reproducible regardless of how the tree tops were ordered.

`variance_thresh` is the main control: it sets how much height variation a single crown may contain, and therefore how readily neighbouring trees merge.

### Relation to PyCrown

pyHRG began as a fork of [PyCrown](https://github.com/manaakiwhenua/pycrown) (Zörner et al. 2018) and keeps its pipeline shape — smooth the CHM, find tree tops as local maxima, delineate crowns. **The Dalponte & Coomes delineation that PyCrown re-implements is not part of pyHRG**; hierarchical region growing is the only method here. If you want Dalponte, use PyCrown, [lidR](https://github.com/r-lidar/lidR) or [itcSegment](https://cran.r-project.org/package=itcSegment) — they do it well and there is no reason to duplicate them.

Differences from that lineage worth knowing:

| | PyCrown | pyHRG |
|---|---|---|
| Delineation | Dalponte & Coomes | Hierarchical Region Growing |
| Over-detected tops | filtered out beforehand | merged by the growing |
| Contested canopy | n/a | arbitrated explicitly, reproducibly |
| Point cloud I/O | yes (laspy) | no — CHM raster in, crowns out |

## Installation

Native dependencies come from conda; pip then installs the package without touching them.

```bash
conda install -c conda-forge numpy numba scipy scikit-image rasterio fiona
pip install --no-deps git+https://github.com/igorpawelec/pyhrg.git
```

The algorithm itself needs only numpy, numba, scipy and scikit-image. `rasterio` and `fiona` are used solely for reading and writing files — the array API works without them.

## Quick start

```python
from pyhrg import CrownDelineator

cd = CrownDelineator.from_file("chm.tif")
cd.smooth(ws=3).detect(hmin=7, ws=5).merge(5.0).screen(10.0)
crowns = cd.delineate(variance_thresh=2.0)

cd.to_raster("crowns.tif")
cd.to_vector("out/", name="crowns")
```

One call, if you do not need the intermediate state:

```python
from pyhrg import delineate_crowns

crowns, tops = delineate_crowns("chm.tif", hmin=7, merge_distance=5.0,
                                variance_thresh=2.0)
```

### Arrays, without touching the disk

Every stage is a plain function. Nothing in the algorithm needs a file path:

```python
import numpy as np
from pyhrg import smooth_chm, detect_tops, as_pixels, HierarchicalRegionGrower

chm = np.load("chm.npy")
smoothed = smooth_chm(chm, ws=3, method="median")
tops = detect_tops(smoothed, hmin=7, ws=5)          # (n, 2) subpixel row/col

grower = HierarchicalRegionGrower(smoothed)
crowns = grower.run_all(as_pixels(tops), variance_thresh=2.0)
print(f"{grower.n_contested} regions were claimed by more than one crown")
```

### Large rasters

Read a window; the geotransform is shifted to match, so exported crowns stay georeferenced.

```python
# 2000 x 2000 px starting at col=1000, row=500
cd = CrownDelineator.from_file("big_chm.tif", window=(1000, 500, 2000, 2000))
```

### Command line

```bash
pyhrg -i chm.tif -o crowns.tif --hmin 7 --variance-thresh 2.0
pyhrg -i chm.tif -o crowns.tif --vector out/ --merge-distance 5 --screen-hmin 10
python -m pyhrg --help
```

## Parameters

**Smoothing** — `smooth_chm(chm, ws, method)`

| Parameter | Default | Description |
|---|---|---|
| `ws` | 3 | Window size (px). Larger = fewer false tops, but merges close crowns |
| `method` | `'median'` | `median`, `mean`, `gaussian`, `maximum`. Median keeps crown edges sharp |

**Tree tops** — `detect_tops`, `merge_tops`, `screen_tops`

| Parameter | Default | Description |
|---|---|---|
| `hmin` | 2.0 | Minimum height (m) for a pixel to be a candidate |
| `ws` | 3 | Local-maximum window (px) = minimum spacing between tops |
| `distance` | 5.0 | Merge radius (px). Grouping is transitive — keep below crown diameter |

**Growing** — `delineate` / `run_all`

| Parameter | Default | Description |
|---|---|---|
| `variance_thresh` | 2.0 | Max height variance (σ²) within a crown. **The main lever** |
| `mask_thresh` | 0.0 | Minimum CHM height treated as canopy (m) |
| `morpho_radius` | 0 | Disk radius for opening/closing the mask. 0 = off |
| `alpha`, `beta`, `gamma` | 1.0, 0.5, 0.1 | Edge weights: mean diff, σ diff, inverse border length |
| `anneal_lambda` | 1.0 | Per-iteration tightening of the threshold. 1.0 = constant |
| `max_iters` | `None` | Cap on grow iterations per seed. `None` grows to natural termination |
| `conflict_rule` | `'height'` | Who wins contested canopy — see below |
| `protect_seeds` | `False` | If True, no tree is ever absorbed; every top yields a crown |
| `retry_rejected` | `False` | Reconsider regions rejected earlier in the same grow |
| `n_jobs` | 1 | Parallel processes. -1 = all cores |

### Conflict rules

Two crowns can claim the same region. `conflict_rule` decides who gets it:

- **`'height'`** (default) — the taller tree wins. Dominant trees overtop their neighbours, so ambiguous canopy goes to the taller crown.
- **`'distance'`** — the nearest seed wins, by distance from the region centroid. Classic ITC behaviour; splits rather than merges.
- **`'similarity'`** — the tree whose apex height best matches the region's mean height wins.

Ties resolve to the lower crown id, so output is fully reproducible. The count of contested regions is available afterwards as `grower.n_contested`.

`protect_seeds=True` disables merging entirely — every input tree top keeps its own crown. Use it when the tops are trusted, e.g. field-measured.

## Notes on behaviour

**The crown count can be lower than the tree-top count.** That is the point: merged trees leave gaps in the id sequence. If you need one crown per top, use `protect_seeds=True`.

**`retry_rejected` only bites in the middle.** A rejected region can become admissible later, because absorbing a large homogeneous neighbour can *lower* a crown's variance. On a 1105-tree synthetic scene this changed nothing at `variance_thresh` 2 or 8 (rejections are final anyway) and nothing at 120 (everything merges regardless), but changed ~17% of the raster at 20. The direction is not predictable — the crown count moved both up and down.

**`n_jobs > 1` is worth it only on large scenes.** The graph is shared per worker rather than per task, but process start-up still costs ~0.1 s each, and growing is cheap per tree (~0.4 s for 1150 trees). Results are identical to sequential either way.

## Performance

Numba compiles the one place that touches every pixel: the single pass that
builds the region adjacency graph and the per-region statistics together.
On a 1000x1000 px CHM with 4900 trees that pass takes **9.7 ms**, against
128 ms for an equivalent written with `numpy.bincount` — a 13x gap, mostly
because it does in one pass what numpy needs several for.

The grow loop is deliberately *not* compiled. It is driven by a heap and set
membership, which Numba does not handle well, and it runs once per tree
rather than once per pixel — roughly 0.4 s for 1150 trees. It is not the
bottleneck.

First call in a session pays the JIT compilation cost (~2 s); `cache=True`
means later runs read the compiled code from disk.

## Testing

```bash
pip install -e ".[test]"
pytest tests/ -v
```

## Repository structure

```
pyhrg/
├── pyhrg/
│   ├── __init__.py      # public API (lazy imports)
│   ├── __main__.py      # python -m pyhrg
│   ├── chm.py           # CHM smoothing
│   ├── treetops.py      # detection, merging, screening
│   ├── hrg.py           # the algorithm: watershed, RAG, growing, arbitration
│   ├── delineate.py     # CrownDelineator pipeline
│   ├── io.py            # raster/vector I/O (rasterio, fiona)
│   └── cli.py           # command line
├── tests/
├── examples/
├── pyproject.toml
├── environment.yaml
├── CITATION.cff
├── CHANGELOG.md
├── CONTRIBUTING.md
└── LICENSE
```

## Requirements

- Python ≥ 3.9
- NumPy ≥ 1.21, Numba ≥ 0.56, SciPy ≥ 1.7, scikit-image ≥ 0.19
- Rasterio ≥ 1.3, Fiona ≥ 1.9 *(only for file I/O)*

## Citation

If you use pyHRG in your research, please cite the software and the work it builds on:

**This implementation**

> Pawelec, I. (2026). *pyHRG: individual tree crown delineation from canopy height models by Hierarchical Region Growing* [Software]. https://github.com/igorpawelec/pyhrg

**Upstream project.** pyHRG is a derivative of PyCrown and keeps its pipeline structure:

> Zörner, J., Dymond, J., Shepherd, J., Jolly, B. (2018). *PyCrown — Fast raster-based individual tree segmentation for LiDAR data.* Landcare Research NZ Ltd. https://doi.org/10.7931/M0SR-DN55
>
> Zörner, J., Dymond, J.R., Shepherd, J.D., Wiser, S.K., Bunting, P., Jolly, B. (2018). LiDAR-based regional inventory of tall trees — Wellington, New Zealand. *Forests* 9(11), 702. https://doi.org/10.3390/f9110702

**Methods used**

> Chan, T.F., Golub, G.H., LeVeque, R.J. (1983). Algorithms for computing the sample variance: analysis and recommendations. *The American Statistician* 37(3), 242–247. — the pairwise formula that merges two regions' statistics in O(1). Note this is *not* Welford's update, which adds one sample at a time; a crown absorbs a whole region at once.
>
> Beucher, S., Meyer, F. (1993). The morphological approach to segmentation: the watershed transformation. In: *Mathematical Morphology in Image Processing*, 433–481. — the marker-based watershed, via scikit-image.
>
> Popescu, S.C., Wynne, R.H. (2004). Seeing the trees in the forest. *Photogrammetric Engineering & Remote Sensing* 70(5), 589–604. — local-maxima tree detection on a smoothed CHM.

See [CITATION.cff](CITATION.cff) for machine-readable metadata.

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

pyHRG derives from PyCrown, which is published under GPLv3; the licence carries over.

## Contributing

Contributions welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).
