"""
pyHRG quickstart — runs on a synthetic stand, no data files needed.

    python examples/quickstart.py

Shows the three things worth understanding about the algorithm:
  1. the pipeline,
  2. how over-detected tree tops get merged back,
  3. what the conflict rules actually do.
"""
import numpy as np

from pyhrg import CrownDelineator, delineate_crowns
from pyhrg.hrg import CONFLICT_RULES


def synthetic_stand(n=140, seed=7):
    """A stand of Gaussian crowns, plus one broad crown with a double peak."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n]
    chm = np.zeros((n, n))

    def bump(cy, cx, h, s):
        return h * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * s ** 2))

    for _ in range(20):
        r, c = rng.integers(12, n - 12, 2)
        chm = np.maximum(chm, bump(r, c, rng.uniform(14, 28), rng.uniform(4, 6)))

    # one wide crown whose surface carries two local maxima
    wide = bump(70, 70, 24, 11) + bump(70, 64, 1.5, 3) + bump(70, 76, 2.0, 3)
    chm = np.maximum(chm, wide)
    return chm.astype(np.float32)


chm = synthetic_stand()

print("=" * 62)
print("1) The pipeline")
print("=" * 62)
cd = CrownDelineator(chm)
cd.smooth(ws=3, method="median")
cd.detect(hmin=5, ws=5)
cd.merge(distance=4.0)
cd.screen(hmin=8.0)
crowns = cd.delineate(variance_thresh=6.0, mask_thresh=1.0)

print()
print("=" * 62)
print("2) Same thing in one call")
print("=" * 62)
crowns2, tops2 = delineate_crowns(
    chm, hmin=5, detect_ws=5, merge_distance=4.0, screen_hmin=8.0,
    variance_thresh=6.0, mask_thresh=1.0,
)
print(f"  {len(np.unique(crowns2)) - 1} crowns from {len(tops2)} tree tops")

print()
print("=" * 62)
print("3) Conflict rules, on a deliberately loose threshold")
print("=" * 62)
print("  A high variance_thresh lets neighbouring trees merge, which is")
print("  when the rules start to matter.\n")
for rule in CONFLICT_RULES:
    cd = CrownDelineator(chm, quiet=True)
    cd.smooth(ws=3).detect(hmin=5, ws=5).merge(4.0)
    out = cd.delineate(variance_thresh=25.0, mask_thresh=1.0,
                       conflict_rule=rule)
    n = len(np.unique(out)) - 1
    print(f"  {rule:11s} -> {n:3d} crowns from {len(cd.tops):3d} tops "
          f"({cd.grower.n_contested} contested regions)")

cd = CrownDelineator(chm, quiet=True)
cd.smooth(ws=3).detect(hmin=5, ws=5).merge(4.0)
out = cd.delineate(variance_thresh=25.0, mask_thresh=1.0, protect_seeds=True)
print(f"  {'protected':11s} -> {len(np.unique(out)) - 1:3d} crowns from "
      f"{len(cd.tops):3d} tops   (merging disabled)")

print()
print("Reproducibility: rerun this file — the numbers will be identical,")
print("and they do not depend on the order the tree tops came in.")
