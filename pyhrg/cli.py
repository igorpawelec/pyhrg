"""
pyhrg.cli — command line interface.

    pyhrg -i chm.tif -o crowns.tif --hmin 7 --variance-thresh 2.0

Copyright (C) 2025 Igor Pawelec
Licence: GPLv3 — see LICENSE.
"""

import argparse
import sys

__all__ = ["build_parser", "main"]


def build_parser():
    p = argparse.ArgumentParser(
        prog="pyhrg",
        description="Delineate tree crowns from a canopy height model "
                    "using Hierarchical Region Growing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, metavar="CHM.tif",
                   help="input canopy height model raster")
    p.add_argument("-o", "--output", required=True, metavar="CROWNS.tif",
                   help="output crown raster (GeoTIFF)")
    p.add_argument("--vector", metavar="DIR",
                   help="also write crown polygons into this directory")
    p.add_argument("--window", nargs=4, type=int, default=None,
                   metavar=("COL", "ROW", "W", "H"),
                   help="process only this sub-region of the raster")

    g = p.add_argument_group("smoothing")
    g.add_argument("--smooth-ws", type=int, default=3,
                   help="smoothing window size (pixels)")
    g.add_argument("--smooth-method", default="median",
                   choices=("median", "mean", "gaussian", "maximum"))

    g = p.add_argument_group("tree tops")
    g.add_argument("--hmin", type=float, default=2.0,
                   help="minimum height for a tree top (m)")
    g.add_argument("--detect-ws", type=int, default=3,
                   help="local-maximum window size (pixels)")
    g.add_argument("--merge-distance", type=float, default=None,
                   help="merge tops closer than this (pixels)")
    g.add_argument("--screen-hmin", type=float, default=None,
                   help="drop tops below this height (m)")

    g = p.add_argument_group("growing")
    g.add_argument("--variance-thresh", type=float, default=2.0,
                   help="maximum height variance within a crown")
    g.add_argument("--mask-thresh", type=float, default=0.0,
                   help="minimum CHM height to consider as canopy (m)")
    g.add_argument("--morpho-radius", type=int, default=0,
                   help="disk radius for cleaning the mask (0 = off)")
    g.add_argument("--conflict-rule", default="height",
                   choices=("height", "distance", "similarity"),
                   help="who wins canopy claimed by two crowns")
    g.add_argument("--protect-seeds", action="store_true",
                   help="never let one tree absorb another")
    g.add_argument("--retry-rejected", action="store_true",
                   help="reconsider regions rejected earlier in a grow")
    g.add_argument("--n-jobs", type=int, default=1,
                   help="parallel processes (-1 = all cores)")

    p.add_argument("-q", "--quiet", action="store_true",
                   help="suppress progress messages")
    p.add_argument("-V", "--version", action="store_true",
                   help="print version and exit")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version:
        from pyhrg import __version__
        print(f"pyhrg {__version__}")
        return 0

    from pyhrg import CrownDelineator

    try:
        cd = CrownDelineator.from_file(args.input, window=args.window,
                                       quiet=args.quiet)
        cd.smooth(ws=args.smooth_ws, method=args.smooth_method)
        cd.detect(hmin=args.hmin, ws=args.detect_ws)
        if args.merge_distance is not None:
            cd.merge(distance=args.merge_distance)
        if args.screen_hmin is not None:
            cd.screen(hmin=args.screen_hmin)
        cd.delineate(
            variance_thresh=args.variance_thresh,
            mask_thresh=args.mask_thresh,
            morpho_radius=args.morpho_radius,
            conflict_rule=args.conflict_rule,
            protect_seeds=args.protect_seeds,
            retry_rejected=args.retry_rejected,
            n_jobs=args.n_jobs,
        )
        cd.to_raster(args.output)
        if args.vector:
            cd.to_vector(args.vector)
    except (OSError, ValueError) as e:
        print(f"pyhrg: error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
