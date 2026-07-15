# Contributing to pyHRG

Bug reports, ideas and pull requests are all welcome.

## Reporting a bug

Please include:
- what you ran (ideally a snippet that reproduces it),
- what you expected and what happened,
- the raster's shape, dtype and value range — most surprises trace back to
  nodata, units or an unexpected CHM range,
- versions: `python -c "import pyhrg, numpy, numba; print(pyhrg.__version__, numpy.__version__, numba.__version__)"`.

A CHM that triggers the problem helps enormously. A synthetic one that
reproduces it helps even more.

## Development setup

```bash
conda env create -f environment.yaml
conda activate pyhrg
pip install --no-deps -e ".[dev]"
pytest tests/ -v
```

The build needs setuptools >= 77: the licence is declared as an SPDX
expression (PEP 639), which older versions do not understand. pip installs
its own build environment, so this only bites if you build with
`--no-build-isolation`.

## Pull requests

- Add a test that fails before your change and passes after.
- Keep `pyflakes pyhrg/*.py` clean.
- The algorithm modules (`chm`, `treetops`, `hrg`) work on arrays and must
  not import rasterio or fiona — file access belongs in `pyhrg.io`.
- If you change what the segmentation produces, say so in the PR and in
  `CHANGELOG.md`. Silent changes to output are the hardest kind to debug for
  anyone with a pipeline in flight.
- Numbers in docstrings should come from a measurement, not an estimate.

## Scope

pyHRG delineates crowns from a CHM by hierarchical region growing. It is
deliberately narrow. Other delineation algorithms, point cloud processing and
CHM generation are out of scope — PyCrown, lidR and itcSegment already cover
that ground well.

## Licence

pyHRG is GPLv3, inherited from PyCrown. Contributions are accepted under the
same licence.
