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

## Releasing

The checklist exists because of a specific failure. `max_iters` changed
default in 0.3.0, and that one change broke CI in two packages at once —
pyHRG with `int | None`, which is a runtime `TypeError` before Python 3.10
while the metadata claims `>=3.9`, and rHRG with a stale `man/` page. Neither
was noticed. **pyHRG then tagged 0.3.0, 0.4.0 and 0.5.0 with the workflow
red**, so three releases could not be imported on the minimum Python they
advertise. rHRG shipped two the same way, rgeoadaptels two more.

Local tests passed in every one of those cases. They were run on one
interpreter, on one operating system, by someone who already knew what the
change was meant to do. The matrix is the part that disagrees.

1. Update `CHANGELOG.md`. If the output changes, say so in those words.
2. Bump the version everywhere it appears. Search for the *old* number and
   read the hits — `grep -rn "0.4.0" --exclude-dir=.git` — rather than
   editing the two or three places you remember.
3. Run the tests locally.
   `pytest tests/` includes `test_python39_compat.py`, which rejects
   syntax the oldest supported Python cannot parse. Development runs on
   3.12, where the defect above was invisible.
4. Commit and push. **Do not tag yet.**
5. **Wait for Actions on the pushed commit and confirm every matrix job is
   green.** Not the previous run, not the branch generally — that commit.
   This is the step that was missing. Either open the Actions tab, or:

   ```bash
   curl -s "https://api.github.com/repos/OWNER/REPO/actions/runs?per_page=1" |
     python -c "import json,sys; r=json.load(sys.stdin)['workflow_runs'][0]; print(r['head_sha'][:7], r['status'], r['conclusion'])"
   ```

   `gh run list` is nicer if the GitHub CLI is installed; it is not
   everywhere, and the curl form needs nothing but a public repo.
6. Only then tag and push the tag:
   `git tag -a vX.Y.Z -m "..." && git push --tags`

The order matters. A tag is what people install and what a DOI points at, so
it should never be the thing that discovers a broken build. If Actions is
red, fix it and release the fix as its own version — the broken tag stays in
history either way.

## Scope

pyHRG delineates crowns from a CHM by hierarchical region growing. It is
deliberately narrow. Other delineation algorithms, point cloud processing and
CHM generation are out of scope — PyCrown, lidR and itcSegment already cover
that ground well.

## Licence

pyHRG is GPLv3, inherited from PyCrown. Contributions are accepted under the
same licence.
