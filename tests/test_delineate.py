"""Tests for pyhrg.delineate and pyhrg.io."""
import numpy as np
import pytest

from pyhrg.delineate import CrownDelineator, delineate_crowns

rasterio = pytest.importorskip("rasterio", reason="file I/O needs rasterio")
from rasterio.transform import from_origin  # noqa: E402


@pytest.fixture
def chm_file(tmp_path, stand):
    chm, tops = stand
    p = tmp_path / "chm.tif"
    with rasterio.open(p, "w", driver="GTiff", height=chm.shape[0],
                       width=chm.shape[1], count=1, dtype="float32",
                       crs="EPSG:2180",
                       transform=from_origin(500000, 600000, 0.5, 0.5)) as dst:
        dst.write(chm, 1)
    return p, chm, tops


class TestPipeline:

    def test_full_chain(self, stand):
        chm, _ = stand
        cd = CrownDelineator(chm, quiet=True)
        cd.smooth(ws=3).detect(hmin=5, ws=5).merge(4.0).screen(8.0)
        crowns = cd.delineate(variance_thresh=6.0, mask_thresh=1.0)
        assert crowns.shape == chm.shape
        assert crowns.dtype == np.int32

    def test_steps_chain(self, stand):
        chm, _ = stand
        cd = CrownDelineator(chm, quiet=True)
        assert cd.smooth().detect() is cd

    def test_delineate_fills_in_missing_steps(self, stand):
        chm, _ = stand
        cd = CrownDelineator(chm, quiet=True)
        crowns = cd.delineate(variance_thresh=6.0, mask_thresh=1.0)
        assert cd.smoothed is not None
        assert cd.tops is not None
        assert crowns.shape == chm.shape

    def test_merge_before_detect_raises(self, stand):
        chm, _ = stand
        with pytest.raises(ValueError, match="detect"):
            CrownDelineator(chm, quiet=True).merge()

    def test_rejects_3d(self):
        with pytest.raises(ValueError, match="2-D"):
            CrownDelineator(np.zeros((2, 10, 10)))

    def test_export_without_geo_raises(self, stand):
        chm, _ = stand
        cd = CrownDelineator(chm, quiet=True)
        cd.delineate(variance_thresh=6.0, mask_thresh=1.0)
        with pytest.raises(ValueError, match="geotransform"):
            cd.to_raster("/tmp/nope.tif")

    def test_export_before_delineate_raises(self, stand):
        chm, _ = stand
        with pytest.raises(ValueError, match="delineate"):
            CrownDelineator(chm, quiet=True).to_raster("/tmp/nope.tif")

    def test_flat_chm_yields_no_crowns(self, flat):
        cd = CrownDelineator(flat, quiet=True)
        cd.smooth(ws=3).detect(hmin=20, ws=5)
        assert len(cd.tops) == 0
        crowns = cd.delineate(mask_thresh=0.5)
        assert crowns.max() == 0

    def test_quiet_is_silent(self, stand, capsys):
        chm, _ = stand
        cd = CrownDelineator(chm, quiet=True)
        cd.smooth().detect(hmin=5, ws=5)
        cd.delineate(variance_thresh=6.0, mask_thresh=1.0)
        assert capsys.readouterr().out == ""

    def test_verbose_reports(self, stand, capsys):
        chm, _ = stand
        CrownDelineator(chm, quiet=False).smooth()
        assert "smooth" in capsys.readouterr().out


class TestOneShot:

    def test_from_array(self, stand):
        chm, _ = stand
        crowns, tops = delineate_crowns(chm, hmin=5, detect_ws=5,
                                        variance_thresh=6.0, mask_thresh=1.0)
        assert crowns.shape == chm.shape
        assert tops.shape[1] == 2

    def test_from_file(self, chm_file):
        path, chm, _ = chm_file
        crowns, tops = delineate_crowns(str(path), hmin=5, detect_ws=5,
                                        variance_thresh=6.0, mask_thresh=1.0)
        assert crowns.shape == chm.shape

    def test_optional_steps(self, stand):
        chm, _ = stand
        crowns, tops = delineate_crowns(chm, hmin=5, detect_ws=5,
                                        merge_distance=4.0, screen_hmin=8.0,
                                        variance_thresh=6.0, mask_thresh=1.0)
        assert len(tops) > 0


class TestIO:

    def test_read_chm(self, chm_file):
        from pyhrg.io import read_chm
        path, chm, _ = chm_file
        arr, transform, crs = read_chm(path)
        assert arr.shape == chm.shape
        assert arr.dtype == np.float32
        assert crs.to_epsg() == 2180

    def test_window_shifts_transform(self, chm_file):
        from pyhrg.io import read_chm
        path, _, _ = chm_file
        full, t_full, _ = read_chm(path)
        win, t_win, _ = read_chm(path, window=(10, 20, 30, 40))
        assert win.shape == (40, 30)
        np.testing.assert_array_equal(win, full[20:60, 10:40])
        # origin must move by the window offset, or exports land in the sea
        assert t_win.c == pytest.approx(t_full.c + 10 * t_full.a)
        assert t_win.f == pytest.approx(t_full.f + 20 * t_full.e)

    def test_roundtrip_raster(self, chm_file, tmp_path):
        path, _, _ = chm_file
        cd = CrownDelineator.from_file(path, quiet=True)
        cd.smooth(ws=3).detect(hmin=5, ws=5)
        crowns = cd.delineate(variance_thresh=6.0, mask_thresh=1.0)
        out = tmp_path / "crowns.tif"
        cd.to_raster(out)
        with rasterio.open(out) as src:
            back = src.read(1)
            assert src.crs.to_epsg() == 2180
        np.testing.assert_array_equal(back, crowns)

    def test_from_file_with_window(self, chm_file):
        path, _, _ = chm_file
        cd = CrownDelineator.from_file(path, window=(10, 10, 50, 50), quiet=True)
        assert cd.chm.shape == (50, 50)


class TestCLI:

    def test_help(self):
        from pyhrg.cli import build_parser
        assert build_parser().prog == "pyhrg"

    def test_runs(self, chm_file, tmp_path):
        from pyhrg.cli import main
        path, _, _ = chm_file
        out = tmp_path / "cli.tif"
        rc = main(["-i", str(path), "-o", str(out), "--hmin", "5",
                   "--detect-ws", "5", "--variance-thresh", "6",
                   "--mask-thresh", "1", "-q"])
        assert rc == 0
        assert out.exists()

    def test_missing_input_reports_error(self, tmp_path):
        from pyhrg.cli import main
        rc = main(["-i", str(tmp_path / "nope.tif"),
                   "-o", str(tmp_path / "o.tif"), "-q"])
        assert rc == 1


class TestCrownCountReporting:
    """The count printed by delineate() must not assume a background exists.

    It was `len(np.unique(crowns)) - 1`, which subtracts one for label 0.
    With mask_thresh below the CHM's minimum every pixel is canopy, no 0
    appears, and two crowns were announced as one. The returned array was
    always right; only the number a user reads was wrong.
    """

    @staticmethod
    def _two_trees(base, gap):
        """`base` lifts the whole scene; `gap` cuts a strip of true zero.

        A Gaussian never reaches zero -- exp(-x) is tiny but positive -- so
        a scene built from Gaussians alone has no pixel at or below
        mask_thresh and therefore no background at all. The strip is what
        makes the with-background case actually have one.
        """
        n = 40
        yy, xx = np.mgrid[0:n, 0:n]
        g = lambda r, c, h, s: h * np.exp(
            -((yy - r) ** 2 + (xx - c) ** 2) / (2 * s ** 2))
        chm = base + np.maximum(g(12, 12, 15, 5), g(28, 28, 14, 5))
        if gap:
            chm[19:22, :] = 0.0
        return chm

    @pytest.mark.parametrize("base,gap,has_background",
                             [(0.0, True, True), (10.0, False, False)])
    def test_reported_count_matches_the_array(self, base, gap, has_background,
                                              capsys):
        from pyhrg import CrownDelineator
        chm = self._two_trees(base, gap)
        d = CrownDelineator(chm, quiet=False)
        crowns = d.smooth(ws=3).detect(hmin=5, ws=5).delineate(
            variance_thresh=8.0, mask_thresh=0.0)

        assert (0 in np.unique(crowns)) is has_background, \
            "the fixture no longer exercises the case it was built for"

        actual = len(set(np.unique(crowns)) - {0})
        line = [l for l in capsys.readouterr().out.splitlines()
                if "crowns" in l][-1]
        assert f"{actual} crowns" in line, \
            f"reported line {line!r} disagrees with {actual} crowns in the array"


class TestTreeTopExportPosition:
    """Exported points must land on the pixel centre, not its corner.

    An affine transform maps grid coordinates whose whole numbers are pixel
    *corners*, while a tree top from center_of_mass is an array index, and
    array indices refer to pixel centres. Without the half-pixel shift every
    exported point sat up and to the left -- 0.25 m on the 0.5 m test
    rasters, systematic, and enough to matter against field-measured stems.
    """

    def test_points_match_rasterio_xy(self, tmp_path, chm_file):
        rasterio = pytest.importorskip("rasterio")
        pytest.importorskip("fiona")
        import json
        from rasterio.transform import xy
        from pyhrg.io import save_tree_tops

        path, _, _ = chm_file
        with rasterio.open(path) as src:
            transform, crs = src.transform, src.crs
            chm = np.nan_to_num(src.read(1), nan=0.0)

        # One whole-pixel top and one subpixel, since center_of_mass
        # produces both and only the second catches a rounding fix.
        tops = np.array([[3.0, 4.0], [5.25, 6.75]])
        save_tree_tops(tops, str(tmp_path), "tt", transform, crs.to_wkt(),
                       chm, driver="GeoJSON")

        g = json.load(open(tmp_path / "tt_treetops.geojson"))
        for feature, (r, c) in zip(g["features"], tops):
            got = feature["geometry"]["coordinates"]
            expected = xy(transform, r, c)
            assert np.allclose(got, expected), \
                f"top ({r}, {c}) written at {got}, pixel centre is {expected}"
