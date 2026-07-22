"""Tests for pyhrg.chm and pyhrg.treetops."""
import numpy as np
import pytest

from pyhrg.chm import smooth_chm, SMOOTHING_METHODS
from pyhrg.treetops import detect_tops, merge_tops, screen_tops, as_pixels


class TestSmoothing:

    def test_all_methods_run(self, stand):
        chm, _ = stand
        for m in SMOOTHING_METHODS:
            out = smooth_chm(chm, ws=3, method=m)
            assert out.shape == chm.shape

    def test_median_reduces_spikes(self, stand):
        chm, _ = stand
        spiky = chm.copy()
        spiky[60, 60] = 200.0
        assert smooth_chm(spiky, ws=3, method="median")[60, 60] < 100.0

    def test_bad_method(self, stand):
        chm, _ = stand
        with pytest.raises(ValueError, match="method must be"):
            smooth_chm(chm, method="nope")

    def test_bad_ws(self, stand):
        chm, _ = stand
        with pytest.raises(ValueError, match="ws must be"):
            smooth_chm(chm, ws=0)

    def test_rejects_3d(self):
        with pytest.raises(ValueError, match="2-D"):
            smooth_chm(np.zeros((3, 10, 10)))


class TestDetectTops:

    def test_finds_the_trees(self, stand):
        chm, expected = stand
        tops = detect_tops(smooth_chm(chm, ws=3), hmin=5, ws=5)
        assert len(tops) == len(expected)

    def test_shape_and_dtype(self, stand):
        chm, _ = stand
        tops = detect_tops(chm, hmin=5, ws=5)
        assert isinstance(tops, np.ndarray)
        assert tops.ndim == 2 and tops.shape[1] == 2
        assert tops.dtype == np.float64

    def test_empty_is_still_2d(self, flat):
        """Regression: an empty (0,) array breaks any [:, 0] downstream."""
        tops = detect_tops(flat, hmin=20, ws=5)
        assert tops.shape == (0, 2)
        tops[:, 0]  # must not raise

    def test_hmin_filters(self, stand):
        chm, _ = stand
        assert len(detect_tops(chm, hmin=100, ws=5)) == 0

    def test_larger_window_detects_no_more(self, stand):
        chm, _ = stand
        sm = smooth_chm(chm, ws=3)
        assert len(detect_tops(sm, hmin=5, ws=9)) <= len(detect_tops(sm, hmin=5, ws=3))


class TestMergeTops:

    def test_merges_close_pair(self):
        tops = np.array([[10.0, 10.0], [10.0, 12.0], [50.0, 50.0]])
        merged = merge_tops(tops, distance=5.0)
        assert len(merged) == 2

    def test_keeps_distant(self):
        tops = np.array([[10.0, 10.0], [50.0, 50.0]])
        assert len(merge_tops(tops, distance=5.0)) == 2

    def test_merged_is_centroid(self):
        tops = np.array([[10.0, 10.0], [10.0, 20.0]])
        merged = merge_tops(tops, distance=20.0)
        assert len(merged) == 1
        np.testing.assert_allclose(merged[0], [10.0, 15.0])

    def test_empty_and_single(self):
        assert merge_tops(np.empty((0, 2))).shape == (0, 2)
        assert merge_tops(np.array([[1.0, 2.0]])).shape == (1, 2)

    def test_accepts_list(self):
        assert merge_tops([(1.0, 2.0), (1.0, 3.0)], distance=5.0).shape == (1, 2)


class TestScreenTops:

    def test_drops_short(self, stand):
        chm, _ = stand
        sm = smooth_chm(chm, ws=3)
        tops = detect_tops(sm, hmin=5, ws=5)
        kept = screen_tops(sm, tops, hmin=100.0)
        assert len(kept) == 0
        assert kept.shape == (0, 2)

    def test_keeps_tall(self, stand):
        chm, _ = stand
        sm = smooth_chm(chm, ws=3)
        tops = detect_tops(sm, hmin=5, ws=5)
        assert len(screen_tops(sm, tops, hmin=0.0)) == len(tops)

    def test_out_of_bounds_top_is_dropped_not_wrapped(self, stand):
        """Negative indices would silently wrap around the raster."""
        chm, _ = stand
        kept = screen_tops(chm, np.array([[-5.0, -5.0], [500.0, 500.0]]), hmin=0.0)
        assert len(kept) == 0


class TestAsPixels:

    def test_floors_to_int(self):
        px = as_pixels(np.array([[3.7, 4.2]]))
        assert px == [(3, 4)]
        assert all(isinstance(v, int) for v in px[0])

    def test_empty(self):
        assert as_pixels(np.empty((0, 2))) == []


class TestEvenWindowRejected:
    """An even window has no centre pixel and sits half a pixel off.

    Measured before the guard went in: smoothing a 40x55 scene and its
    mirror image differed by up to 8.8 m at ws=4, and detect_tops found
    397 tops against 400 on the mirror of chm_150_2023.tif at ws=4, and
    188 against 206 at ws=6. Nobody chooses ws=4 intending a shifted
    window, and nothing in the output shows that it happened.
    """

    @staticmethod
    def _scene():
        rng = np.random.default_rng(4)
        return rng.random((40, 55)) * 20 + 5

    @pytest.mark.parametrize("method", ["median", "mean", "maximum"])
    @pytest.mark.parametrize("ws", [2, 4, 6])
    def test_smooth_rejects_even_ws(self, method, ws):
        with pytest.raises(ValueError, match="odd"):
            smooth_chm(self._scene(), ws=ws, method=method)

    def test_smooth_allows_even_ws_for_gaussian(self):
        """ws only scales sigma there, so the kernel stays symmetric."""
        sc = self._scene()
        a = smooth_chm(sc, ws=4, method="gaussian")
        b = smooth_chm(sc[:, ::-1], ws=4, method="gaussian")[:, ::-1]
        np.testing.assert_allclose(a, b)

    @pytest.mark.parametrize("ws", [2, 4, 6])
    def test_detect_rejects_even_ws(self, ws):
        with pytest.raises(ValueError, match="odd"):
            detect_tops(self._scene(), hmin=5, ws=ws)

    @pytest.mark.parametrize("method", ["median", "mean", "maximum", "gaussian"])
    def test_odd_ws_is_orientation_independent(self, method):
        sc = self._scene()
        a = smooth_chm(sc, ws=5, method=method)
        b = smooth_chm(sc[:, ::-1], ws=5, method=method)[:, ::-1]
        np.testing.assert_allclose(a, b)

    def test_detection_with_odd_ws_is_orientation_independent(self):
        sc = self._scene()
        assert len(detect_tops(sc, hmin=5, ws=5)) == \
            len(detect_tops(sc[:, ::-1], hmin=5, ws=5))
