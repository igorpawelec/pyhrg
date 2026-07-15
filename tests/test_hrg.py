"""Tests for pyhrg.hrg — the algorithm itself."""
import numpy as np
import pytest

from pyhrg.hrg import (
    HierarchicalRegionGrower, resolve_conflicts, CONFLICT_RULES, _merge_stats,
)


class TestMergeStats:

    def test_matches_numpy_on_a_real_split(self):
        rng = np.random.default_rng(0)
        a = rng.normal(20, 3, 400)
        b = rng.normal(25, 2, 150)
        n, m, v = _merge_stats(len(a), a.mean(), a.var(),
                              len(b), b.mean(), b.var())
        both = np.concatenate([a, b])
        assert n == len(both)
        assert m == pytest.approx(both.mean())
        assert v == pytest.approx(both.var())

    def test_empty(self):
        assert _merge_stats(0, 0.0, 0.0, 0, 0.0, 0.0) == (0, 0.0, 0.0)

    def test_variance_can_drop_when_absorbing_a_homogeneous_region(self):
        """Why retry_rejected exists: rejection is not permanent in principle."""
        _, _, v_after = _merge_stats(100, 20.0, 6.0, 400, 20.0, 0.1)
        assert v_after < 6.0


class TestResolveConflicts:
    """Two crowns, each claiming both regions — the common mutual case."""

    @pytest.fixture
    def case(self):
        return dict(
            results={1: {1, 2}, 2: {1, 2}},
            seed_ids=[1, 2],
            seed_rc=[(0, 0), (0, 10)],
            seed_height=np.array([20.0, 30.0]),   # crown 2 is taller
            reg_mean=np.array([0.0, 20.0, 30.0]),
            reg_cy=np.array([0.0, 0.0, 0.0]),
            reg_cx=np.array([0.0, 0.0, 10.0]),
        )

    def test_height_gives_all_to_the_taller(self, case):
        a, n = resolve_conflicts(**case, rule="height")
        assert a[1] == 2 and a[2] == 2
        assert n == 2

    def test_distance_splits_by_nearest_seed(self, case):
        a, _ = resolve_conflicts(**case, rule="distance")
        assert a[1] == 1 and a[2] == 2

    def test_similarity_matches_region_height_to_seed(self, case):
        a, _ = resolve_conflicts(**case, rule="similarity")
        assert a[1] == 1 and a[2] == 2

    def test_protect_seeds_stops_absorption(self, case):
        a, _ = resolve_conflicts(**case, rule="height", protect_seeds=True)
        assert a[1] == 1 and a[2] == 2

    def test_bad_rule(self, case):
        with pytest.raises(ValueError, match="rule must be"):
            resolve_conflicts(**case, rule="bogus")

    def test_ties_break_to_lower_id(self, case):
        case["seed_height"] = np.array([25.0, 25.0])
        a, _ = resolve_conflicts(**case, rule="height")
        assert a[1] == 1 and a[2] == 1

    def test_deterministic_under_shuffling(self, case):
        ref = None
        for i in range(10):
            rng = np.random.default_rng(i)
            shuffled = dict(case)
            shuffled["results"] = {k: set(rng.permutation(list(v)))
                                   for k, v in case["results"].items()}
            a, _ = resolve_conflicts(**shuffled, rule="height")
            if ref is None:
                ref = a
            assert a == ref

    def test_uncontested_region_needs_no_rule(self, case):
        case["results"] = {1: {1}, 2: {2}}
        a, n = resolve_conflicts(**case, rule="height")
        assert a == {1: 1, 2: 2}
        assert n == 0


class TestGrower:

    def test_rejects_3d_input(self):
        with pytest.raises(ValueError, match="2-D"):
            HierarchicalRegionGrower(np.zeros((2, 10, 10)))

    def test_takes_an_array_not_a_path(self, stand):
        """Regression: the algorithm must not require a file on disk."""
        chm, tops = stand
        g = HierarchicalRegionGrower(chm)
        out = g.run_all(tops, variance_thresh=6.0, mask_thresh=1.0)
        assert out.shape == chm.shape

    def test_out_of_bounds_seed_rejected(self, stand):
        chm, tops = stand
        with pytest.raises(ValueError, match="outside the CHM"):
            HierarchicalRegionGrower(chm).run_all(tops + [(999, 999)])

    def test_bad_conflict_rule(self, stand):
        chm, tops = stand
        with pytest.raises(ValueError, match="conflict_rule must be"):
            HierarchicalRegionGrower(chm).run_all(tops, conflict_rule="bogus")

    def test_no_seeds(self, flat):
        out = HierarchicalRegionGrower(flat).run_all([], mask_thresh=0.5)
        assert out.shape == flat.shape
        assert out.max() == 0

    def test_every_rule_runs(self, stand):
        chm, tops = stand
        for rule in CONFLICT_RULES:
            g = HierarchicalRegionGrower(chm)
            out = g.run_all(tops, variance_thresh=6.0, mask_thresh=1.0,
                            conflict_rule=rule)
            assert out.max() <= len(tops)

    def test_reproducible(self, stand):
        chm, tops = stand
        a = HierarchicalRegionGrower(chm).run_all(tops, variance_thresh=6.0,
                                                  mask_thresh=1.0)
        b = HierarchicalRegionGrower(chm).run_all(tops, variance_thresh=6.0,
                                                  mask_thresh=1.0)
        np.testing.assert_array_equal(a, b)

    def test_seed_order_does_not_pick_the_winner(self, two_tops_one_tree):
        """The old bug: the crown with the higher id won, arbitrarily."""
        chm, tops = two_tops_one_tree
        kw = dict(variance_thresh=60.0, mask_thresh=1.0, conflict_rule="height")
        fwd = HierarchicalRegionGrower(chm).run_all(tops, **kw)
        rev = HierarchicalRegionGrower(chm).run_all(tops[::-1], **kw)
        # ids follow input order, so the partition must match with ids swapped
        np.testing.assert_array_equal(fwd == 2, rev == 1)

    def test_over_detected_tops_merge(self, two_tops_one_tree):
        chm, tops = two_tops_one_tree
        g = HierarchicalRegionGrower(chm)
        out = g.run_all(tops, variance_thresh=60.0, mask_thresh=1.0,
                        conflict_rule="height")
        assert len(np.unique(out)) - 1 == 1      # two tops, one crown
        assert g.n_contested > 0

    def test_protect_seeds_keeps_both(self, two_tops_one_tree):
        chm, tops = two_tops_one_tree
        out = HierarchicalRegionGrower(chm).run_all(
            tops, variance_thresh=60.0, mask_thresh=1.0, protect_seeds=True)
        assert len(np.unique(out)) - 1 == 2

    def test_parallel_matches_sequential(self, stand):
        chm, tops = stand
        kw = dict(variance_thresh=6.0, mask_thresh=1.0)
        seq = HierarchicalRegionGrower(chm).run_all(tops, n_jobs=1, **kw)
        par = HierarchicalRegionGrower(chm).run_all(tops, n_jobs=2, **kw)
        np.testing.assert_array_equal(seq, par)

    def test_retry_rejected_runs(self, stand):
        chm, tops = stand
        out = HierarchicalRegionGrower(chm).run_all(
            tops, variance_thresh=6.0, mask_thresh=1.0, retry_rejected=True)
        assert out.shape == chm.shape

    def test_morpho_radius_runs(self, stand):
        chm, tops = stand
        out = HierarchicalRegionGrower(chm).run_all(
            tops, variance_thresh=6.0, mask_thresh=1.0, morpho_radius=2)
        assert out.shape == chm.shape

    def test_background_is_zero(self, stand):
        chm, tops = stand
        out = HierarchicalRegionGrower(chm).run_all(tops, variance_thresh=6.0,
                                                    mask_thresh=5.0)
        assert (out[chm < 1.0] == 0).all()
