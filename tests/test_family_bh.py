"""Full-family BH: eligibility, direct adjustment, exact migration, refusal of filtered families."""

from __future__ import annotations

import numpy as np
import pytest

from mirpathdb_prod.enrichment import family_bh as fb


def _bh(p):
    p = np.asarray(p, dtype=float)
    m = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    return out


def test_eligible_family_counts_non_empty_terms_inside_background():
    t2g = {"A": ["g1", "g2"], "B": ["g9"], "C": [], "D": ["g3", "zz"], "E": ["g1"] * 3}
    m_full, terms = fb.eligible_family(t2g, ["g1", "g2", "g3"])
    assert (m_full, terms) == (3, ["A", "D", "E"])  # B outside, C empty
    assert fb.eligible_family(t2g, ["g1", "g2", "g3"], max_gs_size=1) == (2, ["D", "E"])


def test_full_family_equals_padding_with_ones():
    rng = np.random.default_rng(0)
    p = rng.random(50) ** 2
    direct = fb.full_family_adjust(p, 400)
    padded = _bh(np.concatenate([p, np.ones(350)]))
    assert np.allclose(direct, padded[:50], rtol=1e-13, atol=0)
    assert np.all(padded[50:] == 1.0)
    # With M == m it is plain BH.
    assert np.allclose(fb.full_family_adjust(p, 50), _bh(p), rtol=1e-13)


def test_migration_is_exact_when_m_is_the_complete_family():
    rng = np.random.default_rng(1)
    for m, m_full in ((1, 100), (7, 40), (300, 15046), (66, 15046)):
        p = np.sort(rng.random(m) ** 3)
        p[-1] = 1.0  # a valid P = 1 test inside the overlapping family
        q_overlap = _bh(p)
        migrated = fb.migrate_overlap_adjusted(q_overlap, m, m_full)
        direct = fb.full_family_adjust(p, m_full)
        assert np.allclose(migrated, direct, rtol=1e-13, atol=0)
        rep = fb.verify_migration(p, q_overlap, m_full)
        assert rep["passed"] and rep["m_overlap"] == m


def test_singleton_counterexample():
    """100 singleton terms, one-gene query: overlap-only 0.01, full-family 1.0."""
    p = np.array([0.01])
    assert fb.full_family_adjust(p, 1)[0] == pytest.approx(0.01)
    assert fb.full_family_adjust(p, 100)[0] == pytest.approx(1.0)
    assert fb.migrate_overlap_adjusted([0.01], 1, 100)[0] == pytest.approx(1.0)


def test_migration_refuses_retained_row_count():
    """Using a filtered (retained) subset as m is detected: q_overlap is not BH over it."""
    rng = np.random.default_rng(2)
    p_all = np.sort(rng.random(200) ** 2)
    q_all = _bh(p_all)
    keep = q_all <= 0.25  # the retained subset
    assert 0 < keep.sum() < 200
    rep = fb.verify_migration(p_all[keep], q_all[keep], 15046)
    assert rep["passed"] is False
    with pytest.raises(ValueError, match="refusing to migrate"):
        fb.adjust_cell_results(
            {"mir": {f"t{i}": (float(p_all[i]), float(q_all[i]), 1) for i in np.flatnonzero(keep)}},
            15046,
        )


def test_adjust_cell_results_passes_raw_p_and_overlap_through():
    rng = np.random.default_rng(3)
    p = rng.random(30)
    q = _bh(p)
    res = {
        "mirA": {f"t{i}": (float(p[i]), float(q[i]), int(i % 5) + 1) for i in range(30)},
        "mirB": {},
    }
    out = fb.adjust_cell_results(res, 2000)
    direct = fb.full_family_adjust(p, 2000)
    for i in range(30):
        pp, qo, qf, k = out["mirA"][f"t{i}"]
        assert pp == p[i] and qo == q[i] and k == i % 5 + 1
        assert qf == pytest.approx(direct[i], rel=1e-13)
        assert qf >= qo - 1e-15
    assert out["mirB"] == {}


def test_input_validation():
    with pytest.raises(ValueError):
        fb.full_family_adjust([0.1, 0.2], 1)
    with pytest.raises(ValueError):
        fb.full_family_adjust([1.5], 3)
    with pytest.raises(ValueError):
        fb.migrate_overlap_adjusted([0.1], 0, 10)
    with pytest.raises(ValueError):
        fb.migrate_overlap_adjusted([0.1], 5, 4)
    assert fb.full_family_adjust([], 10).size == 0
