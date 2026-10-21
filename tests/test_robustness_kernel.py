"""Independent kernel: hypergeometric tail, the two BH families, batch semantics."""

from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.sparse as sp

from mirpathdb_prod.validation.robustness.inputs import FamilySpace
from mirpathdb_prod.validation.robustness.kernel import (
    bh_adjust,
    bh_fixed_family,
    fold_enrichment,
    hypergeom_sf,
    hypergeom_sf_unique,
    ora_batch,
    ora_single,
)


# ---- reference implementations (independent of scipy) ----------------------
def _log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def ref_sf(k: int, N: int, K: int, n: int) -> float:
    lo = max(0, n - (N - K))
    hi = min(K, n)
    start = max(k, lo)
    if start > hi:
        return 0.0
    denom = _log_choose(N, n)
    return min(
        1.0,
        sum(
            math.exp(_log_choose(K, x) + _log_choose(N - K, n - x) - denom)
            for x in range(start, hi + 1)
        ),
    )


def ref_bh(p: list[float]) -> list[float]:
    m = len(p)
    order = sorted(range(m), key=lambda i: p[i])
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, p[i] * m / (rank + 1))
        adj[i] = min(1.0, prev)
    return adj


def _space(terms: dict[str, list[str]], genes: list[str]) -> FamilySpace:
    gi = {g: i for i, g in enumerate(genes)}
    names = sorted(terms)
    rows, cols = [], []
    for j, t in enumerate(names):
        for g in terms[t]:
            rows.append(j)
            cols.append(gi[g])
    mat = sp.coo_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(len(names), len(genes))
    ).tocsr()
    return FamilySpace(
        organism="syn", family="SYN", genes=tuple(genes), terms=tuple(names), matrix=mat
    )


# ---- hypergeometric tail -----------------------------------------------------
@pytest.mark.parametrize(
    "N,K,n,k",
    [(50, 10, 5, 1), (50, 10, 5, 3), (21697, 35, 3, 2), (1000, 300, 200, 50), (100, 100, 7, 7)],
)
def test_hypergeom_sf_matches_reference(N, K, n, k):
    got = float(hypergeom_sf(np.array([k]), N, np.array([K]), n)[0])
    assert got == pytest.approx(ref_sf(k, N, K, n), rel=1e-11)


def test_hypergeom_sf_zero_overlap_is_exactly_one():
    assert hypergeom_sf(np.array([0, 0]), 100, np.array([5, 50]), 10).tolist() == [1.0, 1.0]


def test_worked_release_cell_reproduces_stored_raw_p():
    """rno-let-7a-5p x GO:0060252: N=21697, K=35, n=3, k=2 (recorded 7.57614418739981e-06)."""
    got = float(hypergeom_sf(np.array([2]), 21697, np.array([35]), 3)[0])
    assert got == pytest.approx(7.57614418739981e-06, rel=1e-9)


def test_hypergeom_sf_unique_equals_direct():
    rng = np.random.default_rng(0)
    K = rng.integers(1, 200, size=500)
    k = np.minimum(rng.integers(0, 30, size=500), K)
    a = hypergeom_sf_unique(k, K, 40, 2000)
    b = hypergeom_sf(k, 2000, K, 40)
    assert np.allclose(a, b, rtol=0, atol=0)


# ---- BH families ----------------------------------------------------------------
def test_bh_adjust_matches_reference_and_statsmodels():
    rng = np.random.default_rng(1)
    p = np.concatenate([rng.random(40), [1.0, 1.0, 0.5, 0.5, 1e-12]])
    got = bh_adjust(p)
    assert np.allclose(got, ref_bh(p.tolist()), rtol=1e-12)
    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError:  # optional
        return
    assert np.allclose(got, multipletests(p, method="fdr_bh")[1], rtol=1e-12)


def test_bh_fixed_family_equals_padding_with_ones():
    rng = np.random.default_rng(2)
    p = rng.random(30) ** 3
    m_full = 120
    padded = bh_adjust(np.concatenate([p, np.ones(m_full - 30)]))
    assert np.allclose(bh_fixed_family(p, m_full), padded[:30], rtol=1e-12)
    assert np.all(padded[30:] == 1.0)


def test_bh_fixed_family_with_full_listing_equals_bh():
    p = np.array([0.01, 0.2, 0.03, 1.0])
    assert np.allclose(bh_fixed_family(p, 4), bh_adjust(p))
    with pytest.raises(ValueError):
        bh_fixed_family(p, 3)


def test_bh_empty():
    assert bh_adjust(np.zeros(0)).size == 0
    assert bh_fixed_family(np.zeros(0), 10).size == 0


# ---- the labelled synthetic counterexample ----------------------------------------
def test_singleton_counterexample_overlap_only_vs_fixed_family():
    """100 singleton terms over 100 genes, one-gene query.

    Overlap-only BH corrects over the single term hit -> adjusted P = raw P = 1/100.
    Fixed-family BH corrects over all 100 eligible terms -> adjusted P = 1.
    Regression fixture for the two correction families; not a claim about the release.
    """
    genes = [f"g{i}" for i in range(100)]
    terms = {f"t{i}": [genes[i]] for i in range(100)}
    space = _space(terms, genes)
    res = ora_single(space, np.array([7]))
    assert res.m_full == 100
    assert res.hits[0].size == 1
    assert res.p[0][0] == pytest.approx(0.01)
    assert res.padj_overlap[0][0] == pytest.approx(0.01)
    assert res.padj_fixed[0][0] == pytest.approx(1.0)


# ---- batch semantics ---------------------------------------------------------------
def _demo_space() -> FamilySpace:
    genes = [f"g{i:02d}" for i in range(30)]
    terms = {
        "A": genes[0:5],
        "B": genes[3:12],
        "C": genes[10:30],
        "ALL": genes[:],
        "D": genes[28:30],
    }
    return _space(terms, genes + ["orphanX", "orphanY"])


def test_zero_overlap_query_scores_nothing():
    space = _demo_space()
    res = ora_single(space, np.array([30, 31]))  # the two orphans -> no term hit
    assert res.n[0] == 2
    assert res.hits[0].size == 0
    assert res.n_rejections(0.05, "overlap").tolist() == [0]
    assert res.n_rejections(0.05, "fixed").tolist() == [0]


def test_empty_query():
    space = _demo_space()
    res = ora_single(space, np.zeros(0, dtype=np.int64))
    assert res.n[0] == 0 and res.hits[0].size == 0


def test_whole_universe_term_gives_valid_p_one_and_is_in_family():
    genes = [f"g{i:02d}" for i in range(30)]
    space = _space({"A": genes[0:5], "ALL": genes[:]}, genes)  # ALL == the whole background
    res = ora_single(space, np.array([0, 1, 2]))
    d = res.as_dict(0, space.terms)
    assert "ALL" in d
    p_all, padj_all, k_all = d["ALL"]
    assert k_all == 3 and p_all == 1.0 and padj_all == 1.0
    # P = 1 terms count toward m (overlap family) -- 'A' and 'ALL' are hit here.
    assert res.hits[0].size == 2
    # ...and they are NOT excluded from the fixed family either (m_full = 2 here).
    assert res.padj_fixed[0][list(space.terms).index("ALL")] == 1.0


def test_batch_matches_single_and_chunking_is_invariant():
    space = _demo_space()
    rng = np.random.default_rng(3)
    queries = [np.sort(rng.choice(space.N, size=s, replace=False)) for s in (1, 3, 5, 8, 13)]
    b1 = ora_batch(space, queries, chunk=1)
    b2 = ora_batch(space, queries, chunk=512)
    for i, q in enumerate(queries):
        s = ora_single(space, q)
        for a, b, c in zip(
            (b1.hits[i], b1.p[i], b1.padj_overlap[i], b1.padj_fixed[i]),
            (b2.hits[i], b2.p[i], b2.padj_overlap[i], b2.padj_fixed[i]),
            (s.hits[0], s.p[0], s.padj_overlap[0], s.padj_fixed[0]),
            strict=True,
        ):
            assert np.array_equal(a, b) and np.array_equal(a, c)


def test_batch_values_against_reference_by_hand():
    space = _demo_space()
    q = np.array([0, 1, 2, 3, 11, 29])  # hits A(4), B(2), C(2), ALL(6), D(1)
    res = ora_single(space, q)
    d = res.as_dict(0, space.terms)
    N, n = space.N, 6
    sizes = dict(zip(space.terms, space.term_sizes.tolist(), strict=True))
    expect = {"A": 4, "B": 2, "C": 2, "ALL": 6, "D": 1}
    assert {t: v[2] for t, v in d.items()} == expect
    raw = {t: ref_sf(k, N, sizes[t], n) for t, k in expect.items()}
    for t, p in raw.items():
        assert d[t][0] == pytest.approx(p, rel=1e-11)
    adj = dict(zip(expect, ref_bh([raw[t] for t in expect]), strict=True))
    for t in expect:
        assert d[t][1] == pytest.approx(adj[t], rel=1e-11)
    fe = fold_enrichment(np.array([4]), n, np.array([sizes["A"]]), N)[0]
    assert fe == pytest.approx((4 / 6) / (5 / N))


def test_query_outside_background_is_rejected():
    space = _demo_space()
    with pytest.raises(ValueError):
        ora_single(space, np.array([space.N + 5]))


def test_large_overlaps_do_not_overflow():
    """Regression: overlaps above 127 must not wrap (an int8 matrix product did)."""
    genes = [f"g{i:03d}" for i in range(400)]
    space = _space({"BIG": genes[:300], "SMALL": genes[:10]}, genes)
    res = ora_single(space, np.arange(250))
    d = res.as_dict(0, space.terms)
    assert d["BIG"][2] == 250 and d["SMALL"][2] == 10
    assert d["BIG"][0] == pytest.approx(ref_sf(250, 400, 300, 250), rel=1e-11)
