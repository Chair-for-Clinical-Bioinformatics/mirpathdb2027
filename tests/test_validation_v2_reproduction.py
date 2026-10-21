"""Tests for the v2-reproduction validation metrics (engine-isolating comparison)."""

from __future__ import annotations

import numpy as np
import pytest

from mirpathdb_prod.validation.v2_reproduction import (
    NOT_TESTED,
    check_frame_concordance,
    compare_cell,
    pooled_confusion,
    pooled_delta,
    pooled_spearman,
)


def _cell(new_vals, v2_vals, *, pathways=None, mirnas=None):
    """Build a comparison from raw value arrays (default identical frame)."""
    new = np.asarray(new_vals, dtype=float)
    v2 = np.asarray(v2_vals, dtype=float)
    n_m, n_p = new.shape
    mirnas = mirnas or [f"mir-{i}" for i in range(n_m)]
    pathways = pathways or [f"path-{j}" for j in range(n_p)]
    return compare_cell(
        organism="hsa",
        variant="validated_miRTarBase_strong",
        family="GO_BP",
        new_pathways=pathways,
        new_mirnas=mirnas,
        new_values=new,
        v2_pathways=pathways,
        v2_mirnas=mirnas,
        v2_values=v2,
    )


def test_frame_concordance_identical_passes():
    fc = check_frame_concordance(
        organism="hsa", variant="v", family="GO_BP",
        new_pathways=["A", "B"], new_mirnas=["m1", "m2"],
        v2_pathways=["A", "B"], v2_mirnas=["m1", "m2"],
    )
    assert fc.passed
    assert fc.pathways_identical and fc.mirnas_identical
    assert fc.matched_pathways == 2 and fc.matched_mirnas == 2
    assert fc.reasons == []


def test_frame_concordance_order_drift_fails():
    fc = check_frame_concordance(
        organism="hsa", variant="v", family="GO_BP",
        new_pathways=["B", "A"], new_mirnas=["m1"],
        v2_pathways=["A", "B"], v2_mirnas=["m1"],
    )
    assert not fc.passed
    assert not fc.pathways_identical
    # same set, different order
    assert fc.matched_pathways == 2
    assert any("ORDER" in r for r in fc.reasons)


def test_frame_concordance_label_set_drift_fails():
    fc = check_frame_concordance(
        organism="hsa", variant="v", family="GO_BP",
        new_pathways=["A", "C"], new_mirnas=["m1"],
        v2_pathways=["A", "B"], v2_mirnas=["m1"],
    )
    assert not fc.passed
    assert fc.v2_only_pathways == ["B"]
    assert fc.new_only_pathways == ["C"]
    assert any("v2-only" in r for r in fc.reasons)


def test_perfect_monotone_gives_rho_one():
    # new is a strictly monotone transform of v2 -> Spearman rho == 1.
    v2 = np.array([[0.5, 0.1, 0.01], [0.2, 0.001, 0.4]])
    new = v2**2  # monotone increasing transform on (0,1)
    c = _cell(new, v2)
    assert c.concordance.passed
    assert c.n_both_tested == 6
    assert c.spearman_rho == pytest.approx(1.0)


def test_not_tested_excluded_from_metrics():
    # 1.0 sentinel cells are "not tested" and must be dropped from both-tested.
    new = np.array([[0.01, NOT_TESTED], [NOT_TESTED, 0.2]])
    v2 = np.array([[0.02, 0.3], [NOT_TESTED, NOT_TESTED]])
    c = _cell(new, v2)
    # both-tested: only [0,0] (0.01 vs 0.02). new-only: [1,1]. v2-only: [0,1].
    assert c.n_both_tested == 1
    assert c.n_new_only_tested == 1
    assert c.n_v2_only_tested == 1
    assert c.n_tested_new == 2 and c.n_tested_v2 == 2


def test_delta_sign_new_more_significant_is_positive():
    # new p smaller => -log10 larger => Delta = new - v2 > 0.
    new = np.array([[0.001, 0.001]])
    v2 = np.array([[0.01, 0.01]])
    c = _cell(new, v2)
    assert c.delta_mean == pytest.approx(1.0)  # log10(0.01)-log10(0.001) = 1
    assert c.delta_median == pytest.approx(1.0)


def test_directional_confusion_raw_and_bh():
    # Construct a clear raw-significance pattern.
    new = np.array([[0.01, 0.5, 0.5, 0.5]])
    v2 = np.array([[0.01, 0.04, 0.5, 0.5]])
    c = _cell(new, v2)
    raw = next(cf for cf in c.confusions if cf.threshold == "raw_p<=0.05")
    # cell0: both sig; cell1: v2-only sig; cell2/3: both-not.
    assert raw.both_sig == 1
    assert raw.v2_only == 1
    assert raw.new_only == 0
    assert raw.both_not == 2
    assert raw.agreement == pytest.approx(3 / 4)
    # BH axis exists.
    assert any(cf.threshold == "bh_fdr<=0.05" for cf in c.confusions)


def test_bh_excludes_not_tested_from_denominator():
    # Many not-tested sentinels must NOT dilute BH (which would let small p
    # survive). Two tested cells, one tiny p.
    new = np.array([[1e-6] + [NOT_TESTED] * 999])
    v2 = np.array([[1e-6] + [NOT_TESTED] * 999])
    c = _cell(new, v2)
    bh = next(cf for cf in c.confusions if cf.threshold == "bh_fdr<=0.05")
    # only 1 tested cell on each side; with n=1 BH == raw, p=1e-6 is significant.
    assert bh.both_sig == 1
    assert bh.n == 1  # domain = either_tested == the single tested cell


def test_shape_mismatch_raises():
    with pytest.raises(ValueError, match="shape mismatch"):
        compare_cell(
            organism="hsa", variant="v", family="GO_BP",
            new_pathways=["A"], new_mirnas=["m1"],
            new_values=np.array([[0.1]]),
            v2_pathways=["A"], v2_mirnas=["m1"],
            v2_values=np.array([[0.1, 0.2]]),
        )


def test_pooled_helpers():
    c1 = _cell(np.array([[0.01, 0.1]]), np.array([[0.02, 0.2]]))
    c2 = _cell(np.array([[0.001, 0.5]]), np.array([[0.002, 0.6]]))
    rho, p, n = pooled_spearman([c1, c2])
    assert n == 4
    assert not np.isnan(rho)
    pc = pooled_confusion([c1, c2], "raw_p<=0.05")
    assert pc.n == pc.both_sig + pc.v2_only + pc.new_only + pc.both_not
    pd = pooled_delta([c1, c2])
    assert pd["n"] == 4
    assert "median" in pd


def test_all_not_tested_yields_nan_rho():
    new = np.full((2, 2), NOT_TESTED)
    v2 = np.full((2, 2), NOT_TESTED)
    c = _cell(new, v2)
    assert c.n_both_tested == 0
    assert np.isnan(c.spearman_rho)
    assert np.isnan(c.delta_median)
