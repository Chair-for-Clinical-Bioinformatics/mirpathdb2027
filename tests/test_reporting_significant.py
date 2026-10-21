"""Tests for the authoritative significant-subset writer (Deliverable B, M1 fix).

The runners now emit Deliverable B straight from run_ora_batch's authoritative
(p, p_adjust, overlap) tuples via write_significant_long_form, so the published
p_adjust is NOT re-derived from the lossy matrix (the source of the M1 defect).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from mirpathdb_prod.enrichment import family_bh
from mirpathdb_prod.reporting.significant import (
    FDR_CUTOFF,
    SIGNIFICANT_HEADER,
    SIGNIFICANT_HEADER_V2,
    count_retained_rows,
    write_significant_long_form,
)

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _bh(p: np.ndarray) -> np.ndarray:
    m = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    return out


def _enricher_results(
    n_terms: dict[str, int], *, seed: int = 0, with_p1: bool = True
) -> dict[str, dict[str, tuple[float, float, int]]]:
    """Unfiltered ``run_ora_batch``-shaped output: BH over exactly the listed terms."""
    rng = np.random.default_rng(seed)
    out: dict[str, dict[str, tuple[float, float, int]]] = {}
    for mirna, m in n_terms.items():
        p = rng.random(m) ** 3
        if with_p1 and m:
            p[-1] = 1.0  # a GO-root-like term: K = N, exact P = 1, still a tested term
        q = _bh(p)
        out[mirna] = {f"T{i}": (float(p[i]), float(q[i]), int(i % 4) + 1) for i in range(m)}
    return out


def test_writes_only_significant_rows_with_authoritative_values(tmp_path: Path) -> None:
    out = tmp_path / "GO_BP_HIGH.tsv"
    pvals = {
        "hsa-miR-1": {
            "GO:1": (1e-6, 0.001, 12),  # significant -> kept
            "GO:2": (0.30, 0.40, 3),  # padj > cutoff -> dropped
        },
        "hsa-miR-2": {
            "GO:3": (0.02, 0.25, 5),  # padj == cutoff -> kept (<=)
        },
    }
    n = write_significant_long_form(out, ["hsa-miR-1", "hsa-miR-2"], pvals)
    assert n == 2
    rows = _rows(out)
    assert [r["pathway_id"] for r in rows] == ["GO:1", "GO:3"]
    # authoritative padj + overlap are written verbatim (NOT re-derived).
    r0 = rows[0]
    assert float(r0["p_value"]) == 1e-6 and float(r0["p_adjust"]) == 0.001
    assert int(r0["overlap"]) == 12
    assert list(rows[0].keys()) == ["mirna", "pathway_id", "p_value", "p_adjust", "overlap"]


def test_fdr_cutoff_default_is_the_shared_constant(tmp_path: Path) -> None:
    assert FDR_CUTOFF == 0.25
    out = tmp_path / "x.tsv"
    # A row exactly at a custom cutoff boundary is kept; just above is dropped.
    n = write_significant_long_form(
        out,
        ["m"],
        {"m": {"P1": (0.1, 0.10, 2), "P2": (0.1, 0.11, 2)}},
        fdr_cutoff=0.10,
    )
    assert n == 1


def test_empty_input_writes_header_only(tmp_path: Path) -> None:
    out = tmp_path / "empty.tsv"
    n = write_significant_long_form(out, [], {})
    assert n == 0
    assert _rows(out) == []  # header only, no data rows


def test_v2_header_extends_legacy_header_as_prefix() -> None:
    assert SIGNIFICANT_HEADER_V2[: len(SIGNIFICANT_HEADER)] == SIGNIFICANT_HEADER
    assert SIGNIFICANT_HEADER_V2[5:] == ("p_adjust_overlap", "m_overlap", "m_family")


def test_v2_retains_on_full_family_and_writes_both_adjusted_values(tmp_path: Path) -> None:
    res = _enricher_results({"m1": 40, "m2": 120, "m3": 7}, seed=5)
    m_family = 500
    out = tmp_path / "GO_BP_HIGH.tsv"
    n = write_significant_long_form(out, ["m1", "m2", "m3"], res, m_family=m_family)
    rows = _rows(out)
    assert len(rows) == n
    assert list(rows[0].keys()) == list(SIGNIFICANT_HEADER_V2)

    adjusted = family_bh.adjust_cell_results(res, m_family)
    expected = {
        (mirna, t): v
        for mirna, per in adjusted.items()
        for t, v in per.items()
        if v[2] <= FDR_CUTOFF
    }
    assert n == len(expected) > 0
    for r in rows:
        p, q_overlap, q_full, overlap = expected[(r["mirna"], r["pathway_id"])]
        assert float(r["p_value"]) == p and int(r["overlap"]) == overlap
        assert float(r["p_adjust"]) == q_full
        assert float(r["p_adjust_overlap"]) == q_overlap
        assert float(r["p_adjust"]) <= FDR_CUTOFF
        # q_full >= q_overlap on every row.
        assert float(r["p_adjust"]) >= float(r["p_adjust_overlap"])
        # m_overlap = every listed term of that miRNA, the P = 1 term included.
        assert int(r["m_overlap"]) == len(res[r["mirna"]])
        assert int(r["m_family"]) == m_family
    # The P = 1 term is part of m_overlap for each miRNA.
    for mirna, m in (("m1", 40), ("m2", 120), ("m3", 7)):
        assert any(v[0] == 1.0 for v in res[mirna].values())
        assert {int(r["m_overlap"]) for r in rows if r["mirna"] == mirna} <= {m}
    # Floats keep their repr.
    with out.open(encoding="utf-8") as fh:
        first = fh.readlines()[1].rstrip("\r\n").split("\t")
    (p, q_overlap, q_full, _ov) = expected[(first[0], first[1])]
    assert first[2:4] == [repr(p), repr(q_full)] and first[5] == repr(q_overlap)


def test_v2_refuses_a_pre_filtered_family(tmp_path: Path) -> None:
    res = _enricher_results({"m1": 200}, seed=2)
    filtered = {"m1": {t: v for t, v in res["m1"].items() if v[1] <= FDR_CUTOFF}}
    assert 0 < len(filtered["m1"]) < 200
    with pytest.raises(ValueError, match="refusing to migrate"):
        write_significant_long_form(tmp_path / "x.tsv", ["m1"], filtered, m_family=15046)
    with pytest.raises(ValueError, match="refusing to migrate"):
        count_retained_rows(filtered, fdr_cutoff=FDR_CUTOFF, m_family=15046)


def test_singleton_counterexample_dropped_by_v2_kept_by_legacy(tmp_path: Path) -> None:
    """One overlapping term out of 100: overlap-only q = 0.01, full-family q = 1.0."""
    res = {"m": {"T0": (0.01, 0.01, 1)}}
    legacy = tmp_path / "legacy.tsv"
    v2 = tmp_path / "v2.tsv"
    assert write_significant_long_form(legacy, ["m"], res) == 1
    assert write_significant_long_form(v2, ["m"], res, m_family=100) == 0
    assert _rows(v2) == [] and _rows(legacy)[0]["p_adjust"] == "0.01"
    counts = count_retained_rows(res, fdr_cutoff=FDR_CUTOFF, m_family=100)
    assert counts == {
        "rows": 0,
        "rows_overlap_legacy": 1,
        "m_overlap_min": 1,
        "m_overlap_max": 1,
    }


def test_count_retained_rows_matches_writer_and_counts_m_overlap(tmp_path: Path) -> None:
    res = _enricher_results({"a": 30, "b": 300, "c": 0, "d": 12}, seed=9)
    m_family = 2000
    n = write_significant_long_form(
        tmp_path / "x.tsv", ["a", "b", "c", "d"], res, m_family=m_family
    )
    counts = count_retained_rows(res, fdr_cutoff=FDR_CUTOFF, m_family=m_family)
    assert counts["rows"] == n
    assert counts["rows_overlap_legacy"] == write_significant_long_form(
        tmp_path / "legacy.tsv", ["a", "b", "c", "d"], res
    )
    assert counts["rows_overlap_legacy"] >= counts["rows"]
    # min/max over scored miRNAs only (c has no listed term).
    assert (counts["m_overlap_min"], counts["m_overlap_max"]) == (12, 300)
    assert count_retained_rows({"c": {}}, fdr_cutoff=FDR_CUTOFF, m_family=5) == {
        "rows": 0,
        "rows_overlap_legacy": 0,
        "m_overlap_min": None,
        "m_overlap_max": None,
    }
    # Legacy mode counts overlap-only retention in both slots.
    legacy_counts = count_retained_rows(res, fdr_cutoff=FDR_CUTOFF)
    assert legacy_counts["rows"] == legacy_counts["rows_overlap_legacy"]


def test_legacy_mode_is_unchanged_when_m_family_is_none(tmp_path: Path) -> None:
    res = _enricher_results({"m1": 50, "m2": 9}, seed=4)
    out = tmp_path / "legacy.tsv"
    n = write_significant_long_form(out, ["m1", "m2"], res)
    text = out.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[0].split("\t") == list(SIGNIFICANT_HEADER)
    assert len(lines) == n + 1
    expected = [
        f"{mirna}\t{t}\t{p!r}\t{q!r}\t{ov}"
        for mirna in ("m1", "m2")
        for t, (p, q, ov) in res[mirna].items()
        if q <= FDR_CUTOFF
    ]
    assert lines[1:] == expected and text.endswith("\n")


def test_all_runners_emit_authoritative_deliverable_b() -> None:
    """All three production runners call the shared authoritative writer, and the
    main runner defaults long-form ON (--no-long-form to opt out)."""
    src = {
        name: (_SCRIPTS / f"{name}.py").read_text(encoding="utf-8")
        for name in (
            "run_production_pathways",
            "run_production_enrichment",
            "run_macaque_inferred_reactome",
        )
    }
    for name, text in src.items():
        assert "write_significant_long_form" in text, name
    # run_production_pathways: long-form is now default-on via --no-long-form.
    assert "--no-long-form" in src["run_production_pathways"]
    assert "--long-form" not in src["run_production_pathways"].replace("--no-long-form", "")
