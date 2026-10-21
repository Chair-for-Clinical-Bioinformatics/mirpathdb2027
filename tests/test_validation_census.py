"""Tests for the census / matrix-non-degeneracy validator (reconciliation gate)."""

from __future__ import annotations

import json
from pathlib import Path

from mirpathdb_prod.validation import census


def _matrix(path: Path, header: list[str], rows: list[list]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def test_scan_matrix_nondegenerate(tmp_path: Path) -> None:
    p = tmp_path / "m.csv"
    _matrix(p, ["miRNA", "P1", "P2"], [["mir-1", "0.01", "1.0"], ["mir-2", "1.0", "1.0"]])
    r = census.scan_matrix(p)
    assert r["non_degenerate"] is True and r["n_cols"] == 2 and r["n_rows"] == 2


def test_scan_matrix_all_ones_is_degenerate(tmp_path: Path) -> None:
    p = tmp_path / "m.csv"
    _matrix(p, ["miRNA", "P1", "P2"], [["mir-1", "1.0", "1.0"], ["mir-2", "1.0", "1.0"]])
    assert census.scan_matrix(p)["non_degenerate"] is False


def test_scan_matrix_header_only_is_degenerate(tmp_path: Path) -> None:
    p = tmp_path / "m.csv"
    _matrix(p, ["miRNA"], [])
    assert census.scan_matrix(p)["non_degenerate"] is False


def test_validate_matrices_flags_only_degenerate(tmp_path: Path) -> None:
    root = tmp_path / "enrichment"
    _matrix(root / "hsa" / "good.csv", ["miRNA", "P1"], [["m", "0.001"]])
    _matrix(root / "mmu" / "bad.csv", ["miRNA", "P1"], [["m", "1.0"]])
    res = census.validate_matrices_nondegenerate(root)
    assert res["ok"] is False and res["n_matrices"] == 2
    assert [Path(d["path"]).name for d in res["degenerate"]] == ["bad.csv"]


def test_validate_census_anchor_pass_and_fail(tmp_path: Path) -> None:
    prov = tmp_path / "prov"
    prov.mkdir()
    good = {
        "inferred_reactome_summary": {"n_terms": 2777},
        "target_set_census": {
            "HIGH": {"n_mirnas": 1},
            "MEDIUM": {"n_mirnas": 0},
            "LOW": {"n_mirnas": 293},
        },
    }
    (prov / "macaque_inferred_reactome_manifest.json").write_text(json.dumps(good))
    assert census.validate_census(prov)["ok"] is True

    bad = json.loads(json.dumps(good))
    bad["target_set_census"]["HIGH"]["n_mirnas"] = 0  # macaque lone MTI vanished
    (prov / "macaque_inferred_reactome_manifest.json").write_text(json.dumps(bad))
    res = census.validate_census(prov)
    assert res["ok"] is False
    assert any(f["anchor"] == "macaque_HIGH_mirnas" for f in res["failures"])


# --------------------------------------------------------------------------- #
# require_present: the gate must not vacuously pass on absent artefacts.
# --------------------------------------------------------------------------- #
def test_empty_enrichment_dir_vacuous_pass_by_default(tmp_path: Path) -> None:
    """Default (lenient) mode still passes an empty dir -- backward compatible."""
    root = tmp_path / "enrichment"
    root.mkdir()
    res = census.validate_matrices_nondegenerate(root)
    assert res["ok"] is True and res["empty"] is True and res["n_matrices"] == 0


def test_empty_enrichment_dir_fails_under_require_present(tmp_path: Path) -> None:
    """Production gate: an empty enrichment dir is a hard FAILURE, not a pass."""
    root = tmp_path / "enrichment"
    root.mkdir()
    res = census.validate_matrices_nondegenerate(root, require_present=True)
    assert res["ok"] is False and res["empty"] is True


def test_absent_manifest_skips_by_default_but_fails_when_required(tmp_path: Path) -> None:
    """A missing anchor manifest is a skip by default, a FAILURE under the gate."""
    prov = tmp_path / "prov"
    prov.mkdir()  # no manifests at all
    lenient = census.validate_census(prov)
    assert lenient["ok"] is True  # vacuous pass preserved for the default mode
    assert all(c["status"] == "skip" for c in lenient["checks"])

    strict = census.validate_census(prov, require_present=True)
    assert strict["ok"] is False  # absent manifests now fail the gate
    assert strict["failures"]
