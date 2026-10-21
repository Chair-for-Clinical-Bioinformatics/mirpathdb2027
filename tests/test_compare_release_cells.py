"""Tests for scripts/validation/compare_release_cells.py (old-versus-new cell census)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validation" / "compare_release_cells.py"

V1_HEADER = "mirna\tpathway_id\tp_value\tp_adjust\toverlap\n"
V2_HEADER = V1_HEADER.rstrip("\n") + "\tp_adjust_overlap\tm_overlap\tm_family\n"


def _load():
    spec = importlib.util.spec_from_file_location("compare_release_cells", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _payload(root: Path, *, v2: bool) -> None:
    (root / "enrichment_inputs" / "hsa").mkdir(parents=True)
    (root / "enrichment_inputs" / "hsa" / "GO_BP_term2gene.tsv").write_text(
        "pathway_id\tgene\nGO:1\tA\nGO:1\tB\nGO:2\tB\nGO:3\tC\n", encoding="utf-8"
    )
    cell = root / "enrichment_significant" / "hsa"
    cell.mkdir(parents=True)
    if v2:
        rows = [
            ("m1", "GO:1", 1e-4, 0.004, 2, 0.001, 2, 3),
            ("m1", "GO:2", 1e-2, 0.03, 1, 0.01, 2, 3),
            ("m2", "GO:1", 5e-2, 0.2, 1, 0.05, 1, 3),
        ]
        body = "".join("\t".join(map(str, r)) + "\n" for r in rows)
        (cell / "GO_BP_HIGH.tsv").write_text(V2_HEADER + body, encoding="utf-8")
    else:
        rows = [
            ("m1", "GO:1", 1e-4, 0.001, 2),
            ("m1", "GO:2", 1e-2, 0.01, 1),
            ("m2", "GO:1", 5e-2, 0.05, 1),
            ("m3", "GO:3", 1e-1, 0.2, 1),
        ]
        body = "".join("\t".join(map(str, r)) + "\n" for r in rows)
        (cell / "GO_BP_HIGH.tsv").write_text(V1_HEADER + body, encoding="utf-8")


def test_counts_deltas_and_inventory(tmp_path: Path) -> None:
    old, new, out = tmp_path / "old", tmp_path / "new", tmp_path / "out"
    _payload(old, v2=False)
    _payload(new, v2=True)
    mod = _load()
    assert mod.main(["--old", str(old), "--new", str(new), "--out", str(out)]) == 0
    cells = pd.read_csv(out / "cells.tsv", sep="\t")
    assert len(cells) == 1
    r = cells.iloc[0]
    assert (r.old_rows, r.new_rows, r.delta_rows) == (4, 3, -1)
    assert (r.old_rows_significant, r.new_rows_significant) == (3, 2)
    assert (r.old_rows_stringent, r.new_rows_stringent) == (1, 1)
    assert (r.old_pathways_retained, r.new_pathways_retained) == (3, 2)
    assert (r.inventory_pathways_old, r.inventory_pathways_new) == (3, 3)
    # GO:3 lost its only retained row: it is retained-empty, still part of the inventory.
    assert r.new_pathways_retained_empty == 1
    assert r.old_pathways_retained_empty == 0
    assert r.new_m_family == 3
    assert r.new_rows_significant_overlap_family_within_retained == 3
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert summary["inventory_unchanged"] is True
    assert summary["totals"]["delta_rows"] == -1
    assert summary["cells_missing_in_new"] == []


def test_empty_payloads_are_an_error(tmp_path: Path) -> None:
    (tmp_path / "old" / "enrichment_significant").mkdir(parents=True)
    (tmp_path / "new" / "enrichment_significant").mkdir(parents=True)
    mod = _load()
    with pytest.raises(SystemExit):
        mod.main(
            ["--old", str(tmp_path / "old"), "--new", str(tmp_path / "new"), "--out", str(tmp_path)]
        )
