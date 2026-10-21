"""CLI test of scripts/verify_target_pairs.py on a tiny SYNTHETIC payload.

``--contract-dir`` accepts either the release root (containing
``enrichment_significant/``) or that directory itself; the 2026-09-15 candidate
build passed the root and failed with ``FileNotFoundError`` on
``<root>/hsa/GO_BP_HIGH.tsv``. The fixture covers every organism x family x tier
cell the verifier iterates, with stored overlaps computed by hand below, and
checks that a wrong stored overlap is reported.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "vtp_under_test", _REPO / "scripts" / "verify_target_pairs.py"
)
vtp = importlib.util.module_from_spec(_spec)
sys.modules["vtp_under_test"] = vtp
_spec.loader.exec_module(vtp)

# Per organism: one miRNA with one target per atomic band. The Low target is
# spelled upper-case (TargetScan style) and is re-cased onto the family space.
PAIRS = (("Gene1", "High"), ("Gene2", "Medium"), ("GENE3", "Low"))
UNIVERSE = ("Gene1", "Gene2", "Gene3", "Gene4")
TERM2GENE = {"P:1": ("Gene1", "Gene2", "Gene3", "Gene4"), "P:2": ("Gene2",)}
# Stored overlaps = |term & universe & cumulative targets| (0-overlap rows are not
# retained in a contract cell, so HIGH has no P:2 row).
OVERLAPS = {
    "HIGH": {"P:1": 1},
    "MEDIUM": {"P:1": 2, "P:2": 1},
    "LOW": {"P:1": 3, "P:2": 1},
}
HEADER = "organism\tmirna\ttarget_gene\tatomic_evidence\tsource\tsource_annotation"


def _mirna(org: str) -> str:
    return f"{org}-miR-TEST1-3p"


def _write_payload(root: Path) -> Path:
    for org, families in vtp.FAMILIES_BY_ORG.items():
        pairs_dir = root / "target_pairs"
        pairs_dir.mkdir(parents=True, exist_ok=True)
        lines = [HEADER] + [
            f"{org}\t{_mirna(org)}\t{gene}\t{evidence}\tSYNTHETIC\tSYNTHETIC:x"
            for gene, evidence in PAIRS
        ]
        (pairs_dir / f"{org}_target_pairs.tsv").write_text("\n".join(lines) + "\n")
        inputs = root / "enrichment_inputs" / org
        cells = root / "enrichment_significant" / org
        inputs.mkdir(parents=True)
        cells.mkdir(parents=True)
        for family in families:
            (inputs / f"{family}_universe.tsv").write_text("gene\n" + "\n".join(UNIVERSE) + "\n")
            t2g = ["pathway_id\tgene"] + [
                f"{pid}\t{gene}" for pid, genes in TERM2GENE.items() for gene in genes
            ]
            (inputs / f"{family}_term2gene.tsv").write_text("\n".join(t2g) + "\n")
            for tier, overlaps in OVERLAPS.items():
                rows = ["mirna\tpathway_id\tp_value\tp_adjust\toverlap"] + [
                    f"{_mirna(org)}\t{pid}\t0.01\t0.02\t{k}" for pid, k in overlaps.items()
                ]
                (cells / f"{family}_{tier}.tsv").write_text("\n".join(rows) + "\n")
    return root


def _n_cells() -> tuple[int, int]:
    n_cells = sum(len(f) for f in vtp.FAMILIES_BY_ORG.values()) * len(OVERLAPS)
    n_rows = sum(len(f) for f in vtp.FAMILIES_BY_ORG.values()) * sum(
        len(o) for o in OVERLAPS.values()
    )
    return n_cells, n_rows


@pytest.fixture()
def release(tmp_path: Path) -> Path:
    return _write_payload(tmp_path / "release")


@pytest.mark.parametrize("contract", ["root", "enrichment_significant"])
def test_contract_dir_accepts_root_or_enrichment_significant(
    release: Path, tmp_path: Path, contract: str
) -> None:
    contract_dir = release if contract == "root" else release / "enrichment_significant"
    result = tmp_path / f"verify_{contract}.json"
    rc = vtp.main(
        ["--pairs-dir", str(release), "--contract-dir", str(contract_dir), "--result", str(result)]
    )
    assert rc == 0
    report = json.loads(result.read_text())
    n_cells, n_rows = _n_cells()
    assert report["status"] == "done"
    assert report["contract_dir"] == str(release / "enrichment_significant")
    assert len(report["cells"]) == n_cells == 57
    assert report["total_rows"] == n_rows
    assert report["total_mismatches"] == 0


def test_wrong_overlap_is_reported(release: Path, tmp_path: Path) -> None:
    cell = release / "enrichment_significant" / "hsa" / "GO_BP_LOW.tsv"
    text = cell.read_text()
    assert "\tP:1\t0.01\t0.02\t3\n" in text
    cell.write_text(text.replace("\tP:1\t0.01\t0.02\t3\n", "\tP:1\t0.01\t0.02\t2\n"))
    result = tmp_path / "verify_bad.json"
    rc = vtp.main(
        ["--pairs-dir", str(release), "--contract-dir", str(release), "--result", str(result)]
    )
    assert rc == 1
    report = json.loads(result.read_text())
    assert report["status"] == "FAILED"
    assert report["total_mismatches"] == 1
    bad = report["cells"]["hsa/GO_BP_LOW"]
    assert bad["mismatches"] == 1
    assert bad["examples"] == [
        {"mirna": "hsa-miR-TEST1-3p", "pathway_id": "P:1", "stored": 2, "got": 3}
    ]


def test_contract_dir_without_cells_is_a_usage_error(release: Path, tmp_path: Path) -> None:
    empty = tmp_path / "not_a_release"
    empty.mkdir()
    with pytest.raises(SystemExit) as excinfo:
        vtp.main(
            [
                "--pairs-dir",
                str(release),
                "--contract-dir",
                str(empty),
                "--result",
                str(tmp_path / "unused.json"),
            ]
        )
    assert excinfo.value.code == 2
    assert not (tmp_path / "unused.json").exists()
