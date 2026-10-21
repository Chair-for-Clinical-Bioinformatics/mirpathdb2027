"""Tests for scripts/build_release_census.py (REL-001/REL-002)."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "build_release_census",
    Path(__file__).resolve().parents[1] / "scripts" / "build_release_census.py",
)
brc = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(brc)


def _write_cell(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mirna", "pathway_id", "p_value", "p_adjust", "overlap"])
        w.writerows(rows)


@pytest.fixture()
def payload(tmp_path: Path) -> Path:
    p = tmp_path / "mirpathdb-data-test"
    _write_cell(p / "enrichment_significant" / "hsa" / "GO_BP_HIGH.tsv", [
        ["hsa-miR-1", "GO:1", "0.001", "0.004", "3"],   # stringent
        ["hsa-miR-1", "GO:2", "0.01", "0.04", "2"],     # significant
        ["hsa-miR-2", "GO:1", "0.1", "0.2", "1"],       # retained only
    ])
    _write_cell(p / "enrichment_significant" / "hsa" / "GO_BP_LOW.tsv", [
        ["hsa-miR-1", "GO:1", "0.001", "0.03", "5"],
        ["hsa-miR-3", "GO:3", "0.2", "0.24", "1"],
    ])
    _write_cell(p / "enrichment_significant" / "mml" / "REACTOME_LOW.tsv", [
        ["mml-miR-1", "R-MML-1", "0.001", "0.01", "2"],
    ])
    labels = p / "pathway_labels.tsv"
    with labels.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pathway_id", "display_name", "source_db", "family",
                    "source_version", "inferred"])
        w.writerow(["GO:1", "one", "GO", "GO_BP", "v", "false"])
        w.writerow(["GO:2", "two", "GO", "GO_BP", "v", "false"])
        w.writerow(["GO:3", "three", "GO", "GO_BP", "v", "false"])
        w.writerow(["R-MML-1", "inferred one", "Reactome", "REACTOME", "v", "true"])
    (p / "RELEASE_INFO.json").write_text(json.dumps({"release_id": "test-release"}))
    return p


def test_census_counts(payload: Path):
    c = brc.build_census(payload)
    assert c["organisms"] == ["hsa", "mml"]
    assert c["cells_per_organism"] == {"hsa": 2, "mml": 1}
    assert c["cells_total"] == 3
    assert c["associations"]["retained_total"] == 6
    assert c["associations"]["retained_per_organism"] == {"hsa": 5, "mml": 1}
    # rows with p_adjust <= 0.05: 0.004, 0.04, 0.03 (hsa) and 0.01 (mml) -> 4
    assert c["associations"]["significant_total"] == 4
    assert c["associations"]["significant_per_organism"] == {"hsa": 3, "mml": 1}
    # stringent (<= 0.005): 0.004 (hsa) only
    assert c["associations"]["stringent_per_organism"] == {"hsa": 1, "mml": 0}


def test_census_distincts(payload: Path):
    c = brc.build_census(payload)
    # REL-002 definition: distinct miRNAs present in the retained files.
    assert c["mirnas"]["retained_distinct_union"] == 4
    assert c["mirnas"]["retained_distinct_per_organism"] == {"hsa": 3, "mml": 1}
    assert c["mirnas"]["significant_distinct_union"] == 2  # hsa-miR-1, mml-miR-1
    assert c["pathways"]["distinct_accessions_union"] == 4
    assert c["pathways"]["organism_resolved_records"] == 4
    assert c["pathways"]["records_with_significant"] == 3  # GO:1, GO:2, R-MML-1
    assert c["labels"] == {"pathway_label_rows": 4, "inferred_pathway_rows": 1}
    assert c["package"]["data_row_total"] == 6 + 4
    assert c["package"]["enrichment_file_count"] == 3


def test_census_families_and_tiers(payload: Path):
    c = brc.build_census(payload)
    assert c["families"]["retained_rows"] == {"GO_BP": 5, "REACTOME": 1}
    assert c["evidence"]["retained_rows_per_tier"] == {"HIGH": 3, "MEDIUM": 0, "LOW": 3}
    assert c["_meta"]["release_id"] == "test-release"


def test_refuses_non_payload(tmp_path: Path):
    with pytest.raises(SystemExit):
        brc.build_census(tmp_path)


def test_legacy_payload_meta_defaults(payload: Path):
    """A v1 payload (5 columns, no contract_version, no adjustment_family)."""
    c = brc.build_census(payload)
    meta = c["_meta"]
    assert meta["contract_version"] is None
    assert meta["adjustment_family"] == "overlap_family_bh (legacy)"
    assert meta["primary_p_adjust"] == "p_adjust (overlap-family BH, legacy)"
    assert "significant_overlap_family_within_retained_total" not in c["associations"]
    assert "m_family_per_organism_per_family" not in c["families"]


# --------------------------------------------------------------------------- #
# Contract v2 (8-column) payload: counts on p_adjust (= q_full), the legacy
# overlap-family count reported ONLY within the retained set, m_family per cell.
# --------------------------------------------------------------------------- #
V2_HEADER = [
    "mirna", "pathway_id", "p_value", "p_adjust", "overlap",
    "p_adjust_overlap", "m_overlap", "m_family",
]


def _write_v2_cell(path: Path, rows, header=V2_HEADER):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


@pytest.fixture()
def payload_v2(tmp_path: Path) -> Path:
    p = tmp_path / "mirpathdb-data-v2"
    #                mirna        pid     p       q_full   ov  q_overlap m    M
    _write_v2_cell(p / "enrichment_significant" / "hsa" / "GO_BP_HIGH.tsv", [
        ["hsa-miR-1", "GO:1", "0.0001", "0.004", "3", "0.001", "50", "200"],  # stringent both
        ["hsa-miR-1", "GO:2", "0.01", "0.08", "2", "0.02", "50", "200"],      # legacy-only sig
        ["hsa-miR-2", "GO:1", "0.05", "0.2", "1", "0.05", "40", "200"],       # legacy == alpha
        ["hsa-miR-2", "GO:3", "0.06", "0.24", "1", "0.06", "40", "200"],      # retained only
    ])
    _write_v2_cell(p / "enrichment_significant" / "hsa" / "GO_BP_LOW.tsv", [
        ["hsa-miR-3", "GO:1", "0.002", "0.04", "5", "0.01", "100", "200"],    # sig both
    ])
    _write_v2_cell(p / "enrichment_significant" / "mml" / "REACTOME_LOW.tsv", [
        ["mml-miR-1", "R-MML-1", "0.001", "0.01", "2", "0.002", "30", "90"],
    ])
    with (p / "pathway_labels.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pathway_id", "display_name", "source_db", "family",
                    "source_version", "inferred"])
        w.writerow(["GO:1", "one", "GO", "GO_BP", "v", "false"])
        w.writerow(["GO:2", "two", "GO", "GO_BP", "v", "false"])
        w.writerow(["GO:3", "three", "GO", "GO_BP", "v", "false"])
        w.writerow(["R-MML-1", "inferred one", "Reactome", "REACTOME", "v", "true"])
    (p / "contract_version").write_text("2\n")
    (p / "RELEASE_INFO.json").write_text(json.dumps({
        "release_id": "test-v2",
        "adjustment_family": "full_family_bh",
        "primary_p_adjust": "p_adjust (full-family BH over m_family eligible terms)",
    }))
    return p


def test_v2_counts_use_primary_p_adjust(payload_v2: Path):
    c = brc.build_census(payload_v2)
    a = c["associations"]
    assert a["retained_total"] == 6
    # p_adjust <= 0.05: 0.004, 0.04 (hsa), 0.01 (mml) -- NOT the 0.08 row whose
    # legacy value 0.02 would have counted under the old convention.
    assert a["significant_total"] == 3
    assert a["significant_per_organism"] == {"hsa": 2, "mml": 1}
    assert a["stringent_per_organism"] == {"hsa": 1, "mml": 0}
    assert c["families"]["significant_rows"] == {"GO_BP": 2, "REACTOME": 1}
    assert c["evidence"]["significant_rows_per_tier"] == {"HIGH": 1, "MEDIUM": 0, "LOW": 2}


def test_v2_legacy_count_is_within_retained_only(payload_v2: Path):
    c = brc.build_census(payload_v2)
    a = c["associations"]
    # p_adjust_overlap <= 0.05 among retained rows: 0.001, 0.02, 0.05, 0.01 (hsa) + 0.002
    assert a["significant_overlap_family_within_retained_total"] == 5
    assert a["significant_overlap_family_within_retained_per_organism"] == {"hsa": 4, "mml": 1}
    assert a["significant_overlap_family_within_retained_total"] >= a["significant_total"]


def test_v2_m_family_and_meta(payload_v2: Path):
    c = brc.build_census(payload_v2)
    assert c["families"]["m_family_per_organism_per_family"] == {
        "hsa": {"GO_BP": 200},
        "mml": {"REACTOME": 90},
    }
    meta = c["_meta"]
    assert meta["contract_version"] == "2"
    assert meta["release_id"] == "test-v2"
    assert meta["adjustment_family"] == "full_family_bh"
    assert meta["primary_p_adjust"] == "p_adjust (full-family BH over m_family eligible terms)"


def test_v2_meta_derives_primary_from_header_without_release_info(payload_v2: Path):
    (payload_v2 / "RELEASE_INFO.json").unlink()
    meta = brc.build_census(payload_v2)["_meta"]
    assert meta["release_id"] is None
    assert meta["adjustment_family"] == "overlap_family_bh (legacy)"
    assert meta["primary_p_adjust"] == "p_adjust (full-family BH over m_family eligible terms)"


def test_v2_census_body_is_deterministic(payload_v2: Path):
    a = brc.build_census(payload_v2)
    b = brc.build_census(payload_v2)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_v2_refuses_two_m_family_values_in_one_cell(payload_v2: Path):
    _write_v2_cell(payload_v2 / "enrichment_significant" / "mml" / "REACTOME_LOW.tsv", [
        ["mml-miR-1", "R-MML-1", "0.001", "0.01", "2", "0.002", "30", "90"],
        ["mml-miR-2", "R-MML-1", "0.001", "0.01", "2", "0.002", "30", "91"],
    ])
    with pytest.raises(SystemExit, match="not one value per cell"):
        brc.build_census(payload_v2)


def test_v2_refuses_m_family_disagreement_across_tiers(payload_v2: Path):
    _write_v2_cell(payload_v2 / "enrichment_significant" / "hsa" / "GO_BP_LOW.tsv", [
        ["hsa-miR-3", "GO:1", "0.002", "0.04", "5", "0.01", "100", "201"],
    ])
    with pytest.raises(SystemExit, match="differs across tiers"):
        brc.build_census(payload_v2)


def test_refuses_mixed_v1_v2_layouts(payload_v2: Path):
    _write_cell(payload_v2 / "enrichment_significant" / "mml" / "REACTOME_LOW.tsv", [
        ["mml-miR-1", "R-MML-1", "0.001", "0.01", "2"],
    ])
    with pytest.raises(SystemExit, match="mixed contract layouts"):
        brc.build_census(payload_v2)
