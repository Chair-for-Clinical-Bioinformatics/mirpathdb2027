"""Tests for the v2 p-value-matrix reader."""

from __future__ import annotations

import numpy as np
import pytest

from mirpathdb_prod.reproduction.v2_matrix import (
    NOT_TESTED,
    read_v2_matrix,
    v2_matrix_dir,
    v2_matrix_path,
)


def _write_matrix(path, header, rows):
    """Write a v2-shaped CSV: header (no leading label) then miRNA,p,p,... rows."""
    import csv

    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def test_reads_header_data_offset(tmp_path):
    # Header has one fewer field than each data row (the v2 quirk): header col i
    # aligns with data field i+1.
    p = tmp_path / "GO_MF_validated_miRTarBase_strong.csv"
    _write_matrix(
        p,
        ["ATP binding", "DNA binding", "RNA binding"],
        [
            ["hsa-let-7a-5p", "0.0293925", "1.0", "0.5"],
            ["hsa-miR-1", "1.0", "0.01", "1.0"],
        ],
    )
    m = read_v2_matrix(p)
    assert m.pathways == ["ATP binding", "DNA binding", "RNA binding"]
    assert m.mirnas == ["hsa-let-7a-5p", "hsa-miR-1"]
    assert m.values.shape == (2, 3)
    # ATP binding (col 0) for let-7a-5p == 0.0293925 (data field 1).
    assert m.values[0, 0] == pytest.approx(0.0293925)
    assert m.values[1, 1] == pytest.approx(0.01)
    assert m.values[0, 1] == NOT_TESTED


def test_indices(tmp_path):
    p = tmp_path / "GO_MF_validated_miRTarBase_all.csv"
    _write_matrix(p, ["A", "B"], [["mir-x", "1.0", "0.2"]])
    m = read_v2_matrix(p)
    assert m.pathway_index() == {"A": 0, "B": 1}
    assert m.mirna_index() == {"mir-x": 0}


def test_empty_file(tmp_path):
    p = tmp_path / "KEGG_validated_miRTarBase_strong.csv"
    p.write_bytes(b"")
    m = read_v2_matrix(p)
    assert m.pathways == []
    assert m.mirnas == []
    assert m.values.shape == (0, 0)


def test_misshaped_row_raises(tmp_path):
    p = tmp_path / "GO_BP_validated_miRTarBase_strong.csv"
    _write_matrix(p, ["A", "B", "C"], [["mir-x", "0.1", "0.2"]])  # 2 p for 3 cols
    with pytest.raises(ValueError, match="mis-shaped"):
        read_v2_matrix(p)


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_v2_matrix(tmp_path / "nope.csv")


def test_path_helpers():
    root = "/some/v2_root"
    d = v2_matrix_dir(root, "hsa")
    assert d.name == "hsa"
    assert "p_value_matrices" in str(d)
    p = v2_matrix_path(root, "mmu", "GO_MF", "validated_miRTarBase_strong")
    assert p.name == "GO_MF_validated_miRTarBase_strong.csv"
    with pytest.raises(ValueError):
        v2_matrix_dir(root, "rno")


def test_all_not_tested_row(tmp_path):
    p = tmp_path / "GO_CC_validated_miRTarBase_all.csv"
    _write_matrix(p, ["A", "B"], [["mir-x", "1.0", "1.0"]])
    m = read_v2_matrix(p)
    assert np.all(m.values == NOT_TESTED)
