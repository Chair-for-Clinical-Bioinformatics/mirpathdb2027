"""Tests for the controlled-reproduction runner (matrix assembly + universe).

The clusterProfiler R call is monkeypatched so the runner's assembly logic
(pathway/miRNA restriction, universe policy, v2-orientation matrix fill) is
tested without R. One R-backed end-to-end test is included but skipped when no
clusterProfiler Rscript is available.
"""

from __future__ import annotations

import shutil
import subprocess

import numpy as np
import pytest

from mirpathdb_prod.enrichment import clusterprofiler as cp
from mirpathdb_prod.pathways.v2_genesets import V2TermGenes
from mirpathdb_prod.reproduction import runner
from mirpathdb_prod.reproduction.v2_matrix import NOT_TESTED, V2Matrix


def _v2_matrix(mirnas, pathways, values):
    return V2Matrix(pathways=list(pathways), mirnas=list(mirnas), values=np.asarray(values))


def _v2_t2g(term2gene):
    universe = set()
    for g in term2gene.values():
        universe |= set(g)
    return V2TermGenes(
        term2gene={t: set(g) for t, g in term2gene.items()},
        term2name={t: t for t in term2gene},
        universe=universe,
        per_mirna_genes={},
    )


def test_run_cell_assembles_v2_orientation(monkeypatch):
    # 2 miRNAs x 3 pathways v2 matrix; we have a target set for both miRNAs.
    vm = _v2_matrix(
        ["mir-A", "mir-B"],
        ["ATP binding", "DNA binding", "RNA binding"],
        [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
    )
    t2g = _v2_t2g(
        {"ATP binding": {"G1", "G2"}, "DNA binding": {"G2"}, "RNA binding": {"G3"}}
    )
    target_sets = {"mir-A": {"G1", "G2"}, "mir-B": {"G3"}}

    captured = {}

    def fake_batch(gene_sets, t2g_in, universe, *, term2name=None, rscript=None):
        captured["gene_sets"] = {k: set(v) for k, v in gene_sets.items()}
        captured["t2g"] = {k: set(v) for k, v in t2g_in.items()}
        captured["universe"] = set(universe)
        captured["term2name"] = term2name
        # Return a plausible per-set result keyed by display name (TERM2NAME path).
        return (
            {
                "mir-A": {"ATP binding": (0.01, 0.03, 2), "DNA binding": (0.5, 0.5, 1)},
                "mir-B": {"RNA binding": (0.2, 0.2, 1)},
            },
            {"clusterprofiler_version": "X", "r_version": "Y"},
        )

    monkeypatch.setattr(cp, "_run_ora_batch_full", fake_batch)

    res = runner.run_cell(
        organism="hsa",
        variant="validated_miRTarBase_strong",
        family="GO_MF",
        v2_matrix=vm,
        v2_term2gene=t2g,
        target_sets=target_sets,
        universe_choice="mirtarbase_union",
    )
    # Matrix is in v2 orientation, full label space, 1.0 where unscored.
    assert res.values.shape == (2, 3)
    assert res.values[0, 0] == pytest.approx(0.01)  # mir-A / ATP binding
    assert res.values[0, 1] == pytest.approx(0.5)  # mir-A / DNA binding
    assert res.values[0, 2] == NOT_TESTED  # mir-A / RNA binding (unscored)
    assert res.values[1, 2] == pytest.approx(0.2)  # mir-B / RNA binding
    assert res.values[1, 0] == NOT_TESTED
    # mirtarbase_union universe = union of target genes (plus all t2g genes by
    # the defensive top-up only for term2gene_union; here it is the target union).
    assert captured["universe"] == {"G1", "G2", "G3"}
    assert res.provenance["clusterprofiler_version"] == "X"


def test_run_cell_term2gene_union_universe(monkeypatch):
    vm = _v2_matrix(["mir-A"], ["ATP binding"], [[1.0]])
    t2g = _v2_t2g({"ATP binding": {"G1", "G2", "G9"}})
    target_sets = {"mir-A": {"G1"}}

    captured = {}

    def fake_batch(gene_sets, t2g_in, universe, *, term2name=None, rscript=None):
        captured["universe"] = set(universe)
        return ({"mir-A": {}}, {})

    monkeypatch.setattr(cp, "_run_ora_batch_full", fake_batch)
    res = runner.run_cell(
        organism="mmu",
        variant="validated_miRTarBase_strong",
        family="GO_MF",
        v2_matrix=vm,
        v2_term2gene=t2g,
        target_sets=target_sets,
        universe_choice="term2gene_union",
    )
    # term2gene_union universe = all genes in the TERM2GENE (not the query).
    assert captured["universe"] == {"G1", "G2", "G9"}
    assert res.n_universe == 3


def test_run_cell_restricts_to_matrix_rows_and_cols(monkeypatch):
    # v2 matrix has a pathway with no t2g and a miRNA with no target set.
    vm = _v2_matrix(
        ["mir-A", "mir-MISSING"],
        ["ATP binding", "PATH_NO_T2G"],
        [[1.0, 1.0], [1.0, 1.0]],
    )
    t2g = _v2_t2g({"ATP binding": {"G1"}})  # PATH_NO_T2G absent
    target_sets = {"mir-A": {"G1"}}  # mir-MISSING absent

    captured = {}

    def fake_batch(gene_sets, t2g_in, universe, *, term2name=None, rscript=None):
        captured["gene_sets"] = dict(gene_sets)
        captured["t2g_terms"] = set(t2g_in)
        captured["term2name"] = dict(term2name) if term2name else {}
        return ({"mir-A": {"ATP binding": (0.02, 0.02, 1)}}, {})

    monkeypatch.setattr(cp, "_run_ora_batch_full", fake_batch)
    res = runner.run_cell(
        organism="hsa",
        variant="validated_miRTarBase_strong",
        family="GO_MF",
        v2_matrix=vm,
        v2_term2gene=t2g,
        target_sets=target_sets,
        universe_choice="mirtarbase_union",
    )
    # Only mir-A is submitted (mir-MISSING has no target set).
    assert set(captured["gene_sets"]) == {"mir-A"}
    # Only matched pathways drive the TERM2GENE / TERM2NAME.
    assert captured["t2g_terms"] == {"ATP binding"}
    assert set(captured["term2name"]) == {"ATP binding"}
    assert res.n_terms_matched == 1
    # But the OUTPUT matrix still spans the full v2 label space (rows+cols).
    assert res.values.shape == (2, 2)
    assert res.mirnas == ["mir-A", "mir-MISSING"]
    assert res.values[1, 0] == NOT_TESTED  # mir-MISSING unscored
    assert res.values[0, 1] == NOT_TESTED  # PATH_NO_T2G unscored


def test_write_cell_matrix_roundtrip(monkeypatch, tmp_path):
    from mirpathdb_prod.reproduction.v2_matrix import read_v2_matrix

    vm = _v2_matrix(["mir-A"], ["ATP binding", "DNA binding"], [[1.0, 1.0]])
    t2g = _v2_t2g({"ATP binding": {"G1"}, "DNA binding": {"G1"}})
    target_sets = {"mir-A": {"G1"}}
    monkeypatch.setattr(
        cp,
        "_run_ora_batch_full",
        lambda *a, **k: ({"mir-A": {"ATP binding": (0.01, 0.01, 1)}}, {}),
    )
    res = runner.run_cell(
        organism="hsa",
        variant="validated_miRTarBase_strong",
        family="GO_MF",
        v2_matrix=vm,
        v2_term2gene=t2g,
        target_sets=target_sets,
    )
    out = tmp_path / "GO_MF_validated_miRTarBase_strong.csv"
    runner.write_cell_matrix(res, out)
    # Re-read: same orientation/shape, value preserved.
    back = read_v2_matrix(out)
    assert back.pathways == ["ATP binding", "DNA binding"]
    assert back.mirnas == ["mir-A"]
    assert back.values[0, 0] == pytest.approx(0.01)
    assert back.values[0, 1] == NOT_TESTED


def test_variant_token_mapping():
    # The matrix token -> GeneTrail2 subtree token mapping is the verified one.
    assert runner.VARIANTS == {
        "validated_miRTarBase_strong": "miRTarBase_strong",
        "validated_miRTarBase_all": "miRTarBase_all",
    }


# --- R-backed end-to-end (skipped without clusterProfiler) -------------------


def _find_rscript():
    import os

    # $MPD_RSCRIPT (point it at the conda env that has clusterProfiler), else
    # whatever Rscript is on PATH. No machine-specific path is baked in.
    override = os.environ.get("MPD_RSCRIPT")
    if override and shutil.which(override):
        return override
    return shutil.which("Rscript")


def _r_ok():
    rscript = _find_rscript()
    if not rscript:
        return False, None
    try:
        proc = subprocess.run(
            [rscript, "-e", 'library(clusterProfiler); cat("ok")'],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        return ("ok" in proc.stdout), rscript
    except (OSError, subprocess.SubprocessError):
        return False, rscript


_R_OK, _RSCRIPT = _r_ok()
_skip_no_r = pytest.mark.skipif(not _R_OK, reason="Rscript with clusterProfiler not available")


@_skip_no_r
def test_run_cell_end_to_end_nondegenerate():
    vm = _v2_matrix(
        ["mir-A", "mir-B"],
        ["ATP binding", "DNA binding", "RNA binding"],
        [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
    )
    t2g = _v2_t2g(
        {
            "ATP binding": {"G1", "G2", "G3", "G4"},
            "DNA binding": {"G2", "G5", "G6"},
            "RNA binding": {"G7", "G8"},
        }
    )
    target_sets = {"mir-A": {"G1", "G2", "G3"}, "mir-B": {"G7", "G8"}}
    res = runner.run_cell(
        organism="hsa",
        variant="validated_miRTarBase_strong",
        family="GO_MF",
        v2_matrix=vm,
        v2_term2gene=t2g,
        target_sets=target_sets,
        universe_choice="mirtarbase_union",
        rscript=_RSCRIPT,
    )
    # mir-A is strongly enriched for ATP binding (3/3 of its genes are members).
    i = res.mirnas.index("mir-A")
    j = res.pathways.index("ATP binding")
    assert 0.0 < res.values[i, j] < 1.0
    assert res.provenance["clusterprofiler_version"]
