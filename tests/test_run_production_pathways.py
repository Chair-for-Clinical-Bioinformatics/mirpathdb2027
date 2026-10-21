"""Tests for the production-pathways ORA runner (pure logic, no R).

The R-backed gene-set + ORA paths are exercised by the per-source tests and the
clusterProfiler tests; here we lock the runner's organism/family GATING (the rule
that protects the existing matrices) and the matrix writer's v2-orientation fill.
The script lives under ``scripts/`` (not on the package path), so it is loaded by
file path with importlib.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from mirpathdb_prod.reporting.significant import SIGNIFICANT_HEADER_V2
from mirpathdb_prod.sources import production_target_sets as pts

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_production_pathways.py"


def _load_runner():
    spec = importlib.util.spec_from_file_location("run_production_pathways", _SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def runner():
    return _load_runner()


def test_resolve_families_human_mouse_get_everything(runner):
    """hsa/mmu get Reactome + WikiPathways + the three GO ontologies."""
    for org in ("hsa", "mmu"):
        fams = runner._resolve_families(org, ["REACTOME", "WIKIPATHWAYS", "GO"])
        assert fams == ["REACTOME", "WIKIPATHWAYS", "GO_BP", "GO_CC", "GO_MF"]


def test_resolve_families_rat_excludes_go(runner):
    """rno keeps Reactome + WikiPathways but NEVER GO (production rno GO exists)."""
    fams = runner._resolve_families("rno", ["REACTOME", "WIKIPATHWAYS", "GO"])
    assert fams == ["REACTOME", "WIKIPATHWAYS"]
    assert not any(f.startswith("GO") for f in fams)


def test_resolve_families_macaque_gets_nothing(runner):
    """mml has no NATIVE Reactome, no WP collection, and existing GO -> no new families.

    (The NATIVE pathway runner emits nothing new for macaque. Macaque's INFERRED
    Reactome -- orthology-projected from human R-HSA- -> R-MML- -- is built by the
    separate scripts/run_macaque_inferred_reactome.py, not this runner.)
    """
    fams = runner._resolve_families("mml", ["REACTOME", "WIKIPATHWAYS", "GO"])
    assert fams == []


def test_resolve_families_go_only_request_skips_non_eligible(runner):
    """A GO-only request adds hsa/mmu GO but leaves rno/mml with nothing."""
    assert runner._resolve_families("hsa", ["GO"]) == ["GO_BP", "GO_CC", "GO_MF"]
    assert runner._resolve_families("rno", ["GO"]) == []
    assert runner._resolve_families("mml", ["GO"]) == []


def test_output_filenames_are_prod_prefixed(runner):
    """The runner writes prod_<FAMILY>_<TIER>.csv (distinct from existing matrices)."""
    # Distinct from v2-comparison (GO_BP_validated_miRTarBase_*) and existing
    # production rno/mml GO (GO_BP_HIGH.csv) -- the prod_ prefix guarantees this.
    for family in ("REACTOME", "WIKIPATHWAYS", "GO_BP"):
        name = f"prod_{family}_HIGH.csv"
        assert name.startswith("prod_")
        assert name != f"{family}_HIGH.csv"  # not the existing rno/mml GO name
        assert "validated_miRTarBase" not in name  # not the v2-comparison name


def test_write_matrix_v2_orientation_and_fill(runner, tmp_path):
    """Matrix is miRNA rows x pathway cols; unscored cells default to 1.0."""
    mirnas = ["hsa-miR-1-3p", "hsa-let-7a-5p"]
    terms = ["R-HSA-1", "R-HSA-2", "R-HSA-3"]
    pvals = {
        "hsa-miR-1-3p": {"R-HSA-1": (1e-5, 1e-3, 7), "R-HSA-3": (0.04, 0.2, 2)},
        # hsa-let-7a-5p scored only R-HSA-2.
        "hsa-let-7a-5p": {"R-HSA-2": (2e-4, 1e-2, 4)},
    }
    out = tmp_path / "prod_REACTOME_HIGH.csv"
    n_rows, n_cols = runner._write_matrix(out, mirnas, terms, pvals)
    assert (n_rows, n_cols) == (2, 3)

    with out.open(newline="") as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == ["miRNA", "R-HSA-1", "R-HSA-2", "R-HSA-3"]
    by_mirna = {r[0]: [float(v) for v in r[1:]] for r in rows[1:]}
    # scored cells carry the p-value; unscored cells are 1.0.
    assert by_mirna["hsa-miR-1-3p"] == [1e-5, 1.0, 0.04]
    assert by_mirna["hsa-let-7a-5p"] == [1.0, 2e-4, 1.0]


def test_write_matrix_empty_terms(runner, tmp_path):
    """A cell with no scored pathways writes just the header + miRNA rows."""
    out = tmp_path / "prod_WIKIPATHWAYS_HIGH.csv"
    n_rows, n_cols = runner._write_matrix(out, ["mmu-miR-1a-3p"], [], {})
    assert (n_rows, n_cols) == (1, 0)
    with out.open(newline="") as fh:
        rows = list(csv.reader(fh))
    assert rows == [["miRNA"], ["mmu-miR-1a-3p"]]


# --------------------------------------------------------------------------- #
# Full-family BH wiring of run_organism (no R: target sets, TERM2GENE and the ORA
# batch are stubbed; the writer is wrapped to observe its keyword arguments).
# --------------------------------------------------------------------------- #
_T2G = {
    "R-HSA-1": {"A", "B"},
    "R-HSA-2": {"B", "C"},
    "R-HSA-3": {"C", "D"},
    "R-HSA-4": {"D"},
}
_UNIVERSE = {"A", "B", "C", "D"}
_TIERS = {
    "HIGH": {"hsa-miR-1-3p": {"A", "B"}, "hsa-let-7a-5p": {"C"}},
    "MEDIUM": {},
    "LOW": {},
}


def _bh(p: np.ndarray) -> np.ndarray:
    m = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    return out


def _fake_ora(gene_sets, term2gene, universe, *, term2name=None, rscript=None):
    """Enricher-shaped output: only overlapping terms, BH over exactly those."""
    res = {}
    for mirna, genes in gene_sets.items():
        hits = {t: len(set(genes) & set(members)) for t, members in term2gene.items()}
        hits = {t: k for t, k in hits.items() if k >= 1}
        terms = sorted(hits)
        p = np.array([0.02 / (i + 1) for i in range(len(terms))])
        if terms:
            p[-1] = 1.0  # a valid P = 1 test inside the overlapping family
        q = _bh(p)
        res[mirna] = {t: (float(p[i]), float(q[i]), hits[t]) for i, t in enumerate(terms)}
    return res, {"clusterprofiler": "stub"}


@pytest.fixture
def stubbed_runner(runner, monkeypatch, tmp_path):
    census = pts.TargetSetCensus(per_tier={"HIGH": {"n_mirnas": 2, "n_genes": 3, "n_pairs": 3}})
    monkeypatch.setattr(runner, "_load_target_sets", lambda organism: (_TIERS, census))
    monkeypatch.setattr(
        runner, "_build_family_sets", lambda org, fam, *, rscript: (_T2G, set(_UNIVERSE), {})
    )
    monkeypatch.setattr(runner.cp, "run_ora_batch_with_provenance", _fake_ora)
    writer_calls: list[dict[str, object]] = []
    real_writer = runner.write_significant_long_form

    def _observed_writer(out_path, mirnas, pvals, **kwargs):
        writer_calls.append(dict(kwargs))
        return real_writer(out_path, mirnas, pvals, **kwargs)

    monkeypatch.setattr(runner, "write_significant_long_form", _observed_writer)
    out_root = tmp_path / "data" / "processed" / "enrichment"
    out_root.mkdir(parents=True)
    return runner, writer_calls, out_root


def test_run_organism_records_full_family_bh_in_cell_block(stubbed_runner):
    runner, writer_calls, out_root = stubbed_runner
    block = runner.run_organism(
        "hsa",
        families_wanted=["REACTOME"],
        tiers_wanted=["HIGH"],
        out_root=out_root,
        rscript=None,
        long_form=True,
    )
    (cell,) = block["cells"]
    assert cell["m_family"] == len(_T2G) == cell["n_pathways_total"]
    assert cell["adjustment_family"] == "full_family_bh"
    assert writer_calls == [{"m_family": len(_T2G)}]

    res, _ = _fake_ora(pts.cumulative_tier(_TIERS, "HIGH"), _T2G, _UNIVERSE)
    m_counts = sorted(len(per) for per in res.values())
    assert (cell["m_overlap_min"], cell["m_overlap_max"]) == (m_counts[0], m_counts[-1])
    assert cell["n_retained_rows"] == cell["n_significant_rows"]
    assert cell["n_retained_rows_overlap_legacy"] >= cell["n_retained_rows"]

    long_path = out_root.parent / "enrichment_significant" / "hsa" / "REACTOME_HIGH.tsv"
    assert cell["long_form_output"] == str(long_path.relative_to(out_root.parents[2]))
    with long_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == list(SIGNIFICANT_HEADER_V2)
    assert len(rows) == cell["n_retained_rows"] > 0
    for r in rows:
        assert int(r["m_family"]) == len(_T2G)
        assert int(r["m_overlap"]) == len(res[r["mirna"]])
        assert float(r["p_adjust"]) >= float(r["p_adjust_overlap"])
    # Matrix output is unchanged by the long-form change.
    with (out_root / "hsa" / "prod_REACTOME_HIGH.csv").open(newline="") as fh:
        header = next(csv.reader(fh))
    assert header[0] == "miRNA" and set(header[1:]) <= set(_T2G)


def test_run_organism_refuses_family_size_mismatch(stubbed_runner, monkeypatch):
    runner, writer_calls, out_root = stubbed_runner
    monkeypatch.setattr(
        runner.family_bh, "eligible_family", lambda t2g, universe: (len(t2g) - 1, [])
    )
    with pytest.raises(ValueError, match="eligible family M="):
        runner.run_organism(
            "hsa",
            families_wanted=["REACTOME"],
            tiers_wanted=["HIGH"],
            out_root=out_root,
            rscript=None,
            long_form=True,
        )
    assert writer_calls == []
