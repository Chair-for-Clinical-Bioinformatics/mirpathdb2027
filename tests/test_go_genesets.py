"""Tests for the production GO gene-set source (org.*.eg.db -> TERM2GENE).

Pure-Python coverage: the long-form TSV parse/round-trip, the family/universe/
summary accessors, the org-db taxid mapping, and the GoGeneSets dataclass. The
R-backed materialisation (and the GO-symbol-space sanity) is exercised by an
integration test that runs only when the clusterProfiler conda env's Rscript +
the rat/macaque annotation packages are available.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from mirpathdb_prod.pathways import go_genesets as gg


def _rscript_available() -> str | None:
    """Return the Rscript binary to use, or None if no R is reachable."""
    return os.environ.get(gg.RSCRIPT_ENV) or shutil.which("Rscript")


def test_org_db_for_taxid() -> None:
    assert gg.org_db_for_taxid(10116) == "org.Rn.eg.db"
    assert gg.org_db_for_taxid(9544) == "org.Mmu.eg.db"
    assert gg.org_db_for_taxid(9606) == "org.Hs.eg.db"
    with pytest.raises(ValueError, match="no org.*.eg.db package"):
        gg.org_db_for_taxid(9031)  # chicken: excluded, no production org.db


def test_go_families_constant() -> None:
    assert gg.GO_FAMILIES == ("GO_BP", "GO_CC", "GO_MF")


def _toy_genesets() -> gg.GoGeneSets:
    return gg.GoGeneSets(
        taxid=10116,
        org_db_package="org.Rn.eg.db",
        term2gene={
            "GO_BP": {"GO:0001": {"Actb", "Gapdh"}, "GO:0002": {"Nppa"}},
            "GO_CC": {"GO:0003": {"Actb"}},
            "GO_MF": {},
        },
        term2name={"GO:0001": "process one", "GO:0002": "process two", "GO:0003": "comp"},
        provenance={"org_db_package": "org.Rn.eg.db", "org_db_version": "3.22.0"},
    )


def test_family_universe_summary() -> None:
    gs = _toy_genesets()
    assert gs.family("GO_BP") == {"GO:0001": {"Actb", "Gapdh"}, "GO:0002": {"Nppa"}}
    assert gs.universe("GO_BP") == {"Actb", "Gapdh", "Nppa"}
    assert gs.universe("GO_MF") == set()
    summ = gs.summary()
    assert summ["GO_BP"] == {"n_terms": 2, "n_genes": 3, "n_pairs": 3}
    assert summ["GO_CC"] == {"n_terms": 1, "n_genes": 1, "n_pairs": 1}


def test_family_unknown_raises() -> None:
    with pytest.raises(KeyError, match="no GO family"):
        _toy_genesets().family("GO_REACTOME")


def test_write_read_round_trip(tmp_path: Path) -> None:
    gs = _toy_genesets()
    out = tmp_path / "go.tsv"
    gg.write_go_genesets_tsv(gs, out)
    assert out.exists()
    back = gg.read_go_genesets_tsv(out, taxid=10116, org_db_package="org.Rn.eg.db")
    assert back.family("GO_BP") == gs.family("GO_BP")
    assert back.family("GO_CC") == gs.family("GO_CC")
    assert back.term2name["GO:0001"] == "process one"
    assert back.summary()["GO_BP"] == gs.summary()["GO_BP"]


# --------------------------------------------------------------------------- #
# integration (real org.*.eg.db via R; skipped if Rscript/packages absent)
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(_rscript_available() is None, reason="no Rscript on PATH/$MPD_RSCRIPT")
def test_integration_rat_go_genesets_symbol_space() -> None:
    """Build real rat GO sets and confirm they are non-trivial + symbol-spaced.

    Slow (shells to R + org.Rn.eg.db); guarded by Rscript presence. Confirms each
    GO family has thousands of terms, the universe is large, the provenance pins
    the org.db version, and the symbols are rat-cased (e.g. ``Actb`` not ``ACTB``).
    """
    try:
        gs = gg.build_go_genesets(10116)
    except RuntimeError as exc:  # org.Rn.eg.db not installed in this R
        pytest.skip(f"org.Rn.eg.db unavailable: {exc}")
    for fam in gg.GO_FAMILIES:
        assert len(gs.family(fam)) > 1000, f"{fam} suspiciously small"
    assert len(gs.universe("GO_BP")) > 10000
    assert gs.provenance.get("org_db_package") == "org.Rn.eg.db"
    assert gs.provenance.get("org_db_version")
    # rat symbols are title-cased (Actb), not human UPPERCASE (ACTB).
    all_syms = gs.universe("GO_BP")
    assert "Actb" in all_syms
    assert "ACTB" not in all_syms


# --------------------------------------------------------------------------- #
# Review remediation M2: GO ancestor propagation (GOALL) config flag.
# --------------------------------------------------------------------------- #
def test_go_propagation_columns_and_resolver() -> None:
    assert gg.GO_PROPAGATION_COLUMNS["all"] == ("GOALL", "ONTOLOGYALL")
    assert gg.GO_PROPAGATION_COLUMNS["direct"] == ("GO", "ONTOLOGY")
    assert gg.DEFAULT_GO_PROPAGATION == "all"
    # Default when unset -> canonical propagated GOALL.
    assert gg.resolve_go_propagation({}) == "all"
    assert gg.resolve_go_propagation({"go_propagation": "direct"}) == "direct"


def test_go_propagation_resolver_rejects_unknown() -> None:
    import pytest

    with pytest.raises(ValueError, match="unknown go_propagation"):
        gg.resolve_go_propagation({"go_propagation": "propagated"})


def test_go_driver_uses_column_variables_not_hardcoded() -> None:
    """The R driver indexes go[[go_col]] (parametrized) so it serves both modes;
    it no longer hardcodes the direct GO/ONTOLOGY columns."""
    d = gg._R_GO_DRIVER
    assert "go_col  <- args[[4]]" in d
    assert "ont_col <- args[[5]]" in d
    assert "go[[go_col]]" in d and "go[[ont_col]]" in d
