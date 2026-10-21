"""Tests for the production Reactome gene-set source (reactome.db -> TERM2GENE).

Pure-Python coverage: the long-form TSV parse/round-trip, the family/universe/
summary accessors, the species-code taxid mapping (incl. macaque exclusion), and
the ReactomeGeneSets dataclass. The R-backed materialisation (reactome.db x
org.*.eg.db, and the symbol-space sanity) is exercised by an integration test
that runs only when the clusterProfiler conda env's Rscript + reactome.db +
the organism annotation packages are available.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from mirpathdb_prod.pathways import reactome_genesets as rx


def _rscript_available() -> str | None:
    return os.environ.get(rx.RSCRIPT_ENV) or shutil.which("Rscript")


def test_species_code_for_taxid() -> None:
    assert rx.species_code_for_taxid(9606) == "HSA"
    assert rx.species_code_for_taxid(10090) == "MMU"
    assert rx.species_code_for_taxid(10116) == "RNO"
    # macaque (9544): Reactome ships ZERO curated R-MML- pathways -> excluded from
    # the NATIVE map (inferred macaque Reactome is built separately).
    with pytest.raises(ValueError, match="no curated R-MML-|taxid 9544"):
        rx.species_code_for_taxid(9544)


def test_reactome_family_constant() -> None:
    assert rx.REACTOME_FAMILY == "REACTOME"
    # macaque is intentionally absent from the NATIVE species map (no curated
    # R-MML-); inferred macaque Reactome is built by inferred_reactome_genesets.
    assert 9544 not in rx.REACTOME_SPECIES_CODE
    assert set(rx.REACTOME_SPECIES_CODE) == {9606, 10090, 10116}


def _toy_genesets() -> rx.ReactomeGeneSets:
    return rx.ReactomeGeneSets(
        taxid=9606,
        species_code="HSA",
        org_db_package="org.Hs.eg.db",
        term2gene={
            "R-HSA-1": {"ACTB", "GAPDH"},
            "R-HSA-2": {"TP53"},
        },
        term2name={"R-HSA-1": "pathway one", "R-HSA-2": "pathway two"},
        provenance={"reactome_db_version": "1.95.0", "org_db_version": "3.22.0"},
    )


def test_family_universe_summary() -> None:
    gs = _toy_genesets()
    assert gs.family() == {"R-HSA-1": {"ACTB", "GAPDH"}, "R-HSA-2": {"TP53"}}
    assert gs.universe() == {"ACTB", "GAPDH", "TP53"}
    assert gs.summary() == {"n_terms": 2, "n_genes": 3, "n_pairs": 3}


def test_family_unknown_raises() -> None:
    with pytest.raises(KeyError, match="no Reactome family"):
        _toy_genesets().family("GO_BP")


def test_write_read_round_trip(tmp_path: Path) -> None:
    gs = _toy_genesets()
    out = tmp_path / "reactome.tsv"
    rx.write_reactome_genesets_tsv(gs, out)
    assert out.exists()
    back = rx.read_reactome_genesets_tsv(
        out, taxid=9606, species_code="HSA", org_db_package="org.Hs.eg.db"
    )
    assert back.family() == gs.family()
    assert back.term2name["R-HSA-1"] == "pathway one"
    assert back.summary() == gs.summary()


# --------------------------------------------------------------------------- #
# integration (real reactome.db x org.*.eg.db via R; skipped if R/pkgs absent)
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(_rscript_available() is None, reason="no Rscript on PATH/$MPD_RSCRIPT")
def test_integration_mouse_reactome_symbol_space() -> None:
    """Build real mouse Reactome sets and confirm non-trivial + symbol-spaced.

    Slow (shells to R + reactome.db + org.Mm.eg.db); guarded by Rscript presence.
    Mouse is used (smaller than human) to keep the test quick. Confirms hundreds
    of pathways, a large symbol universe, version pins, and rat/mouse-cased
    symbols (``Actb`` not ``ACTB``).
    """
    try:
        gs = rx.build_reactome_genesets(10090)
    except RuntimeError as exc:  # reactome.db / org.Mm.eg.db not installed
        pytest.skip(f"reactome.db/org.Mm.eg.db unavailable: {exc}")
    assert len(gs.family()) > 500, "mouse Reactome suspiciously small"
    assert len(gs.universe()) > 3000
    assert gs.provenance.get("reactome_db_version")
    assert gs.species_code == "MMU"
    uni = gs.universe()
    assert "Actb" in uni
    assert "ACTB" not in uni  # mouse symbols are title-cased, not human UPPERCASE
    # pathway ids carry the mouse stable-id prefix.
    assert all(pid.startswith("R-MMU-") for pid in gs.term2gene)
