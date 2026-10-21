"""Tests for the single-pass multi-cell TERM2GENE collector.

Builds a tiny synthetic GeneTrail2-shaped ``.tar.gz`` (the real layout:
``.../miRPathDB2_<org>_genetrail2_results/<variant>/<MIMAT>/<FAMILY>.txt`` with
the 8-column result body) and checks that the streaming multi-cell collector
returns exactly what the single-cell loader does.
"""

from __future__ import annotations

import io
import tarfile

import pytest

from mirpathdb_prod.pathways.v2_genesets import load_v2_term2gene
from mirpathdb_prod.reproduction.v2_term2gene_stream import load_v2_term2gene_multi

_HEADER = "#Name\tReference\tHits\tScore\tExpected Score\tP-value\tInfo\tRegulation_direction"


def _result_body(rows):
    """rows = list[(term_name, info_genes_csv, pvalue)] -> GeneTrail2 file body."""
    lines = [_HEADER]
    for name, info, p in rows:
        hits = len([g for g in info.split(",") if g])
        lines.append(f"{name}\t{{}}\t{hits}\t{hits}\t0.1\t{p}\t{info}\t1")
    return "\n".join(lines) + "\n"


def _build_tarball(path, organism, entries):
    """entries = {(variant, mimat, family): result_rows}."""
    marker = f"miRPathDB2_{organism}_genetrail2_results"
    with tarfile.open(path, "w:gz") as tf:
        for (variant, mimat, family), rows in entries.items():
            body = _result_body(rows).encode("utf-8")
            name = f"local/x/{marker}/{variant}/{mimat}/{family}.txt"
            info = tarfile.TarInfo(name=name)
            info.size = len(body)
            tf.addfile(info, io.BytesIO(body))


@pytest.fixture
def synthetic_tarball(tmp_path):
    p = tmp_path / "miRPathDB2_mmu_genetrail2_results.tar.gz"
    entries = {
        ("miRTarBase_strong", "MIMAT0001", "GO_MF"): [
            ("ATP binding", "Gja1,Myc", "0.05"),
            ("DNA binding", "Myc", "0.2"),
        ],
        ("miRTarBase_strong", "MIMAT0002", "GO_MF"): [
            ("ATP binding", "Casp3", "0.1"),
        ],
        ("miRTarBase_all", "MIMAT0001", "GO_MF"): [
            ("ATP binding", "Gja1,Myc,Foo", "0.03"),
        ],
        ("miRTarBase_strong", "MIMAT0001", "GO_BP"): [
            ("apoptosis", "Casp3", "0.01"),
        ],
        # A prediction variant + KEGG family that must NOT leak into results.
        ("prediction_union", "MIMAT0009", "GO_MF"): [("ATP binding", "Zzz", "0.5")],
        ("miRTarBase_strong", "MIMAT0001", "KEGG"): [("pathwayX", "Myc", "0.4")],
    }
    _build_tarball(p, "mmu", entries)
    return p


def test_multi_collects_requested_cells(synthetic_tarball):
    out = load_v2_term2gene_multi(
        synthetic_tarball,
        organism="mmu",
        variants=["miRTarBase_strong", "miRTarBase_all"],
        families=["GO_MF", "GO_BP"],
    )
    # All 4 requested (variant, family) combos present.
    assert set(out) == {
        ("miRTarBase_strong", "GO_MF"),
        ("miRTarBase_strong", "GO_BP"),
        ("miRTarBase_all", "GO_MF"),
        ("miRTarBase_all", "GO_BP"),
    }
    mf_strong = out[("miRTarBase_strong", "GO_MF")]
    # ATP binding unions Gja1,Myc (MIMAT0001) + Casp3 (MIMAT0002).
    assert mf_strong.term2gene["ATP binding"] == {"Gja1", "Myc", "Casp3"}
    assert mf_strong.term2gene["DNA binding"] == {"Myc"}
    assert mf_strong.universe == {"Gja1", "Myc", "Casp3"}
    assert set(mf_strong.per_mirna_genes) == {"MIMAT0001", "MIMAT0002"}
    # GO_BP only has the strong/apoptosis entry.
    bp_strong = out[("miRTarBase_strong", "GO_BP")]
    assert bp_strong.term2gene == {"apoptosis": {"Casp3"}}


def test_prediction_and_kegg_excluded_when_not_requested(synthetic_tarball):
    out = load_v2_term2gene_multi(
        synthetic_tarball,
        organism="mmu",
        variants=["miRTarBase_strong"],
        families=["GO_MF"],
    )
    # Only the one requested cell; prediction_union + KEGG never appear.
    assert set(out) == {("miRTarBase_strong", "GO_MF")}
    assert "Zzz" not in out[("miRTarBase_strong", "GO_MF")].universe


def test_matches_single_cell_loader(synthetic_tarball):
    multi = load_v2_term2gene_multi(
        synthetic_tarball,
        organism="mmu",
        variants=["miRTarBase_strong", "miRTarBase_all"],
        families=["GO_MF", "GO_BP"],
    )
    for variant in ("miRTarBase_strong", "miRTarBase_all"):
        for family in ("GO_MF", "GO_BP"):
            single = load_v2_term2gene(
                synthetic_tarball, organism="mmu", variant=variant, family=family
            )
            cell = multi[(variant, family)]
            assert single.term2gene == cell.term2gene
            assert single.universe == cell.universe
            assert single.per_mirna_genes == cell.per_mirna_genes


def test_unknown_family_raises(synthetic_tarball):
    with pytest.raises(ValueError, match="unknown families"):
        load_v2_term2gene_multi(
            synthetic_tarball, organism="mmu", variants=["miRTarBase_strong"], families=["BOGUS"]
        )


def test_missing_tarball_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_v2_term2gene_multi(
            tmp_path / "nope.tar.gz",
            organism="mmu",
            variants=["miRTarBase_strong"],
            families=["GO_MF"],
        )
