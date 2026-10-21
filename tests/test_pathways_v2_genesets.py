"""Tests for v2 gene-level TERM2GENE reconstruction (pathways/v2_genesets.py).

Covers the GeneTrail2-result tarball reconstruction with a synthetic tarball
that mimics the real ``<variant>/<MIMAT>/<FAMILY>.txt`` layout and column
schema, plus a real-data smoke test gated on the offline v2 mmu tarball being
present.
"""

from __future__ import annotations

import io
import os
import tarfile
from pathlib import Path

import pytest

from mirpathdb_prod.pathways import v2_genesets

# Real offline v2 ground-truth tarball (mmu is small enough to stream in CI).
# The miRPathDB 2.0 reference archive lives outside the repository; point
# MPD_V2_REFERENCE_DIR at its version_2 directory to enable these checks
# (they skip cleanly when it is absent).
_V2_REFERENCE_DIR = os.environ.get("MPD_V2_REFERENCE_DIR", "")
V2_MMU_GENETRAIL2 = (
    Path(_V2_REFERENCE_DIR) / "miRPathDB2_mmu_genetrail2_results.tar.gz"
    if _V2_REFERENCE_DIR
    else None
)


def _v2_ground_truth_available() -> bool:
    if V2_MMU_GENETRAIL2 is None:
        return False
    try:
        return V2_MMU_GENETRAIL2.exists()
    except OSError:  # unreadable mount/permission: same as absent
        return False

_GT2_HEADER = "#Name\tReference\tHits\tScore\tExpected Score\tP-value\tInfo\tRegulation_direction"


def _gt2_row(name: str, hits: int, pvalue: float, info: str) -> str:
    # Reference/Score/Expected Score are unused by the loader -> placeholders.
    return f"{name}\t1\t{hits}\t0\t0\t{pvalue}\t{info}\t1"


def _make_synthetic_tarball(path: Path, organism: str) -> None:
    marker = f"miRPathDB2_{organism}_genetrail2_results"

    def _add(tf: tarfile.TarFile, arc: str, text: str) -> None:
        data = text.encode("utf-8")
        info = tarfile.TarInfo(name=arc)
        info.size = len(data)
        tf.addfile(info, io.BytesIO(data))

    with tarfile.open(path, "w:gz") as tf:
        # variant=miRTarBase_strong, two miRNAs, family GO_MF
        _add(
            tf,
            f"prefix/{marker}/miRTarBase_strong/MIMAT0000001/GO_MF.txt",
            "\n".join(
                [
                    _GT2_HEADER,
                    _gt2_row("ATP binding", 2, 0.0146, "Gja1,Pdgfb"),
                    _gt2_row("DNA binding", 1, 0.5, "Gja1"),
                    _gt2_row("empty term", 0, 1.0, ""),  # no overlap -> ignored
                ]
            ),
        )
        _add(
            tf,
            f"prefix/{marker}/miRTarBase_strong/MIMAT0000002/GO_MF.txt",
            "\n".join(
                [
                    _GT2_HEADER,
                    _gt2_row("ATP binding", 1, 0.7, "Zfp36l1"),
                    _gt2_row("kinase activity", 1, 0.2, "Mapk1"),
                ]
            ),
        )
        # a different family in the same variant (must NOT leak into GO_MF)
        _add(
            tf,
            f"prefix/{marker}/miRTarBase_strong/MIMAT0000001/GO_BP.txt",
            "\n".join([_GT2_HEADER, _gt2_row("apoptosis", 1, 0.3, "Casp3")]),
        )
        # a different variant (must NOT leak into miRTarBase_strong)
        _add(
            tf,
            f"prefix/{marker}/miRTarBase_all/MIMAT0000001/GO_MF.txt",
            "\n".join([_GT2_HEADER, _gt2_row("ATP binding", 1, 0.9, "Xyz1")]),
        )


def test_reconstruct_term2gene_unions_across_mirnas(tmp_path) -> None:
    tb = tmp_path / "miRPathDB2_mmu_genetrail2_results.tar.gz"
    _make_synthetic_tarball(tb, "mmu")

    res = v2_genesets.load_v2_term2gene(
        tb, organism="mmu", variant=v2_genesets.VARIANT_STRONG, family="GO_MF"
    )

    # ATP binding gene set = union of overlaps across both miRNAs.
    assert res.term2gene["ATP binding"] == {"Gja1", "Pdgfb", "Zfp36l1"}
    assert res.term2gene["DNA binding"] == {"Gja1"}
    assert res.term2gene["kinase activity"] == {"Mapk1"}
    # zero-overlap term contributed no genes -> absent.
    assert "empty term" not in res.term2gene
    # the GO_BP family did not leak in.
    assert "apoptosis" not in res.term2gene
    # the miRTarBase_all variant did not leak in.
    assert "Xyz1" not in res.universe

    # universe = union of all genes seen for this (variant, family).
    assert res.universe == {"Gja1", "Pdgfb", "Zfp36l1", "Mapk1"}
    # term2name is identity on v2 display names.
    assert res.term2name["ATP binding"] == "ATP binding"
    # per-miRNA target sets recovered.
    assert res.per_mirna_genes["MIMAT0000001"] == {"Gja1", "Pdgfb"}
    assert res.per_mirna_genes["MIMAT0000002"] == {"Zfp36l1", "Mapk1"}


def test_reconstruct_rejects_unknown_family(tmp_path) -> None:
    tb = tmp_path / "miRPathDB2_mmu_genetrail2_results.tar.gz"
    _make_synthetic_tarball(tb, "mmu")
    with pytest.raises(ValueError, match="unknown family"):
        v2_genesets.load_v2_term2gene(
            tb, organism="mmu", variant=v2_genesets.VARIANT_STRONG, family="PANTHER"
        )


def test_reconstruct_missing_tarball(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        v2_genesets.load_v2_term2gene(
            tmp_path / "nope.tar.gz",
            organism="mmu",
            variant=v2_genesets.VARIANT_STRONG,
            family="GO_MF",
        )


@pytest.mark.skipif(
    not _v2_ground_truth_available(),
    reason="offline v2 mmu GeneTrail2 tarball not available",
)
def test_real_v2_mmu_strong_go_mf_has_known_gene() -> None:
    # Real ground truth: MIMAT0000513 (mmu) GO_MF strong includes 'ATP binding'
    # with overlap gene 'Gja1' (verified directly from the tarball).
    res = v2_genesets.load_v2_term2gene(
        V2_MMU_GENETRAIL2,
        organism="mmu",
        variant=v2_genesets.VARIANT_STRONG,
        family="GO_MF",
    )
    assert res.term2gene, "no terms reconstructed from real v2 mmu GO_MF strong"
    assert "ATP binding" in res.term2gene
    assert "Gja1" in res.term2gene["ATP binding"]
    assert "Gja1" in res.universe
    assert "MIMAT0000513" in res.per_mirna_genes
