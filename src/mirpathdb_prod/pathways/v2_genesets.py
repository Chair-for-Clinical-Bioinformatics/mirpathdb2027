"""Reconstruct v2 gene-level TERM2GENE from the GeneTrail2 result trees.

The miRPathDB2 ``*.gmt`` files ship pathway -> **miRNA** membership, which is
the wrong space for driving ORA. The gene-level pathway membership v2 actually
tested against lives in the GeneTrail2 result tarballs:

    <variant>/<MIMAT_id>/<FAMILY>.txt

Each ``<FAMILY>.txt`` is a tab-separated GeneTrail2 over-representation result
for ONE miRNA against ONE family, with columns::

    #Name  Reference  Hits  Score  Expected Score  P-value  Info  Regulation_direction

- ``#Name``  = the term name (e.g. ``ATP binding``)
- ``Hits``   = overlap count (genes of this miRNA's target set in the term)
- ``Info``   = the overlap genes themselves, **comma-separated symbols**
               (e.g. ``Pdgfb,Zfp36l1``)
- ``P-value`` = the raw ORA p-value (the same number that lands in the v2
               p-value matrices)

Unioning the ``Info`` column across every miRNA of a variant reconstructs the
gene-set membership v2 used **within that variant's target universe** — exactly
the TERM2GENE needed to reproduce v2's ORA with clusterProfiler. The universe is
the union of all genes that appear in any term's ``Info`` for that variant.

This module does NOT run R or the network; it only reads (optionally streams)
the tarballs read-only.
"""

from __future__ import annotations

import io
import tarfile
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from mirpathdb_prod.pathways.gmt import KNOWN_FAMILIES

# GeneTrail2 variant subtree names inside the tarballs.
VARIANT_STRONG = "miRTarBase_strong"
VARIANT_ALL = "miRTarBase_all"
VARIANT_PRED_UNION = "prediction_union"
VARIANT_PRED_INTERSECTION = "prediction_intersection"

#: The two controlled (input-reproducible) validated variants.
CONTROLLED_VARIANTS: tuple[str, ...] = (VARIANT_STRONG, VARIANT_ALL)

# GeneTrail2 result column indices (0-based) in the <FAMILY>.txt files.
_COL_NAME = 0
_COL_HITS = 2
_COL_PVALUE = 5
_COL_INFO = 6


@dataclass
class V2TermGenes:
    """Reconstructed gene-level membership for one variant+family.

    Attributes:
        term2gene: ``term name -> set(gene_symbols)`` (the TERM2GENE for ORA).
        term2name: ``term name -> display name`` (identity here — v2 term names
            are already display strings; kept for the enricher's TERM2NAME slot).
        universe: Union of all genes appearing in any term (the background).
        per_mirna_genes: ``MIMAT_id -> set(gene_symbols)`` — each miRNA's target
            set as v2 saw it (union of that miRNA's per-term Info genes). Useful
            as the query gene sets when reproducing v2 cell-for-cell.
    """

    term2gene: dict[str, set[str]] = field(default_factory=dict)
    term2name: dict[str, str] = field(default_factory=dict)
    universe: set[str] = field(default_factory=set)
    per_mirna_genes: dict[str, set[str]] = field(default_factory=dict)


def _split_info_genes(info: str) -> list[str]:
    """Split a GeneTrail2 ``Info`` cell into gene symbols (comma-separated)."""
    info = info.strip()
    if not info:
        return []
    return [g.strip() for g in info.split(",") if g.strip()]


def _parse_result_text(text: str) -> Iterator[tuple[str, list[str]]]:
    """Yield ``(term_name, overlap_genes)`` from one GeneTrail2 result file body.

    Lines that are blank, the header (``#Name...``), or malformed are skipped.
    Only rows with at least one overlap gene contribute.
    """
    for raw in text.splitlines():
        line = raw.rstrip("\r")
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) <= _COL_INFO:
            continue
        name = cols[_COL_NAME].strip()
        if not name:
            continue
        genes = _split_info_genes(cols[_COL_INFO])
        if genes:
            yield name, genes


def _member_relpath_parts(
    member_name: str, prefix_depth_marker: str
) -> tuple[str, str, str] | None:
    """Parse ``.../<variant>/<MIMAT>/<FAMILY>.txt`` -> (variant, mimat, family).

    ``prefix_depth_marker`` is the directory name (e.g.
    ``miRPathDB2_mmu_genetrail2_results``) under which the variant subtrees sit;
    everything before it is ignored so the function is path-prefix agnostic.
    Returns ``None`` for members that are not a ``<FAMILY>.txt`` leaf.
    """
    parts = Path(member_name).parts
    if prefix_depth_marker not in parts:
        return None
    idx = parts.index(prefix_depth_marker)
    tail = parts[idx + 1 :]
    if len(tail) != 3:
        return None
    variant, mimat, leaf = tail
    if not leaf.endswith(".txt"):
        return None
    family = leaf[: -len(".txt")]
    if family not in KNOWN_FAMILIES:
        return None
    return variant, mimat, family


def load_v2_term2gene(
    genetrail2_tarball: str | Path,
    *,
    organism: str,
    variant: str,
    family: str,
) -> V2TermGenes:
    """Reconstruct the gene-level TERM2GENE for one (variant, family) from a tarball.

    Streams the tarball and reads only the ``<variant>/*/<family>.txt`` members,
    so the 72 GB hsa tarball is processed without bulk extraction. The mmu
    tarball (136 MB) streams just as well.

    Args:
        genetrail2_tarball: Path to ``miRPathDB2_<org>_genetrail2_results.tar.gz``.
        organism: ``"hsa"`` or ``"mmu"`` (selects the in-tar directory marker
            ``miRPathDB2_<org>_genetrail2_results``).
        variant: One of :data:`CONTROLLED_VARIANTS` (``miRTarBase_strong`` /
            ``miRTarBase_all``); the prediction variants are out of scope for the
            controlled reproduction but are accepted if present.
        family: One of :data:`mirpathdb_prod.pathways.gmt.KNOWN_FAMILIES`.

    Returns:
        A :class:`V2TermGenes` with ``term2gene``, ``term2name``, ``universe``,
        and ``per_mirna_genes``.

    Raises:
        FileNotFoundError: if the tarball is missing.
        ValueError: if ``family`` is not a known family.
    """
    tar_path = Path(genetrail2_tarball)
    if not tar_path.exists():
        raise FileNotFoundError(f"GeneTrail2 tarball not found: {tar_path}")
    if family not in KNOWN_FAMILIES:
        raise ValueError(f"unknown family {family!r}; known: {KNOWN_FAMILIES}")

    marker = f"miRPathDB2_{organism}_genetrail2_results"
    out = V2TermGenes()

    with tarfile.open(tar_path, mode="r:gz") as tf:
        for member in tf:
            if not member.isfile():
                continue
            parsed = _member_relpath_parts(member.name, marker)
            if parsed is None:
                continue
            m_variant, mimat, m_family = parsed
            if m_variant != variant or m_family != family:
                continue
            extracted = tf.extractfile(member)
            if extracted is None:
                continue
            body = io.TextIOWrapper(extracted, encoding="utf-8", errors="replace").read()
            mirna_genes = out.per_mirna_genes.setdefault(mimat, set())
            for term_name, genes in _parse_result_text(body):
                bucket = out.term2gene.setdefault(term_name, set())
                bucket.update(genes)
                out.term2name.setdefault(term_name, term_name)
                out.universe.update(genes)
                mirna_genes.update(genes)
    return out
