"""Authoritative significant-subset long-form (Deliverable B) writer.

The production ORA kernel (:func:`mirpathdb_prod.enrichment.clusterprofiler.run_ora_batch`)
returns ``{mirna: {pathway_id: (p_value, p_adjust, overlap)}}`` per scored cell --
clusterProfiler computes the BH ``p.adjust`` and the ``Count`` overlap inline. The
matrix writer keeps only the raw ``p`` (untested cells become ``1.0``).

This module writes the significant-subset long-form contract
(``data/processed/enrichment_significant/<org>/<FAMILY>_<TIER>.tsv``) straight from
those authoritative tuples, so the published raw ``p_value`` + ``overlap`` are the
build's own values -- NOT re-derived downstream from the lossy matrix.

Adjustment convention (contract v2): the published ``p_adjust`` is the
**full-family** BH value -- BH over the ``M`` eligible terms of the cell (every term
of the complete TERM2GENE input with >= 1 gene in the universe; see
:mod:`mirpathdb_prod.enrichment.family_bh`), a family that is fixed by the inputs
before any query is seen. The enricher only lists the ``m`` terms overlapping a
query (including valid ``P = 1`` tests) and adjusts over those ``m``; that
query-dependent value is kept as the LEGACY column ``p_adjust_overlap``. Both are
computed from the complete listed family of each miRNA -- ``m_overlap`` is
``len(results[mirna])`` and must never be a retained-row count -- and the exact
migration identity ``q_full = min(1, q_overlap * M / m)`` is verified per miRNA
before anything is written. Retention (``FDR_CUTOFF``) applies to the full-family
value. ``q_full >= q_overlap`` always, so v2 retains a subset of the legacy rows.

When ``m_family`` is not given the writer keeps the legacy 5-column layout, where
``p_adjust`` is the overlap-only value (contract v1).
"""

from __future__ import annotations

import csv
from pathlib import Path

from mirpathdb_prod.enrichment.family_bh import adjust_cell_results

# The production significant-subset convention (ADR-0001): keep cells with
# BH-FDR <= this cutoff. Shared so every emitter/consumer uses one value.
FDR_CUTOFF = 0.25

# The production SIGNIFICANT criterion (ADR-0001): a retained cell is called
# significant at BH-FDR <= this value. Distinct from FDR_CUTOFF on purpose (cells
# are RETAINED at <= 0.25 and SIGNIFICANT at <= 0.05).
# ``scripts/export_frontend_contracts.py`` must agree -- pinned by tests.
FDR_SIGNIFICANT = 0.05

SIGNIFICANT_HEADER = ("mirna", "pathway_id", "p_value", "p_adjust", "overlap")
# Contract v2: the first five columns keep name/order/type (p_adjust = full-family
# value); the legacy overlap-only value and both family sizes follow.
SIGNIFICANT_HEADER_V2 = SIGNIFICANT_HEADER + ("p_adjust_overlap", "m_overlap", "m_family")


def write_significant_long_form(
    out_path: str | Path,
    mirnas: list[str],
    pvals: dict[str, dict[str, tuple[float, float, int]]],
    *,
    fdr_cutoff: float = FDR_CUTOFF,
    m_family: int | None = None,
) -> int:
    """Write the authoritative significant-subset long-form; return the row count.

    With ``m_family`` (the cell's eligible-family size ``M``) this emits the v2
    columns ``mirna, pathway_id, p_value, p_adjust, overlap, p_adjust_overlap,
    m_overlap, m_family`` for every scored cell whose full-family BH value
    ``p_adjust`` is ``<= fdr_cutoff``; ``p_adjust_overlap`` is the enricher's
    overlap-only value and ``m_overlap`` the complete number of terms the enricher
    listed for that miRNA. ``pvals`` must be the UNFILTERED ``run_ora_batch``
    output -- :func:`~mirpathdb_prod.enrichment.family_bh.adjust_cell_results`
    raises on a filtered family. Raw ``p_value`` and ``overlap`` are written
    verbatim (no ORA re-run).

    Without ``m_family`` the legacy 5-column layout is written, retaining on the
    overlap-only ``p_adjust``.

    Args:
        out_path: Destination TSV (parent dirs are created).
        mirnas: Row order (typically ``sorted(gene_sets)``).
        pvals: ``{mirna: {pathway_id: (p_value, p_adjust, overlap)}}`` from the ORA.
        fdr_cutoff: Keep cells with ``p_adjust <= fdr_cutoff``.
        m_family: Eligible-family size ``M`` of the cell; ``None`` for the legacy layout.

    Returns:
        Number of retained rows written.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    if m_family is None:
        with out_path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(SIGNIFICANT_HEADER)
            for mirna in mirnas:
                for pid, (p, padj, overlap) in pvals.get(mirna, {}).items():
                    if padj > fdr_cutoff:
                        continue
                    writer.writerow([mirna, pid, repr(float(p)), repr(float(padj)), int(overlap)])
                    n_rows += 1
        return n_rows

    adjusted = adjust_cell_results(pvals, m_family)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(SIGNIFICANT_HEADER_V2)
        for mirna in mirnas:
            rows = adjusted.get(mirna, {})
            m_overlap = len(rows)
            for pid, (p, q_overlap, q_full, overlap) in rows.items():
                if q_full > fdr_cutoff:
                    continue
                writer.writerow(
                    [
                        mirna,
                        pid,
                        repr(float(p)),
                        repr(float(q_full)),
                        int(overlap),
                        repr(float(q_overlap)),
                        m_overlap,
                        int(m_family),
                    ]
                )
                n_rows += 1
    return n_rows


def count_retained_rows(
    pvals: dict[str, dict[str, tuple[float, float, int]]],
    *,
    fdr_cutoff: float = FDR_CUTOFF,
    m_family: int | None = None,
) -> dict[str, int | None]:
    """Count the rows each convention retains at ``fdr_cutoff`` over ALL of ``pvals``.

    Returns ``rows`` (full-family retention when ``m_family`` is given, else the
    legacy count), ``rows_overlap_legacy`` (overlap-only ``p_adjust <= fdr_cutoff``),
    and ``m_overlap_min`` / ``m_overlap_max`` -- the complete listed-term count over
    the miRNAs with >= 1 listed term (``None`` when no miRNA was scored).
    """
    m_counts = [len(rows) for rows in pvals.values() if rows]
    legacy = sum(
        1 for rows in pvals.values() for (_p, padj, _ov) in rows.values() if padj <= fdr_cutoff
    )
    if m_family is None:
        full = legacy
    else:
        adjusted = adjust_cell_results(pvals, m_family)
        full = sum(
            1
            for rows in adjusted.values()
            for (_p, _q_overlap, q_full, _ov) in rows.values()
            if q_full <= fdr_cutoff
        )
    return {
        "rows": full,
        "rows_overlap_legacy": legacy,
        "m_overlap_min": min(m_counts) if m_counts else None,
        "m_overlap_max": max(m_counts) if m_counts else None,
    }


__all__ = [
    "FDR_CUTOFF",
    "FDR_SIGNIFICANT",
    "SIGNIFICANT_HEADER",
    "SIGNIFICANT_HEADER_V2",
    "count_retained_rows",
    "write_significant_long_form",
]
