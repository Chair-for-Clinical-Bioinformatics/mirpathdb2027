"""End-to-end controlled v2-reproduction runner (the workhorse glue).

For each comparison **cell** = (organism, matrix-variant, family) this:

1. reads the frozen v2 matrix to learn which miRNAs (rows) and pathways
   (columns) must line up for a cell-for-cell comparison;
2. takes the reconstructed v2 TERM2GENE (pathway -> gene symbols) for the cell,
   restricted to the pathways present in the v2 matrix;
3. builds each miRNA's target gene set from miRTarBase (gene-SYMBOL space — the
   v2 GeneTrail2 ID space), restricted to the matrix's miRNAs;
4. chooses the ORA universe (documented; primary + a sensitivity variant);
5. runs :func:`mirpathdb_prod.enrichment.clusterprofiler.run_ora_batch` over ALL
   the cell's miRNAs in ONE parallel R session;
6. assembles our ``(miRNA x pathway)`` p-value matrix in the SAME orientation as
   v2 (miRNA rows, pathway columns; not-tested cells filled with ``1.0`` to match
   v2's sentinel) and writes it as CSV.

Variant token mapping (verified): the v2 **matrix** files are named with
``validated_miRTarBase_strong`` / ``validated_miRTarBase_all``; the GeneTrail2
**subtree** uses ``miRTarBase_strong`` / ``miRTarBase_all``; and
:func:`mirpathdb_prod.sources.target_sets.build_target_sets` expects the matrix
token. :data:`VARIANTS` carries the three aligned spellings per variant.

Universe (the #1 divergence driver — documented choice + sensitivity):
    * **primary** ``term2gene_union`` — the union of all genes appearing in the
      cell's reconstructed TERM2GENE. This is the most v2-faithful background:
      clusterProfiler's effective universe for a term is ``universe ∩
      genes(TERM2GENE)``, so this exactly mirrors the gene space GeneTrail2 had
      available for that variant.
    * **sensitivity** ``mirtarbase_union`` — the union of all target genes across
      the cell's miRNAs (the "all genes any miRNA could hit" background). Because
      the effective N is intersected with the TERM2GENE space, this differs from
      the primary only through genes that are in some miRNA's targets but in no
      term (and vice-versa); the validator quantifies the resulting Δ.
"""

from __future__ import annotations

import csv
import time
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np

from mirpathdb_prod.enrichment import clusterprofiler
from mirpathdb_prod.pathways.v2_genesets import V2TermGenes
from mirpathdb_prod.reproduction.v2_matrix import NOT_TESTED, V2Matrix

#: (matrix_token, genetrail2_token) for each controlled variant.
VARIANTS: dict[str, str] = {
    "validated_miRTarBase_strong": "miRTarBase_strong",
    "validated_miRTarBase_all": "miRTarBase_all",
}

UniverseChoice = Literal["term2gene_union", "mirtarbase_union"]
DEFAULT_UNIVERSE: UniverseChoice = "term2gene_union"


@dataclass
class CellResult:
    """Outcome of running ONE (organism, variant, family) cell.

    Attributes:
        organism / variant / family: The cell identity.
        universe_choice: Which background was used.
        mirnas / pathways: The matched row/column labels (v2 order).
        values: Our ``(n_mirnas, n_pathways)`` p-value matrix (``1.0`` ==
            not tested, matching v2's sentinel).
        n_universe: Size of the universe handed to enricher.
        n_terms_t2g: Number of terms in the reconstructed TERM2GENE for the cell.
        n_terms_matched: Pathways from the v2 matrix that exist in the TERM2GENE.
        elapsed_s: Wall-time for the enricher batch call.
        provenance: clusterProfiler/R version pins from the run.
    """

    organism: str
    variant: str
    family: str
    universe_choice: str
    mirnas: list[str]
    pathways: list[str]
    values: np.ndarray
    n_universe: int
    n_terms_t2g: int
    n_terms_matched: int
    elapsed_s: float
    provenance: dict[str, str] = field(default_factory=dict)


def _restrict_term2gene(
    term2gene: Mapping[str, set[str]], pathways: Iterable[str]
) -> dict[str, set[str]]:
    """Keep only the TERM2GENE terms that are columns of the v2 matrix."""
    want = set(pathways)
    return {t: set(g) for t, g in term2gene.items() if t in want}


def _universe_for(
    choice: UniverseChoice,
    t2g_genes: set[str],
    target_sets: Mapping[str, set[str]],
) -> set[str]:
    """Return the background gene universe for the chosen policy."""
    if choice == "term2gene_union":
        return set(t2g_genes)
    if choice == "mirtarbase_union":
        genes: set[str] = set()
        for s in target_sets.values():
            genes |= s
        return genes
    raise ValueError(f"unknown universe choice {choice!r}")


def run_cell(
    *,
    organism: str,
    variant: str,
    family: str,
    v2_matrix: V2Matrix,
    v2_term2gene: V2TermGenes,
    target_sets: Mapping[str, set[str]],
    universe_choice: UniverseChoice = DEFAULT_UNIVERSE,
    rscript: str | None = None,
) -> CellResult:
    """Run one controlled (organism, variant, family) cell and assemble its matrix.

    Args:
        organism: ``"hsa"`` / ``"mmu"``.
        variant: A matrix-token variant (key of :data:`VARIANTS`).
        family: Pathway family (``GO_BP`` / ``GO_CC`` / ``GO_MF`` / ``REACTOME``
            / ``WIKIPATHWAYS``).
        v2_matrix: The parsed v2 matrix (defines the rows/cols to align to).
        v2_term2gene: Reconstructed TERM2GENE for the cell (pathway -> genes).
        target_sets: ``miRNA_display -> set(target_gene)`` (miRTarBase, symbol).
        universe_choice: Background policy (see module docstring).
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).

    Returns:
        A :class:`CellResult` whose ``values`` matrix has the SAME orientation as
        v2 (miRNA rows, pathway cols), every cell filled (``1.0`` where enricher
        scored no term, matching v2's not-tested sentinel).
    """
    # Restrict pathways to those present in BOTH the v2 matrix and the TERM2GENE.
    matrix_pathways = list(v2_matrix.pathways)
    t2g = _restrict_term2gene(v2_term2gene.term2gene, matrix_pathways)
    matched_pathways = [p for p in matrix_pathways if p in t2g]

    # Restrict miRNAs to the v2 matrix rows for which we have a target set.
    matrix_mirnas = list(v2_matrix.mirnas)
    cell_target_sets = {m: target_sets[m] for m in matrix_mirnas if m in target_sets}

    universe = _universe_for(universe_choice, set(v2_term2gene.universe), cell_target_sets)
    if universe_choice == "term2gene_union":
        # Defensive: ensure every TERM2GENE gene is in the universe (so the
        # effective N == |TERM2GENE genes|, the v2-faithful background).
        for g in t2g.values():
            universe |= g

    term2name = {p: p for p in matched_pathways}  # v2 term names are display names

    t0 = time.time()
    raw, prov = clusterprofiler.run_ora_batch_with_provenance(
        cell_target_sets,
        t2g,
        universe,
        term2name=term2name,
        rscript=rscript,
    )
    elapsed = time.time() - t0

    # Assemble in v2 orientation over the FULL v2 row/col label space so the
    # matrix lines up cell-for-cell with v2. Unscored / unmatched cells = 1.0.
    n_m, n_p = len(matrix_mirnas), len(matrix_pathways)
    values = np.full((n_m, n_p), NOT_TESTED, dtype=np.float64)
    col_idx = {p: j for j, p in enumerate(matrix_pathways)}
    for i, mirna in enumerate(matrix_mirnas):
        per_term = raw.get(mirna)
        if not per_term:
            continue
        for term, (pval, _padj, _ov) in per_term.items():
            j = col_idx.get(term)
            if j is not None:
                values[i, j] = pval

    return CellResult(
        organism=organism,
        variant=variant,
        family=family,
        universe_choice=universe_choice,
        mirnas=matrix_mirnas,
        pathways=matrix_pathways,
        values=values,
        n_universe=len(universe),
        n_terms_t2g=len(t2g),
        n_terms_matched=len(matched_pathways),
        elapsed_s=elapsed,
        provenance=prov,
    )


def write_cell_matrix(result: CellResult, out_path: str | Path) -> Path:
    """Write a cell's matrix as CSV in the v2 orientation (miRNA rows, path cols).

    The header row is the quoted pathway names with NO leading label cell (to
    match the v2 matrix file shape exactly); each data row is the miRNA name
    followed by its p-values.

    Returns the path written.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(result.pathways)
        for i, mirna in enumerate(result.mirnas):
            writer.writerow([mirna, *(repr(float(v)) for v in result.values[i])])
    return out
