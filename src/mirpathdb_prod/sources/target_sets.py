"""miRNA -> target-gene set builders (v2-comparable validated variants).

This module turns raw source rows (loaded by the per-source adapters in this
package) into the ``miRNA_display_name -> set(target_gene)`` maps the controlled
v2 reproduction consumes. It implements the two **validated_miRTarBase** variants
the validation requires. (The production HIGH/MEDIUM/LOW evidence rollup is a
separate, implemented module, :mod:`mirpathdb_prod.sources.production_target_sets`.)

Gene-ID space (the #1 silent-failure risk — see RECON):
    The v2 GeneTrail2 ORA ran in **gene-SYMBOL** space (mmu ``Info`` symbols such
    as ``H2-Q1``/``Smad7``; hsa symbols). miRTarBase ships both a ``Target Gene``
    SYMBOL and a ``Target Gene (Entrez ID)``. The default and v2-matching choice
    is ``id_space="symbol"`` (use the symbol column directly — no ``org.*.eg.db``
    mapping is needed because clusterProfiler ``enricher`` is ID-agnostic and the
    TERM2GENE drives the space). ``id_space="entrez"`` is kept as a sensitivity
    axis / fallback.

The two validated variants (controlled inputs — no predictor-version ambiguity):
    * ``validated_miRTarBase_strong`` — only rows whose ``Support Type`` equals
      the single STRONG functional value
      (:data:`mirpathdb_prod.sources.mirtarbase.STRONG_SUPPORT_TYPE`, i.e.
      reporter assay / Western blot / qPCR-grade functional evidence).
    * ``validated_miRTarBase_all`` — every miRTarBase MTI row (strong + weak),
      regardless of support type.

Pure transform: no I/O, no normalization of the miRNA label (it is already in
v2 display-name form), no judgment.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Literal

from mirpathdb_prod.sources import mirtarbase

# The two v2-comparable validated variants this builder produces.
MirtarbaseVariant = Literal["validated_miRTarBase_strong", "validated_miRTarBase_all"]
VALIDATED_VARIANTS: tuple[MirtarbaseVariant, ...] = (
    "validated_miRTarBase_strong",
    "validated_miRTarBase_all",
)

# Gene-ID spaces a target set can be emitted in.
IdSpace = Literal["symbol", "entrez"]

# miRTarBase column the miRNA display name comes from (already v2-row-label form).
_MIRNA_COL = "miRNA"
_SUPPORT_COL = "Support Type"
_SYMBOL_COL = "Target Gene"
_ENTREZ_COL = "Target Gene (Entrez ID)"


def _gene_value(row: dict[str, Any], id_space: IdSpace) -> str | None:
    """Extract the target-gene identifier in the requested space (or None)."""
    if id_space == "symbol":
        raw = row.get(_SYMBOL_COL)
    elif id_space == "entrez":
        raw = row.get(_ENTREZ_COL)
    else:  # pragma: no cover - guarded by Literal + explicit check below
        raise ValueError(f"unknown id_space {id_space!r}; use 'symbol' or 'entrez'")
    if raw is None:
        return None
    text = str(raw).strip()
    return text or None


def build_target_sets(
    mti_rows: Iterable[dict[str, Any]],
    *,
    variant: MirtarbaseVariant,
    id_space: IdSpace = "symbol",
) -> dict[str, set[str]]:
    """Build ``miRNA_display -> set(target_gene)`` from raw miRTarBase MTI rows.

    Args:
        mti_rows: Raw MTI row dicts as returned by
            :func:`mirpathdb_prod.sources.mirtarbase.load_mti` (keyed by the
            official column names).
        variant: Which validated variant to build. ``..._strong`` keeps only
            rows whose ``Support Type`` equals
            :data:`~mirpathdb_prod.sources.mirtarbase.STRONG_SUPPORT_TYPE`;
            ``..._all`` keeps every row.
        id_space: ``"symbol"`` (default, v2-matching) uses the ``Target Gene``
            symbol column; ``"entrez"`` uses ``Target Gene (Entrez ID)``.

    Returns:
        A dict mapping each miRNA display name to the (de-duplicated) set of its
        target genes in the requested ID space. miRNAs whose only rows lack a
        gene value in that space do not appear.

    Raises:
        ValueError: if ``variant`` or ``id_space`` is not recognised.
    """
    if variant not in VALIDATED_VARIANTS:
        raise ValueError(
            f"unknown variant {variant!r}; expected one of {VALIDATED_VARIANTS}"
        )
    if id_space not in ("symbol", "entrez"):
        raise ValueError(f"unknown id_space {id_space!r}; use 'symbol' or 'entrez'")

    strong_only = variant == "validated_miRTarBase_strong"
    result: dict[str, set[str]] = {}

    for row in mti_rows:
        if strong_only and row.get(_SUPPORT_COL) != mirtarbase.STRONG_SUPPORT_TYPE:
            continue
        mirna = row.get(_MIRNA_COL)
        if mirna is None:
            continue
        mirna_name = str(mirna).strip()
        if not mirna_name:
            continue
        gene = _gene_value(row, id_space)
        if gene is None:
            continue
        result.setdefault(mirna_name, set()).add(gene)

    return result


def summarize_target_sets(target_sets: dict[str, set[str]]) -> dict[str, int]:
    """Return small count summary of a built target-set map (for reports/tests).

    Keys: ``n_mirnas`` (miRNAs with ≥1 target), ``n_unique_genes`` (distinct
    target genes across all miRNAs), ``n_pairs`` (total miRNA-gene pairs),
    ``max_targets`` / ``min_targets`` (largest/smallest per-miRNA target count).
    """
    if not target_sets:
        return {
            "n_mirnas": 0,
            "n_unique_genes": 0,
            "n_pairs": 0,
            "max_targets": 0,
            "min_targets": 0,
        }
    sizes = [len(v) for v in target_sets.values()]
    genes: set[str] = set()
    for v in target_sets.values():
        genes |= v
    return {
        "n_mirnas": len(target_sets),
        "n_unique_genes": len(genes),
        "n_pairs": sum(sizes),
        "max_targets": max(sizes),
        "min_targets": min(sizes),
    }
