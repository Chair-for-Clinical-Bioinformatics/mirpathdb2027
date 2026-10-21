"""Pathway-source layer.

GO is the PRIMARY pathway source for the published DB. KEGG is EXCLUDED from the
published bundle (redistribution license) but kept as an OPTIONAL, license-gated,
USER-RUNNABLE source in the released open pipeline. Reactome and WikiPathways are
under evaluation (flagged for Tuesday). See ``docs/pathway_sources.md``.

This package ships the registry that reads ``config/pathway_sources.yaml`` and
exposes the role / provenance / version / license + redistribution flag of each
source. No pathway data is ingested here.

Per-organism pathway coverage is native (``org.*.eg.db`` GO + ``reactome.db``
Reactome + WikiPathways GMT) for hsa/mmu/rno. **Macaque (mml) is the exception:**
Reactome ships no curated ``R-MML-`` pathways, so macaque Reactome is **INFERRED**
(orthology-projected from human ``R-HSA-``) by
:mod:`mirpathdb_prod.pathways.inferred_reactome_genesets`; macaque WikiPathways is
NONE.
"""

from __future__ import annotations

from mirpathdb_prod.pathways.inferred_reactome_genesets import (
    InferredReactomeGeneSets,
    build_inferred_macaque_reactome,
)
from mirpathdb_prod.pathways.registry import (
    PathwaySource,
    PathwaySourceRegistry,
    load_pathway_sources,
)

__all__ = [
    "InferredReactomeGeneSets",
    "PathwaySource",
    "PathwaySourceRegistry",
    "build_inferred_macaque_reactome",
    "load_pathway_sources",
]
