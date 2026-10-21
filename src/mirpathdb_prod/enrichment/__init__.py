"""Enrichment engine (over-representation analysis).

PRODUCTION DECISION (refined plan): the prototype's custom ORA/BH kernel is
RETIRED. Production enrichment is delegated to established external engines:

* **clusterProfiler** (Bioconductor / R) — over-representation (GO + Reactome +
  WikiPathways) against the ``org.*.eg.db`` / GO.db / reactome.db annotations.
  This is the **sole production engine** for the current build: it generated all
  production matrices (clusterProfiler 4.18.4 / R 4.5.3) over hsa/mmu/rno/mml —
  see ``provenance/production_pathways_manifest.json``.
* **g:Profiler** — multi-source functional enrichment against an archived
  g:Profiler datasource database snapshot. **Implemented but secondary**: a
  working client, **not run in the current production build** (preview).

Both are exposed behind ONE common interface,
:func:`mirpathdb_prod.enrichment.base.run_enrichment`, so callers never depend
on which engine ran. The dispatch layer here holds no R or network logic of its
own; the engines do. See ``docs/enrichment.md`` for the full design and the
version-pinning provenance requirement.

Consequence to remember: with a different engine and no KEGG, miRPathDB2
*numeric* reproduction becomes **DIRECTIONAL** (concordance of trends, not exact
p-values); organism concordance stays an **EXACT** gate.
"""

from __future__ import annotations

from mirpathdb_prod.enrichment.base import (
    EnrichmentEngine,
    EnrichmentRequest,
    EnrichmentResult,
    get_engine,
    run_enrichment,
)

__all__ = [
    "EnrichmentEngine",
    "EnrichmentRequest",
    "EnrichmentResult",
    "get_engine",
    "run_enrichment",
]
