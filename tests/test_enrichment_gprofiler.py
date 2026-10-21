"""Tests for the g:Profiler engine seam (enrichment/gprofiler.py).

g:Profiler is the secondary production-preview engine; the real call hits the
network, so these tests cover the offline contracts only: organism-code mapping,
the locked GO source set, and the reproducibility guard that REFUSES to run
against an unpinned 'current' datasource DB unless explicitly opted in.
"""

from __future__ import annotations

import pytest

from mirpathdb_prod.enrichment import gprofiler as gp
from mirpathdb_prod.enrichment.base import EnrichmentRequest


def test_organism_codes_cover_production_set() -> None:
    assert set(gp.GPROFILER_ORGANISM_CODES) == {9606, 10090, 10116, 9544}
    assert gp.GPROFILER_ORGANISM_CODES[9606] == "hsapiens"
    assert gp.GPROFILER_ORGANISM_CODES[10090] == "mmusculus"


def test_organism_code_helper() -> None:
    assert gp.GProfilerEngine.organism_code(9606) == "hsapiens"
    with pytest.raises(ValueError, match="not covered"):
        gp.GProfilerEngine.organism_code(99999)


def test_go_sources_locked() -> None:
    assert gp.GO_SOURCES == ("GO:BP", "GO:CC", "GO:MF")


def test_run_rejects_unknown_organism() -> None:
    engine = gp.GProfilerEngine()
    req = EnrichmentRequest(genes=("A",), organism=99999, pathway_source="go")
    with pytest.raises(ValueError, match="not covered by g:Profiler"):
        engine.run(req)


def test_run_refuses_unpinned_current_db() -> None:
    # The reproducibility guard: no base_url and no allow_current -> refuse,
    # BEFORE any network/client import.
    engine = gp.GProfilerEngine()
    req = EnrichmentRequest(genes=("A",), organism=9606, pathway_source="go")
    with pytest.raises(ValueError, match="pinned archived datasource"):
        engine.run(req)
