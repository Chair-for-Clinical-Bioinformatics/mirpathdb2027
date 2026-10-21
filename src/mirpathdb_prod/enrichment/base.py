"""Common enrichment interface (engine-agnostic dispatch layer).

Defines the single contract every enrichment engine must satisfy and the
request/result containers shared across engines. This is the interface the rest
of the build codes against; the concrete engines it dispatches to are
implemented — :mod:`~.clusterprofiler` (the R-backed clusterProfiler engine that
produced the production matrices) and :mod:`~.gprofiler` (an implemented
g:Profiler client, secondary / not run in the current production build).

The whole point of this interface is that callers express *what* they want
enriched (a gene set, an organism, a pathway source) without knowing *which*
engine answers — so the engine can be swapped or pinned without touching call
sites. No engine logic, R, or network lives in this module itself; it only
assembles the request and dispatches.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

# Engine identifiers (mirror the keys accepted by get_engine / run_enrichment).
ENGINE_CLUSTERPROFILER = "clusterprofiler"
ENGINE_GPROFILER = "gprofiler"
KNOWN_ENGINES: tuple[str, ...] = (ENGINE_CLUSTERPROFILER, ENGINE_GPROFILER)


@dataclass(frozen=True)
class EnrichmentRequest:
    """An over-representation request, independent of the engine that runs it.

    Attributes:
        genes: The query gene set (e.g. a miRNA's predicted/validated targets),
            as source identifiers to be mapped by the engine.
        organism: Canonical NCBI taxid of the organism the gene set belongs to
            (selects the correct annotation, e.g. org.Hs/Mm/Rn/Cf/Gg.eg.db).
        pathway_source: Registry id of the pathway library to test against
            (e.g. "go" — primary; see :mod:`mirpathdb_prod.pathways.registry`).
        background: Optional explicit background/universe gene set; ``None``
            lets the engine use its default universe (recorded in provenance).
        params: Engine-agnostic knobs (e.g. p-value cutoff, min/max set size).
    """

    genes: Sequence[str]
    organism: int
    pathway_source: str
    background: Sequence[str] | None = None
    params: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class EnrichmentResult:
    """The outcome of an enrichment run, plus the provenance to reproduce it.

    The ``provenance`` mapping is REQUIRED to pin every version that affects the
    numbers (engine version; for clusterProfiler the ``org.*.eg.db`` annotation
    package versions; for g:Profiler the archived datasource DB version). See
    ``docs/enrichment.md``.

    Attributes:
        engine: Which engine produced this result (an entry of KNOWN_ENGINES).
        request: The request that was run (echoed for auditability).
        rows: One dict per scored term (term id, name, counts, p, p.adjust,
            mapped genes).
        provenance: Version pins required for reproducibility (engine +
            annotation/datasource DB versions + retrieval date).
    """

    engine: str
    request: EnrichmentRequest
    rows: list[dict[str, object]] = field(default_factory=list)
    provenance: dict[str, object] = field(default_factory=dict)


@runtime_checkable
class EnrichmentEngine(Protocol):
    """Protocol every enrichment engine implements.

    Concrete engines (clusterProfiler, g:Profiler) provide ``name`` and
    ``run``; callers reach them through :func:`run_enrichment` and need not know
    the concrete type.
    """

    #: Engine identifier (one of KNOWN_ENGINES).
    name: str

    def run(self, request: EnrichmentRequest) -> EnrichmentResult:
        """Run an over-representation analysis for ``request``."""
        ...


def get_engine(engine: str) -> EnrichmentEngine:
    """Return the implemented engine adapter for ``engine``.

    Args:
        engine: One of :data:`KNOWN_ENGINES`
            (``"clusterprofiler"`` or ``"gprofiler"``).

    Returns:
        An :class:`EnrichmentEngine` adapter instance.

    Raises:
        ValueError: if ``engine`` is not a known engine.
    """
    # Imported lazily so the (heavy, later) engine modules are not pulled in
    # just to define the interface.
    if engine == ENGINE_CLUSTERPROFILER:
        from mirpathdb_prod.enrichment.clusterprofiler import ClusterProfilerEngine

        return ClusterProfilerEngine()
    if engine == ENGINE_GPROFILER:
        from mirpathdb_prod.enrichment.gprofiler import GProfilerEngine

        return GProfilerEngine()
    raise ValueError(f"unknown enrichment engine {engine!r}; known: {KNOWN_ENGINES}")


def run_enrichment(
    genes: Sequence[str],
    organism: int,
    pathway_source: str,
    *,
    engine: str = ENGINE_CLUSTERPROFILER,
    background: Sequence[str] | None = None,
    **params: object,
) -> EnrichmentResult:
    """Run an over-representation analysis through the chosen engine.

    The one common entrypoint the rest of the build uses. It assembles an
    :class:`EnrichmentRequest` and dispatches to the engine selected by
    ``engine`` (default clusterProfiler, the R-backed engine that produced the
    production matrices; ``"gprofiler"`` selects the implemented but secondary
    g:Profiler client).

    Args:
        genes: Query gene set (source identifiers).
        organism: Canonical NCBI taxid (selects annotation / datasource).
        pathway_source: Registry id of the pathway library (e.g. "go").
        engine: Engine to use (default: clusterProfiler).
        background: Optional explicit universe gene set.
        **params: Engine-agnostic parameters (cutoffs, set-size bounds, ...);
            the clusterProfiler engine also reads ``term2gene`` / ``term2name``
            here.

    Returns:
        An :class:`EnrichmentResult` (one row per scored term + pinned
        provenance).

    Raises:
        ValueError: if ``engine`` is not a known engine.
        NotImplementedError: from the clusterProfiler engine **only** when no
            ``term2gene`` is supplied in ``params`` (materialising GO
            ``TERM2GENE`` from ``org.*.eg.db`` inside this entrypoint is
            deferred; production callers pass ``term2gene`` explicitly).
    """
    request = EnrichmentRequest(
        genes=tuple(genes),
        organism=organism,
        pathway_source=pathway_source,
        background=tuple(background) if background is not None else None,
        params=dict(params),
    )
    return get_engine(engine).run(request)
