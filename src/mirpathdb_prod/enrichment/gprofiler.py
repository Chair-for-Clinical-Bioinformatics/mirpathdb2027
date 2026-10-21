"""g:Profiler over-representation engine (secondary; current GO).

g:Profiler is the production *secondary* functional-enrichment engine,
complementing clusterProfiler's GO ORA. It is used here for the
production-preview comparison (current GO), not the controlled v2 reproduction
(that is clusterProfiler with v2's reconstructed gene sets).

Reproducibility (REQUIRED, see ``docs/enrichment.md``): g:Profiler advances its
underlying datasource DB on a rolling schedule, so a reproducible build MUST pin
the **archived datasource DB version** rather than querying "current". g:Profiler
exposes a stable archived-version selector via ``base_url``
(e.g. ``https://biit.cs.ut.ee/gprofiler_archive3/e<ENSEMBL>_eg<EG>_pXX``). The
engine refuses to run "current" unless the caller explicitly opts in.

The ``gprofiler-official`` client is imported lazily inside :meth:`run` so this
module imports with no network and no hard dependency at import time.
"""

from __future__ import annotations

from mirpathdb_prod.enrichment.base import (
    ENGINE_GPROFILER,
    EnrichmentRequest,
    EnrichmentResult,
)

# Canonical taxid -> g:Profiler organism code.
GPROFILER_ORGANISM_CODES: dict[int, str] = {
    9606: "hsapiens",
    10090: "mmusculus",
    10116: "rnorvegicus",
    9544: "mmulatta",  # Macaca mulatta (macaque — predicted-dominant; lone validated MTI marked)
}

# Provenance pin slots a real run must populate.
PROVENANCE_PINS: tuple[str, ...] = (
    "gprofiler_version",
    "datasource_db_version",  # the archived g:Profiler DB snapshot — pin, never "current"
    "archive_base_url",
)

# GO sources g:Profiler exposes (the production-preview comparison restricts to
# GO to line up with clusterProfiler's GO families).
GO_SOURCES: tuple[str, ...] = ("GO:BP", "GO:CC", "GO:MF")


class GProfilerEngine:
    """g:Profiler engine adapter (implements :class:`EnrichmentEngine`).

    ``run`` requires either an explicit archived ``base_url`` in
    ``request.params['base_url']`` (the reproducible path — pins the datasource
    DB snapshot) or an explicit ``allow_current=True`` opt-in (non-reproducible,
    for ad-hoc previews only).
    """

    name = ENGINE_GPROFILER

    def run(self, request: EnrichmentRequest) -> EnrichmentResult:
        """Run functional enrichment via g:Profiler (secondary engine).

        Args:
            request: The engine-agnostic request. ``request.organism`` selects
                the g:Profiler organism code (see :data:`GPROFILER_ORGANISM_CODES`).
                ``request.params`` knobs:

                - ``base_url`` (str): archived g:Profiler endpoint pinning the
                    datasource DB version (REQUIRED unless ``allow_current``).
                - ``allow_current`` (bool): opt in to the rolling "current" DB
                    (non-reproducible). Default ``False``.
                - ``sources`` (list[str]): g:Profiler sources to query
                    (default :data:`GO_SOURCES`).
                - ``no_iea`` (bool): exclude electronic GO annotations.

        Returns:
            An :class:`EnrichmentResult` with one row per enriched term
            (``term`` / ``name`` / ``source`` / ``pvalue`` / ``overlap``) and
            provenance pinning the gprofiler client + archived DB base URL.

        Raises:
            ValueError: if the organism is not a g:Profiler-covered production
                organism, or if neither a ``base_url`` nor ``allow_current`` is
                given (refuses an unpinned, non-reproducible run by default).
        """
        organism_code = GPROFILER_ORGANISM_CODES.get(request.organism)
        if organism_code is None:
            raise ValueError(
                f"organism taxid {request.organism} not covered by g:Profiler "
                f"production set; known: {sorted(GPROFILER_ORGANISM_CODES)}"
            )

        base_url = request.params.get("base_url")
        allow_current = bool(request.params.get("allow_current", False))
        if not base_url and not allow_current:
            raise ValueError(
                "g:Profiler requires a pinned archived datasource DB: pass "
                "request.params['base_url'] (an archived gprofiler_archive3 "
                "endpoint) for a reproducible run, or set "
                "request.params['allow_current']=True to explicitly accept the "
                "rolling 'current' DB (non-reproducible)"
            )

        sources = request.params.get("sources") or list(GO_SOURCES)
        no_iea = bool(request.params.get("no_iea", False))

        # Lazy import so the module loads (and tests collect) without the client
        # or any network present.
        import gprofiler as _gprofiler_pkg
        from gprofiler import GProfiler

        client_kwargs: dict[str, object] = {"return_dataframe": True}
        if base_url:
            client_kwargs["base_url"] = str(base_url)
        gp = GProfiler(**client_kwargs)  # type: ignore[arg-type]

        background = list(request.background) if request.background is not None else None
        profile = gp.profile(
            organism=organism_code,
            query=list(request.genes),
            sources=list(sources),
            background=background,
            no_evidences=False,
            no_iea=no_iea,
            user_threshold=1.0,  # return all terms; we keep the raw p
            significance_threshold_method="fdr",
        )

        rows: list[dict[str, object]] = []
        for record in profile.to_dict(orient="records"):
            rows.append(
                {
                    "term": record.get("native"),
                    "name": record.get("name"),
                    "source": record.get("source"),
                    "pvalue": float(record.get("p_value", 1.0)),
                    "overlap": int(record.get("intersection_size", 0) or 0),
                }
            )

        provenance: dict[str, object] = dict.fromkeys(PROVENANCE_PINS)
        provenance["gprofiler_version"] = getattr(_gprofiler_pkg, "__version__", None)
        provenance["archive_base_url"] = str(base_url) if base_url else "current"
        provenance["datasource_db_version"] = (
            str(base_url) if base_url else "current (UNPINNED — non-reproducible)"
        )
        provenance["organism_code"] = organism_code

        return EnrichmentResult(
            engine=self.name,
            request=request,
            rows=rows,
            provenance=provenance,
        )

    @staticmethod
    def organism_code(taxid: int) -> str:
        """Return the g:Profiler organism code for a production taxid."""
        code = GPROFILER_ORGANISM_CODES.get(taxid)
        if code is None:
            raise ValueError(
                f"organism taxid {taxid} not covered by g:Profiler production set"
            )
        return code
