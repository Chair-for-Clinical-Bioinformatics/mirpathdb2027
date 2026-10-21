"""Full-family Benjamini-Hochberg adjustment for the production ORA.

The production kernel (``clusterProfiler::enricher``) returns, per query, only the
pathways that share at least one gene with the query (raw ``P = Pr[X >= k]``,
``k >= 1``) and adjusts over exactly those ``m`` tests. This module defines the
alternative, query-independent family:

* ``M`` = number of *eligible* pathway terms of a cell = terms of the complete
  annotation input whose membership inside the background is non-empty and within
  the production size limits (``minGSSize = 1``, ``maxGSSize = 1e8``). ``M`` is
  fixed by the inputs before any query is seen.
* every eligible term is one test; a term with zero overlap has the exact raw
  ``P = Pr[X >= 0] = 1``;
* Benjamini-Hochberg is applied across the ``M`` tests.

Because the ``M - m`` unlisted tests all have ``P = 1`` they adjust to exactly 1
and never lower a listed adjusted value, so the full-family adjusted value of a
listed term is ``min(1, min_{j >= rank} (M / j) * p_(j))`` over the ``m`` listed
tests -- i.e. the overlap-only BH computation with numerator ``M`` instead of
``m``. Equivalently, an already computed overlap-only adjusted value migrates
exactly as ``q_full = min(1, q_overlap * M / m)`` **provided** ``m`` is the
complete overlap-tested count (every returned term, including valid ``P = 1``
tests). ``m`` must never be taken from a retained/filtered row set or inferred
from ``P < 1``; :func:`verify_migration` checks the identity against the direct
computation and is the guard the export path uses.

Raw P-values and overlaps are never modified here.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import numpy as np

# Production enricher size limits (see enrichment/clusterprofiler.py, _MAX_GS_SIZE).
PRODUCTION_MIN_GS_SIZE = 1
PRODUCTION_MAX_GS_SIZE = 100_000_000


def eligible_family(
    term2gene: Mapping[str, Iterable[str]],
    universe: Iterable[str],
    *,
    min_gs_size: int = PRODUCTION_MIN_GS_SIZE,
    max_gs_size: int = PRODUCTION_MAX_GS_SIZE,
) -> tuple[int, list[str]]:
    """Return ``(M, eligible_terms)`` for one cell from its complete inputs.

    A term is eligible when its membership restricted to ``universe`` has a size
    inside ``[min_gs_size, max_gs_size]``. The universe is the cell's background
    exactly as handed to the enricher; genes outside it never count.
    """
    uni = set(str(g).strip() for g in universe if str(g).strip())
    eligible: list[str] = []
    for term in sorted(term2gene):
        size = len({str(g).strip() for g in term2gene[term] if str(g).strip()} & uni)
        if min_gs_size <= size <= max_gs_size:
            eligible.append(term)
    return len(eligible), eligible


def full_family_adjust(pvalues: Sequence[float], m_full: int) -> np.ndarray:
    """BH-adjusted values over a family of ``m_full`` tests, ``len(pvalues)`` of them listed.

    The ``m_full - len(pvalues)`` unlisted tests have ``P = 1``. Result is in the
    input order. ``m_full`` must be >= the number of listed tests.
    """
    p = np.asarray(pvalues, dtype=np.float64)
    m = p.size
    if m_full < m:
        raise ValueError(f"family size M={m_full} smaller than the {m} listed tests")
    if m == 0:
        return np.zeros(0, dtype=np.float64)
    if np.any(p < 0) or np.any(p > 1) or np.any(np.isnan(p)):
        raise ValueError("raw p-values must lie in [0, 1]")
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m_full / np.arange(1, m + 1)
    adj_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m, dtype=np.float64)
    out[order] = np.minimum(adj_sorted, 1.0)
    return out


def migrate_overlap_adjusted(q_overlap: Sequence[float], m_overlap: int, m_full: int) -> np.ndarray:
    """Exact migration ``q_full = min(1, q_overlap * M / m)``.

    ``m_overlap`` MUST be the complete number of overlap-tested terms of that
    query (all terms the enricher returned, including valid ``P = 1`` tests),
    never a retained-row count. Values already clipped at 1 stay 1.
    """
    if m_overlap <= 0:
        raise ValueError("m_overlap must be the positive complete overlap-tested count")
    if m_full < m_overlap:
        raise ValueError(f"family size M={m_full} smaller than m={m_overlap}")
    q = np.asarray(q_overlap, dtype=np.float64)
    return np.minimum(1.0, q * (m_full / m_overlap))


def verify_migration(
    pvalues: Sequence[float],
    q_overlap: Sequence[float],
    m_full: int,
    *,
    rtol: float = 1e-12,
) -> dict[str, object]:
    """Check ``migrate_overlap_adjusted`` against the direct full-family computation.

    ``pvalues``/``q_overlap`` must be the COMPLETE listed family of one query
    (``len == m``). Returns a report with the maximum relative difference and a
    ``passed`` flag; a failure means ``m`` is not the complete tested count (or
    ``q_overlap`` was not produced by BH over exactly these ``m`` tests).
    """
    p = np.asarray(pvalues, dtype=np.float64)
    q = np.asarray(q_overlap, dtype=np.float64)
    if p.size != q.size:
        raise ValueError("pvalues and q_overlap must have the same length (the complete family)")
    m = p.size
    direct = full_family_adjust(p, m_full)
    migrated = migrate_overlap_adjusted(q, m, m_full) if m else np.zeros(0)
    scale = np.maximum(np.maximum(np.abs(direct), np.abs(migrated)), 1e-300)
    rel = np.abs(direct - migrated) / scale
    # Also confirm q_overlap itself is BH over exactly these m tests.
    q_recomputed = full_family_adjust(p, m) if m else np.zeros(0)
    rel_q = np.abs(q_recomputed - q) / np.maximum(
        np.maximum(np.abs(q), np.abs(q_recomputed)), 1e-300
    )
    return {
        "m_overlap": int(m),
        "m_full": int(m_full),
        "max_rel_diff_direct_vs_migrated": float(rel.max()) if m else 0.0,
        "max_rel_diff_q_overlap_recomputed": float(rel_q.max()) if m else 0.0,
        "passed": bool((rel.max() <= rtol if m else True) and (rel_q.max() <= rtol if m else True)),
    }


def adjust_cell_results(
    results: Mapping[str, Mapping[str, tuple[float, float, int]]],
    m_full: int,
) -> dict[str, dict[str, tuple[float, float, float, int]]]:
    """Add the full-family adjusted value to unfiltered per-query enricher results.

    ``results`` is ``{set_name: {term: (p, q_overlap, overlap)}}`` exactly as
    :func:`mirpathdb_prod.enrichment.clusterprofiler.run_ora_batch` returns it
    (every scored term, no filtering). Returns ``{set_name: {term: (p, q_overlap,
    q_full, overlap)}}``; raw P and overlaps are passed through unchanged. The
    migration identity is verified for every set and a violation raises, because
    it would mean the listed family is not the complete tested family.
    """
    out: dict[str, dict[str, tuple[float, float, float, int]]] = {}
    for name, rows in results.items():
        terms = list(rows)
        p = np.array([rows[t][0] for t in terms], dtype=np.float64)
        q = np.array([rows[t][1] for t in terms], dtype=np.float64)
        if terms:
            check = verify_migration(p, q, m_full)
            if not check["passed"]:
                raise ValueError(
                    f"{name}: overlap-only adjusted values are not BH over the {len(terms)} "
                    f"listed tests (max rel diff "
                    f"{check['max_rel_diff_q_overlap_recomputed']:.2e}); "
                    f"refusing to migrate a filtered family"
                )
        q_full = full_family_adjust(p, m_full) if terms else np.zeros(0)
        out[name] = {
            t: (float(rows[t][0]), float(rows[t][1]), float(q_full[i]), int(rows[t][2]))
            for i, t in enumerate(terms)
        }
    return out
