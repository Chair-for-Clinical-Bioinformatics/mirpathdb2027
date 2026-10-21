"""Independent hypergeometric ORA kernel + the two BH correction families.

Model (identical to the production ``clusterProfiler::enricher`` call with
``pvalueCutoff = qvalueCutoff = 1``, ``minGSSize = 1``, ``maxGSSize = 1e8``,
``pAdjustMethod = "BH"``; see ``DOSE:::enricher_internal``):

* ``N`` = background size (genes of the family universe),
* ``K`` = membership size of a term within the background,
* ``n`` = query genes inside the background,
* ``k`` = query genes inside the term,
* raw ``P = Pr[X >= k]`` for ``X ~ Hypergeom(N, K, n)`` (R: ``phyper(k-1, K, N-K,
  n, lower.tail = FALSE)``).

Correction families:

* **overlap-only BH** (the production procedure): BH over exactly the terms with
  ``k >= 1`` (``m`` = their number); terms with ``k = 0`` are never scored.
* **fixed-family BH** (comparator): BH over every eligible term of the cell
  (``m_full`` = eligible terms, defined from the background alone), assigning the
  exact ``P = Pr[X >= 0] = 1`` to zero-overlap terms.

This module is an *independent* implementation (numpy + scipy). Its numerical
agreement with the pinned R kernel is measured by :mod:`.r_kernel`; until that
agreement is on record for a case class, results from here are labelled
"independent implementation", not "production-kernel run".
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.stats import hypergeom

from mirpathdb_prod.validation.robustness.inputs import FamilySpace


def hypergeom_sf(k: np.ndarray, N: int, K: np.ndarray, n: np.ndarray | int) -> np.ndarray:
    """``Pr[X >= k]`` for ``X ~ Hypergeom(N, K, n)`` (vectorised, exact tail)."""
    k = np.asarray(k, dtype=np.int64)
    K = np.asarray(K, dtype=np.int64)
    n = np.asarray(n, dtype=np.int64)
    out = hypergeom.sf(k - 1, N, K, n)
    # sf() of k-1 already returns exactly 1.0 for k == 0; guard the rounding edge.
    return np.clip(np.asarray(out, dtype=np.float64), 0.0, 1.0)


def hypergeom_sf_unique(k: np.ndarray, K: np.ndarray, n: int, N: int) -> np.ndarray:
    """Same as :func:`hypergeom_sf` for a fixed ``n``, evaluated once per unique (K, k)."""
    k = np.asarray(k, dtype=np.int64)
    K = np.asarray(K, dtype=np.int64)
    if k.size == 0:
        return np.zeros(0, dtype=np.float64)
    pairs = np.stack([K, k], axis=1)
    uniq, inverse = np.unique(pairs, axis=0, return_inverse=True)
    p_uniq = hypergeom_sf(uniq[:, 1], N, uniq[:, 0], n)
    return p_uniq[np.asarray(inverse).ravel()]


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values over ``p`` (``m = len(p)``), input order."""
    p = np.asarray(p, dtype=np.float64)
    m = p.size
    if m == 0:
        return np.zeros(0, dtype=np.float64)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    adj_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m, dtype=np.float64)
    adj[order] = np.minimum(adj_sorted, 1.0)
    return adj


def bh_fixed_family(p: np.ndarray, m_full: int) -> np.ndarray:
    """BH over a family of ``m_full`` tests where the ``m_full - len(p)`` unlisted tests have P = 1.

    Equivalent to ``bh_adjust(concat(p, ones(m_full - len(p))))`` restricted to
    the listed entries: the trailing P = 1 tests adjust to exactly 1 and never
    lower any listed adjusted value.
    """
    p = np.asarray(p, dtype=np.float64)
    m = p.size
    if m_full < m:
        raise ValueError(f"m_full={m_full} < number of listed p-values {m}")
    if m == 0:
        return np.zeros(0, dtype=np.float64)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m_full / np.arange(1, m + 1)
    adj_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m, dtype=np.float64)
    adj[order] = np.minimum(adj_sorted, 1.0)
    return adj


def fold_enrichment(k: np.ndarray, n: int, K: np.ndarray, N: int) -> np.ndarray:
    """``(k/n) / (K/N)`` (clusterProfiler's ``FoldEnrichment``)."""
    k = np.asarray(k, dtype=np.float64)
    K = np.asarray(K, dtype=np.float64)
    if n == 0:
        return np.full(k.shape, np.nan)
    return (k / n) / (K / N)


@dataclass
class OraBatch:
    """ORA results for a batch of queries against ONE :class:`FamilySpace`.

    Per-query results are ragged (only terms with ``k >= 1`` are listed, exactly
    as the production kernel reports them); ``m_full`` is the fixed-family size.
    """

    N: int
    m_full: int
    n: np.ndarray  # (B,) query size inside the background
    hits: list[np.ndarray]  # per query: indices (into space.terms) with k >= 1
    k: list[np.ndarray]
    K: list[np.ndarray]
    p: list[np.ndarray]
    padj_overlap: list[np.ndarray]
    padj_fixed: list[np.ndarray]

    @property
    def size(self) -> int:
        return len(self.hits)

    def m_overlap(self) -> np.ndarray:
        return np.array([h.size for h in self.hits], dtype=np.int64)

    def n_rejections(self, alpha: float, family: str = "overlap") -> np.ndarray:
        arrs = self.padj_overlap if family == "overlap" else self.padj_fixed
        return np.array([int((a <= alpha).sum()) for a in arrs], dtype=np.int64)

    def rejected_terms(self, i: int, alpha: float, family: str = "overlap") -> np.ndarray:
        arr = (self.padj_overlap if family == "overlap" else self.padj_fixed)[i]
        return self.hits[i][arr <= alpha]

    def as_dict(self, i: int, terms: Sequence[str]) -> dict[str, tuple[float, float, int]]:
        """``{term: (p, padj_overlap, k)}`` in the production result shape."""
        return {
            terms[t]: (float(p), float(q), int(kk))
            for t, p, q, kk in zip(
                self.hits[i], self.p[i], self.padj_overlap[i], self.k[i], strict=True
            )
        }


def _query_matrix(N: int, queries: Sequence[np.ndarray]) -> sp.csc_matrix:
    """``genes x B`` 0/1 CSC matrix of the query indicator vectors."""
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    for b, q in enumerate(queries):
        q = np.asarray(q, dtype=np.int64)
        if q.size and (q.min() < 0 or q.max() >= N):
            raise ValueError(f"query {b}: gene index outside the background (N={N})")
        rows.append(q)
        cols.append(np.full(q.size, b, dtype=np.int64))
    if rows:
        r = np.concatenate(rows)
        c = np.concatenate(cols)
    else:
        r = np.zeros(0, dtype=np.int64)
        c = np.zeros(0, dtype=np.int64)
    # int32 counts: an int8 product would wrap at k = 128 (overlaps of several
    # hundred genes are routine for large queries against broad GO terms).
    mat = sp.coo_matrix((np.ones(r.size, dtype=np.int32), (r, c)), shape=(N, len(queries))).tocsc()
    mat.data[:] = 1
    return mat


def ora_batch(space: FamilySpace, queries: Sequence[np.ndarray], *, chunk: int = 512) -> OraBatch:
    """Score every query of ``queries`` (gene-index arrays) against ``space``.

    Queries are expected to lie inside the background (indices into
    ``space.genes``); a query is the unique gene set, duplicates are ignored.
    An empty query yields no scored term (``m = 0``), matching the production
    runner which never submits empty sets.
    """
    N = space.N
    K_all = space.term_sizes
    m_full = space.n_terms
    n = np.array([np.unique(np.asarray(q, dtype=np.int64)).size for q in queries], dtype=np.int64)
    hits: list[np.ndarray] = []
    ks: list[np.ndarray] = []
    Ks: list[np.ndarray] = []
    ps: list[np.ndarray] = []
    padj_o: list[np.ndarray] = []
    padj_f: list[np.ndarray] = []
    A = space.matrix.astype(np.int32).tocsr()  # never let counts overflow a narrow dtype
    for start in range(0, len(queries), chunk):
        block = [np.unique(np.asarray(q, dtype=np.int64)) for q in queries[start : start + chunk]]
        if not block:
            continue
        Q = _query_matrix(N, block)
        C = (A @ Q).tocsc()  # terms x B counts
        for j, q in enumerate(block):
            col = C.getcol(j)
            hit = col.indices.astype(np.int64)
            kk = col.data.astype(np.int64)
            order = np.argsort(hit, kind="mergesort")
            hit, kk = hit[order], kk[order]
            KK = K_all[hit]
            n_q = int(q.size)
            p = hypergeom_sf_unique(kk, KK, n_q, N) if hit.size else np.zeros(0)
            hits.append(hit)
            ks.append(kk)
            Ks.append(KK)
            ps.append(p)
            padj_o.append(bh_adjust(p))
            padj_f.append(bh_fixed_family(p, m_full))
    return OraBatch(
        N=N, m_full=m_full, n=n, hits=hits, k=ks, K=Ks, p=ps, padj_overlap=padj_o, padj_fixed=padj_f
    )


def ora_single(space: FamilySpace, query: np.ndarray) -> OraBatch:
    """Convenience wrapper: one query."""
    return ora_batch(space, [np.asarray(query, dtype=np.int64)])
