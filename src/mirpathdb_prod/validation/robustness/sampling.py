"""Deterministic samplers: uniform-gene null and degree-binned (target-frequency-aware) null.

Seeds are derived from the study seed and the configuration coordinates via
``numpy.random.SeedSequence`` so every (stratum, query size, replicate block)
is reproducible and independent of execution order or resume state.
"""

from __future__ import annotations

import zlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from mirpathdb_prod.validation.robustness.inputs import FamilySpace


def _crc(part: object) -> int:
    return zlib.crc32(str(part).encode("utf-8")) & 0xFFFFFFFF


def config_rng(base_seed: int, *parts: object) -> np.random.Generator:
    """A generator seeded from ``base_seed`` + the CRC32 of each coordinate."""
    ss = np.random.SeedSequence([int(base_seed), *(_crc(p) for p in parts)])
    return np.random.Generator(np.random.PCG64(ss))


def uniform_samples(rng: np.random.Generator, N: int, n: int, reps: int) -> list[np.ndarray]:
    """``reps`` uniform samples of ``n`` distinct gene indices from ``range(N)`` (sorted)."""
    if n < 0 or n > N:
        raise ValueError(f"cannot draw n={n} genes without replacement from N={N}")
    return [np.sort(rng.choice(N, size=n, replace=False)) for _ in range(reps)]


def uniform_pool_samples(
    rng: np.random.Generator, pool: np.ndarray, n: int, reps: int
) -> list[np.ndarray]:
    """``reps`` uniform samples of ``n`` distinct genes drawn from ``pool`` only (sorted).

    Control for the degree-matched null: same restriction to the tier's target-union
    pool (genes with degree >= 1) but no degree matching, so pool restriction and
    degree matching can be told apart.
    """
    pool = np.unique(np.asarray(pool, dtype=np.int64))
    if n < 0 or n > pool.size:
        raise ValueError(f"cannot draw n={n} genes without replacement from a pool of {pool.size}")
    return [np.sort(rng.choice(pool, size=n, replace=False)) for _ in range(reps)]


def network_degrees(space: FamilySpace, queries: Mapping[str, np.ndarray]) -> np.ndarray:
    """Gene degree = number of eligible miRNAs (of one tier) targeting the gene.

    Computed on the complete organism/tier target network inside the family
    background (all eligible miRNAs, including any later focal miRNA -- the
    focal miRNA's own edges are NOT removed).
    """
    deg = np.zeros(space.N, dtype=np.int64)
    for idx in queries.values():
        deg[np.unique(np.asarray(idx, dtype=np.int64))] += 1
    return deg


@dataclass
class DegreeBins:
    """Quantile bins of gene degree, defined without any pathway label.

    Genes with degree 0 (never targeted in the tier) get bin ``-1`` and are never
    sampled: a degree-matched draw can only contain targeted genes, which makes
    this null a target-network-restricted null by construction.
    """

    degrees: np.ndarray  # (N,)
    edges: np.ndarray  # increasing interior cut points on degree
    bin_of_gene: np.ndarray  # (N,) int, -1 for degree 0
    members: list[np.ndarray]  # gene indices per bin

    @classmethod
    def from_degrees(cls, degrees: np.ndarray, n_bins: int = 10) -> DegreeBins:
        degrees = np.asarray(degrees, dtype=np.int64)
        pos = degrees[degrees > 0]
        if pos.size == 0:
            return cls(
                degrees=degrees,
                edges=np.zeros(0),
                bin_of_gene=np.full(degrees.size, -1),
                members=[],
            )
        qs = np.linspace(0, 100, n_bins + 1)[1:-1]
        edges = np.unique(np.percentile(pos, qs, method="nearest")).astype(np.float64)
        bins = np.full(degrees.size, -1, dtype=np.int64)
        mask = degrees > 0
        bins[mask] = np.searchsorted(edges, degrees[mask], side="left")
        n_b = int(bins.max()) + 1 if mask.any() else 0
        members = [np.flatnonzero(bins == b) for b in range(n_b)]
        return cls(degrees=degrees, edges=edges, bin_of_gene=bins, members=members)

    @property
    def n_bins(self) -> int:
        return len(self.members)

    def counts_for(self, query: np.ndarray) -> np.ndarray:
        q = np.unique(np.asarray(query, dtype=np.int64))
        b = self.bin_of_gene[q]
        if (b < 0).any():
            raise ValueError("query contains genes with degree 0 (not in the tier network)")
        return np.bincount(b, minlength=self.n_bins).astype(np.int64)

    def bin_ranges(self) -> list[tuple[int, int]]:
        """(min degree, max degree) per bin."""
        out = []
        for m in self.members:
            d = self.degrees[m]
            out.append((int(d.min()), int(d.max())) if d.size else (0, 0))
        return out

    def diagnostics(self, query: np.ndarray) -> dict[str, object]:
        """Feasibility / randomisation diagnostics for one query."""
        q = np.unique(np.asarray(query, dtype=np.int64))
        counts = self.counts_for(q)
        pops = np.array([m.size for m in self.members], dtype=np.int64)
        occ = np.where(pops > 0, counts / np.maximum(pops, 1), 0.0)
        # Expected fraction of a resample that coincides with the query's own genes.
        expected_self = float((counts**2 / np.maximum(pops, 1)).sum() / max(q.size, 1))
        ranges = self.bin_ranges()
        return {
            "n_query": int(q.size),
            "n_bins": self.n_bins,
            "bin_degree_ranges": ranges,
            "bin_populations": pops.tolist(),
            "query_counts_per_bin": counts.tolist(),
            "max_bin_occupancy": float(occ.max()) if occ.size else 0.0,
            "n_bins_nearly_fixed": int((occ > 0.5).sum()),
            "expected_self_overlap_fraction": expected_self,
            "mean_log10_degree_query": float(np.log10(self.degrees[q]).mean())
            if q.size
            else float("nan"),
        }


def degree_matched_samples(
    rng: np.random.Generator, bins: DegreeBins, query: np.ndarray, reps: int
) -> list[np.ndarray]:
    """Size-matched draws preserving the query's count in every degree bin (no replacement)."""
    counts = bins.counts_for(query)
    out: list[np.ndarray] = []
    for _ in range(reps):
        parts = [
            rng.choice(bins.members[b], size=int(c), replace=False)
            for b, c in enumerate(counts)
            if c > 0
        ]
        out.append(np.sort(np.concatenate(parts)) if parts else np.zeros(0, dtype=np.int64))
    return out


def achieved_degree_balance(
    bins: DegreeBins, query: np.ndarray, samples: Sequence[np.ndarray]
) -> dict[str, float]:
    """Mean log10 degree of samples vs. the query (per-bin matching is exact by construction)."""
    q = np.unique(np.asarray(query, dtype=np.int64))
    qd = float(np.log10(bins.degrees[q]).mean()) if q.size else float("nan")
    sd = [float(np.log10(bins.degrees[s]).mean()) for s in samples if s.size]
    distinct = np.unique(np.concatenate(list(samples))).size if samples else 0
    return {
        "query_mean_log10_degree": qd,
        "samples_mean_log10_degree": float(np.mean(sd)) if sd else float("nan"),
        "samples_sd_log10_degree": float(np.std(sd)) if sd else float("nan"),
        "distinct_genes_sampled": int(distinct),
        "distinct_genes_over_query_size": float(distinct / max(q.size, 1)),
    }
