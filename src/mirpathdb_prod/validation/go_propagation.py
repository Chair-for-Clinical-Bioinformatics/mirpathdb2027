"""GO ``all`` vs ``direct`` propagation sensitivity -- pure numerical core.

Production GO gene sets are built from the PROPAGATED ``GOALL``/``ONTOLOGYALL``
columns of ``org.*.eg.db`` (``config.yaml`` ``go_propagation: all``, the GO
true-path rule: a gene annotated to a child term also counts toward every
ancestor). The ``direct`` mode (``GO``/``ONTOLOGY``) keeps only the gene's own
annotations. This module quantifies how much that choice changes (a) the gene-set
STRUCTURE and (b) the ENRICHMENT CALLS, so the owner can see the sensitivity of
the published database to the propagation decision WITHOUT changing the
production default.

Everything here is pure Python over plain dicts / lists so the metrics are
testable on tiny synthetic fixtures; the R / filesystem edges (building gene sets,
running clusterProfiler, dumping the GO DAG, locating production outputs) live in
:mod:`mirpathdb_prod.validation.go_propagation_driver`.

Definitions (also written verbatim into the human report)
---------------------------------------------------------
* **Association** = one ``(term, gene SYMBOL)`` pair of a family's TERM2GENE.
  ``n_assoc_propagated = n_assoc_all - n_assoc_direct``; an association present in
  ``direct`` but NOT in ``all`` is an anomaly (the true-path rule guarantees
  ``direct`` is a subset of ``all``) and is REPORTED, never dropped silently.
* **GO depth** of a term = the MINIMUM number of parent edges from the ontology
  root (``GO:0008150`` BP, ``GO:0005575`` CC, ``GO:0003674`` MF) to the term, over
  ALL relationship types in ``GO.db``'s ``GO{BP,CC,MF}PARENTS`` (``is_a``,
  ``part_of``, ``regulates``, ``positively_regulates``, ``negatively_regulates``).
  The root has depth 0; a term unreachable from its root has depth ``None`` and
  is counted separately (``unreachable``). Computed by :func:`compute_depths`, a
  deterministic BFS over the edge dump.
* **Cell** = one ``(organism, mirna, family, term, tier)`` ORA test. A cell is
  TESTED when clusterProfiler ``enricher`` scored it (overlap >= 1). A cell that
  was not tested is NEVER treated as ``p = 1``: it is simply absent from the
  :class:`CellTable`.
* **Tested proxy for production matrices**: the production matrix CSV stores only
  the raw ``p`` and writes ``1.0`` for untested cells, which is indistinguishable
  from a genuinely computed ``p == 1.0`` (the documented M1 lossiness, exactly as
  ``scripts/export_frontend_contracts.py`` handles it). So when a production
  matrix is reused, ``tested := p < 1.0`` and ``p_adjust`` is RE-COMPUTED by
  step-up BH per miRNA over its tested terms (:func:`bh_adjust`) -- EXCEPT for
  cells present in the authoritative significant rows (``p_adjust <= 0.25``),
  whose ``p_adjust`` + ``overlap`` are taken verbatim. The authoritative rows of
  a production cell live in TWO files once the contract export ran: the shipped
  ``enrichment_significant/<org>/<FAM>_<TIER>.tsv`` (foreign-species miRNAs
  removed by the export-time species policy) and the verbatim quarantine
  ``species_excluded/<org>/<FAM>_<TIER>.tsv``; BOTH are overlaid (the matrices
  themselves keep every miRNA, and so does the isolated ``direct`` side).
* **padj floor on a reused production side**: a cell absent from the overlaid
  authoritative rows is KNOWN to have an authoritative ``p_adjust > 0.25`` only
  when those rows are COMPLETE for its miRNA; then its significance calls use
  ``max(recomputed BH, nextafter(0.25))`` so no threshold ``<= 0.25`` can call it
  significant (:meth:`CellTable.effective_padj`). Completeness is decided PER
  miRNA from the data at hand, never from contract bookkeeping (the export's
  ``rows_excluded_species`` count drops to 0 on an idempotent re-export): the
  shipped TSV always vouches for the miRNAs that PASS the export-time species
  policy (:func:`mirna_retained_by_species_policy`, a pinned mirror of
  ``scripts/export_frontend_contracts.retain_row``); foreign-species /
  unparseable miRNAs are vouched for only by the quarantine or by a TSV that
  still carries such rows (never filtered). Otherwise those miRNAs are EXEMPT
  from the floor and their non-overlaid cells use the recomputed BH as is
  (approximate, counted in ``n_padj_floor_exempt``). The floor is dropped
  entirely (``None``) when the row files are inconsistent with the build
  manifest / matrix / emitted contract (:func:`production_padj_floor`). The padj
  correlations / delta quantiles ALWAYS use the raw recomputed values, never
  the floored ones.
* **Significance**: a cell is SIGNIFICANT at threshold ``t`` when
  ``p_adjust <= t``. The primary threshold is the production criterion
  (``reporting.significant.FDR_SIGNIFICANT`` = 0.05, overridable by config / CLI,
  see :func:`resolve_significance_threshold`); any secondary thresholds are
  labelled ``"secondary sensitivity threshold, NOT the production criterion"``.
* **-log10 transform**: ``p`` (or ``p_adjust``) is clamped to ``10**-cap``
  (``cap = 350``) before the log so a zero never becomes ``inf``.

Ranking formulas (verbatim, no manual cherry-picking)
-----------------------------------------------------
* ``go_terms_gained.tsv``: sort by ``gained`` desc, then
  ``median_neglog10_padj_shift`` desc, then term id.
* ``go_terms_lost.tsv``: sort by ``lost`` desc, then
  ``median_neglog10_padj_shift`` asc, then term id.
* "largest median shift": top by ``abs(median_neglog10_padj_shift)``.
* "largest membership expansion": top by ``size_delta``.
* "broad parent terms": terms with ``depth <= 2`` sorted by ``gained`` desc.
* per-organism/tier affected summary: totals of gained / lost by
  ``(organism, tier)`` over the compared families.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
from array import array
from collections import deque
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any

from mirpathdb_prod.reporting.significant import (
    FDR_CUTOFF,
    FDR_SIGNIFICANT,
    SIGNIFICANT_HEADER,
)

# --------------------------------------------------------------------------- #
# Constants / definitions
# --------------------------------------------------------------------------- #
# GO ontology roots (depth 0) per family token.
GO_ROOTS: dict[str, str] = {
    "GO_BP": "GO:0008150",
    "GO_CC": "GO:0005575",
    "GO_MF": "GO:0003674",
}
# GO.db ontology code (edge dump column) -> family token.
ONTOLOGY_TO_FAMILY: dict[str, str] = {"BP": "GO_BP", "CC": "GO_CC", "MF": "GO_MF"}

GO_DEPTH_DEFINITION = (
    "minimum number of parent edges from the ontology root (GO:0008150 BP, "
    "GO:0005575 CC, GO:0003674 MF) to the term, over ALL relationship types in "
    "GO.db's GO{BP,CC,MF}PARENTS (is_a, part_of, regulates, positively_regulates, "
    "negatively_regulates); root depth 0; terms unreachable from the root have "
    "depth None and are counted separately"
)
UNIVERSE_DEFINITION = (
    "union of the family's TERM2GENE member symbols; clusterProfiler intersects "
    "the query with it"
)
MULTIPLE_TESTING = 'BH (clusterProfiler pAdjustMethod="BH", per gene set)'
TESTED_DEFINITION_AUTHORITATIVE = (
    "every cell scored by clusterProfiler enricher (overlap >= 1), read from the "
    "full long-form cells TSV with exact p, p_adjust, overlap"
)
TESTED_DEFINITION_PRODUCTION = (
    "proxy: matrix p < 1.0 (the production matrix writes 1.0 for untested cells, "
    "indistinguishable from a computed p == 1.0 -- documented M1 lossiness)"
)

# padj provenance of a CellTable.
PADJ_AUTHORITATIVE = "authoritative"
PADJ_RECOMPUTED = "recomputed_bh_from_matrix"

# Why a reused production cell did (not) get the > 0.25 padj floor.
FLOOR_REASON_QUARANTINE = (
    "significant TSV + species_excluded quarantine overlaid: authoritative rows complete"
)
FLOOR_REASON_UNFILTERED = (
    "significant TSV overlaid; it still carries foreign-species rows, so the export-time "
    "species policy never filtered it: authoritative rows complete"
)
FLOOR_REASON_NATIVE_ONLY = (
    "significant TSV overlaid; every matrix miRNA passes the species policy, so the "
    "export-time filter cannot have removed any of its rows: authoritative rows complete"
)
FLOOR_REASON_FOREIGN_EXEMPT = (
    "significant TSV overlaid WITHOUT the species_excluded quarantine: floor applied to the "
    "miRNAs that pass the species policy only; foreign-species / unparseable miRNAs use the "
    "recomputed BH as is (approximate, see n_padj_floor_exempt)"
)
FLOOR_REASON_INCONSISTENT = (
    "NO floor: the authoritative row files disagree with the build manifest / matrix / "
    "emitted contract -- recomputed BH used as is (approximate)"
)
FLOOR_REASON_NO_TSV = "NO floor: no significant TSV -- recomputed BH used as is (approximate)"

# Export-time species policy, mirrored from scripts/export_frontend_contracts.py
# (``mirna_prefix`` / ``retain_row``; pinned by tests). The contract export moves
# every significant row whose miRNA FAILS this predicate from the shipped TSV to
# the species_excluded quarantine, so a reused production cell must know which
# matrix miRNAs the shipped TSV alone can no longer vouch for.
VIRAL_PREFIXES: frozenset[str] = frozenset({"ebv", "kshv", "hcmv"})
VIRAL_HOST_ORGANISM = "hsa"

# Secondary thresholds: a sensitivity check alongside the production criterion.
DEFAULT_SECONDARY_FDRS: tuple[float, ...] = (0.005,)
SECONDARY_FDR_LABEL = "secondary sensitivity threshold, NOT the production criterion"


def threshold_above_floor_warning(thresholds: Sequence[float]) -> str | None:
    """Warning text when a threshold exceeds :data:`FDR_CUTOFF` on a reused production side.

    The ``> 0.25`` floor makes calls exact only up to 0.25; above it every
    non-overlaid cell is called on the recomputed BH, whose ``n`` is the matrix's
    ``p < 1.0`` proxy (smaller than clusterProfiler's), so those calls are
    approximate. ``None`` when every threshold is ``<= FDR_CUTOFF``.
    """
    above = sorted({float(t) for t in thresholds if float(t) > FDR_CUTOFF})
    if not above:
        return None
    return (
        f"threshold(s) {above} exceed FDR_CUTOFF = {FDR_CUTOFF}: on a reused production "
        f"all side the significance calls above {FDR_CUTOFF} use the recomputed BH "
        f"(approximate; the matrix's p < 1.0 tested proxy shrinks BH's n)"
    )

# -log10 floor for p == 0 (mirrors validation.v2_reproduction._NEG_LOG10_CAP).
NEG_LOG10_CAP = 350.0

# Where the comparison threshold may come from, in precedence order.
THRESHOLD_SOURCE_CLI = "cli:--fdr"
THRESHOLD_SOURCE_CONFIG = "config:validation.go_propagation.fdr_significant"
THRESHOLD_SOURCE_CODE = "code:mirpathdb_prod.reporting.significant.FDR_SIGNIFICANT"

RANKING_FORMULAS: dict[str, str] = {
    "go_terms_gained": (
        "sorted by gained desc, then median_neglog10_padj_shift desc, then term id"
    ),
    "go_terms_lost": (
        "sorted by lost desc, then median_neglog10_padj_shift asc, then term id"
    ),
    "largest_median_shift": "top by abs(median_neglog10_padj_shift)",
    "largest_membership_expansion": "top by size_delta (= size_all - size_direct)",
    "broad_parent_terms": "terms with depth <= 2 sorted by gained desc",
    "affected_summary": "totals of gained / lost by (organism, tier) over the compared families",
}

# Term-status vocabulary of the per-term structural table.
TERM_ONLY_ALL = "only_all"
TERM_ONLY_DIRECT = "only_direct"
TERM_BOTH_IDENTICAL = "both_identical"
TERM_BOTH_CHANGED = "both_changed_membership"

DEPTH_UNREACHABLE = "unreachable"

# Column orders of the deterministic TSV artefacts.
STRUCTURAL_TERMS_HEADER = (
    "organism",
    "family",
    "term",
    "name",
    "size_direct",
    "size_all",
    "size_delta",
    "depth",
    "status",
)
STRUCTURAL_TSV_HEADER = (
    "organism",
    "family",
    "n_genes_direct",
    "n_genes_all",
    "n_terms_direct",
    "n_terms_all",
    "n_assoc_direct",
    "n_assoc_all",
    "n_assoc_propagated",
    "pct_assoc_propagated",
    "n_direct_only_assoc",
    "term_size_median_direct",
    "term_size_median_all",
    "n_terms_only_in_all",
    "n_terms_only_in_direct",
    "n_terms_both_changed_membership",
    "n_terms_both_identical",
    "n_terms_depth_unreachable",
    "term_size_min_direct",
    "term_size_min_all",
    "term_size_q25_direct",
    "term_size_q25_all",
    "term_size_median_direct_q",
    "term_size_median_all_q",
    "term_size_q75_direct",
    "term_size_q75_all",
    "term_size_max_direct",
    "term_size_max_all",
    "depth_distribution_gained_terms",
    "depth_distribution_all_terms",
)
TERM_RECORD_HEADER = (
    "organism",
    "family",
    "term",
    "name",
    "n_sig_direct",
    "n_sig_all",
    "gained",
    "lost",
    "median_neglog10_padj_shift",
    "n_shared_tested",
    "size_direct",
    "size_all",
    "size_delta",
    "depth",
)

# How many rows the human-sized ranked tables keep (full versions are gzipped).
TOP_N_RANKED = 200


# --------------------------------------------------------------------------- #
# Small numeric helpers
# --------------------------------------------------------------------------- #
def bh_adjust(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg adjust a list of p-values (returns same order).

    Standard step-up BH with the monotone (cumulative-min from the largest rank)
    enforcement: ``padj_(i) = min_{j>=i} ( p_(j) * n / rank_j )``, clamped to 1.
    ``n`` = ``len(pvalues)`` = the number of TESTED terms for this miRNA. This is
    the same function as ``scripts/export_frontend_contracts.bh_adjust`` (pinned
    by a test) and reproduces R's ``p.adjust(method="BH")``.
    """
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    adj = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        val = pvalues[idx] * n / rank
        prev = min(prev, val)
        adj[idx] = min(prev, 1.0)
    return adj


def neg_log10(p: float, *, cap: float = NEG_LOG10_CAP) -> float:
    """``-log10(p)`` with ``p`` clamped to ``10**-cap`` so ``p == 0`` stays finite.

    Implemented on the log side (``min(-log10(p), cap)``; ``p <= 0`` -> ``cap``)
    because ``10.0 ** -350`` underflows to ``0.0`` in a double, which would make a
    p-side clamp a no-op and let ``log10(0)`` raise.
    """
    cap = float(cap)
    if p <= 0.0:
        return cap
    return min(-math.log10(p), cap)


def quantiles(values: Sequence[float], probs: Sequence[float]) -> list[float]:
    """Linear-interpolation quantiles (numpy's default) of an unsorted sequence."""
    xs = sorted(values)
    n = len(xs)
    if n == 0:
        raise ValueError("quantiles of an empty sequence")
    out: list[float] = []
    for q in probs:
        pos = q * (n - 1)
        lo = math.floor(pos)
        hi = math.ceil(pos)
        if lo == hi:
            out.append(float(xs[lo]))
        else:
            frac = pos - lo
            out.append(float(xs[lo] + (xs[hi] - xs[lo]) * frac))
    return out


def _five_number(values: Sequence[float]) -> dict[str, float] | None:
    if not values:
        return None
    q = quantiles(values, (0.0, 0.25, 0.5, 0.75, 1.0))
    return {"min": q[0], "q25": q[1], "median": q[2], "q75": q[3], "max": q[4]}


def _seven_number(values: Sequence[float]) -> dict[str, float] | None:
    if not values:
        return None
    q = quantiles(values, (0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0))
    return {
        "min": q[0],
        "q05": q[1],
        "q25": q[2],
        "median": q[3],
        "q75": q[4],
        "q95": q[5],
        "max": q[6],
    }


def _median(values: Sequence[float]) -> float | None:
    if not values:
        return None
    return quantiles(values, (0.5,))[0]


def _spearman(x: Sequence[float], y: Sequence[float]) -> float | None:
    """Spearman rho; ``None`` for < 3 points or a constant side (undefined)."""
    if len(x) < 3 or len(x) != len(y):
        return None
    if min(x) == max(x) or min(y) == max(y):
        return None
    from scipy.stats import spearmanr

    rho = spearmanr(x, y)[0]
    if rho is None or (isinstance(rho, float) and math.isnan(rho)):
        return None
    return float(rho)


def _pearson(x: Sequence[float], y: Sequence[float]) -> float | None:
    """Pearson r; ``None`` for < 3 points or a constant side (undefined)."""
    if len(x) < 3 or len(x) != len(y):
        return None
    if min(x) == max(x) or min(y) == max(y):
        return None
    import numpy as np

    r = float(np.corrcoef(np.asarray(x, dtype=float), np.asarray(y, dtype=float))[0, 1])
    return None if math.isnan(r) else r


def _fmt(value: object) -> str:
    """TSV cell formatting: floats via ``repr`` (round-trip exact), ``None`` empty."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _depth_key(depth: int | None) -> str:
    return DEPTH_UNREACHABLE if depth is None else str(depth)


def _histogram_str(hist: Mapping[str, int]) -> str:
    """Serialise a depth histogram as ``0:1;1:5;unreachable:2`` (numeric order)."""

    def _sort_key(k: str) -> tuple[int, int]:
        return (1, 0) if k == DEPTH_UNREACHABLE else (0, int(k))

    return ";".join(f"{k}:{hist[k]}" for k in sorted(hist, key=_sort_key))


# --------------------------------------------------------------------------- #
# 2a. GO depth + structural comparison
# --------------------------------------------------------------------------- #
def compute_depths(edges: Iterable[tuple[str, str]], root: str) -> dict[str, int]:
    """Minimum edge distance from ``root`` for every term reachable from it.

    ``edges`` are ``(child, parent)`` pairs over ALL relationship types. The BFS
    walks parent -> child from the root; a term reached through several parents
    keeps the minimum distance. Terms not reachable from ``root`` are absent from
    the result (callers map them to ``None`` / ``unreachable``). Deterministic:
    BFS distances do not depend on edge order.
    """
    children: dict[str, list[str]] = {}
    for child, parent in edges:
        children.setdefault(parent, []).append(child)
    depth: dict[str, int] = {root: 0}
    queue: deque[str] = deque([root])
    while queue:
        node = queue.popleft()
        d = depth[node] + 1
        for ch in children.get(node, ()):
            if ch not in depth:
                depth[ch] = d
                queue.append(ch)
    return depth


def read_go_parent_edges(path: str | Path) -> dict[str, list[tuple[str, str]]]:
    """Read the R edge dump (``child, parent, relationship, ontology``) by family.

    The ``all``/root pseudo-parent that GO.db reports for each ontology root is
    dropped so the root really is the BFS origin.
    """
    out: dict[str, list[tuple[str, str]]] = {fam: [] for fam in GO_ROOTS}
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fam = ONTOLOGY_TO_FAMILY.get((row.get("ontology") or "").strip())
            child = (row.get("child") or "").strip()
            parent = (row.get("parent") or "").strip()
            if fam is None or not child or not parent or parent == "all":
                continue
            out[fam].append((child, parent))
    return out


def depths_by_family(
    edges_by_family: Mapping[str, Iterable[tuple[str, str]]],
) -> dict[str, dict[str, int]]:
    """Run :func:`compute_depths` from each family's root."""
    return {
        fam: compute_depths(edges_by_family.get(fam, ()), root) for fam, root in GO_ROOTS.items()
    }


@dataclass(frozen=True)
class TermSize:
    """Per-term membership sizes under both modes (one row of the structural table)."""

    term: str
    size_direct: int
    size_all: int
    depth: int | None
    status: str

    @property
    def size_delta(self) -> int:
        return self.size_all - self.size_direct


@dataclass
class StructuralStats:
    """Structural ``direct`` vs ``all`` comparison of one (organism, family)."""

    organism: str
    family: str
    n_genes_direct: int
    n_genes_all: int
    n_terms_direct: int
    n_terms_all: int
    n_assoc_direct: int
    n_assoc_all: int
    n_assoc_propagated: int
    pct_assoc_propagated: float | None
    n_direct_only_assoc: int
    term_size_median_direct: float | None
    term_size_median_all: float | None
    term_size_quantiles_direct: dict[str, float] | None
    term_size_quantiles_all: dict[str, float] | None
    n_terms_only_in_all: int
    n_terms_only_in_direct: int
    n_terms_both_changed_membership: int
    n_terms_both_identical: int
    n_terms_depth_unreachable: int
    depth_distribution_gained_terms: dict[str, int]
    depth_distribution_all_terms: dict[str, int]
    terms: list[TermSize] = field(default_factory=list, repr=False)

    def to_dict(self) -> dict[str, Any]:
        """Summary dict (per-term rows excluded; they go to the gzipped TSV)."""
        out = {f.name: getattr(self, f.name) for f in fields(self) if f.name != "terms"}
        return out

    def to_row(self) -> dict[str, str]:
        """Flat TSV row (quantiles + histograms serialised)."""
        d = self.to_dict()
        qd = d.pop("term_size_quantiles_direct") or {}
        qa = d.pop("term_size_quantiles_all") or {}
        gained = d.pop("depth_distribution_gained_terms")
        all_terms = d.pop("depth_distribution_all_terms")
        row = {k: _fmt(v) for k, v in d.items()}
        for key in ("min", "q25", "q75", "max"):
            row[f"term_size_{key}_direct"] = _fmt(qd.get(key))
            row[f"term_size_{key}_all"] = _fmt(qa.get(key))
        # The median is already a top-level column; the quantile block repeats it
        # under a distinct name so the five-number summary stays contiguous.
        row["term_size_median_direct_q"] = _fmt(qd.get("median"))
        row["term_size_median_all_q"] = _fmt(qa.get("median"))
        row["depth_distribution_gained_terms"] = _histogram_str(gained)
        row["depth_distribution_all_terms"] = _histogram_str(all_terms)
        return row


def structural_comparison(
    direct: Mapping[str, set[str]],
    all_: Mapping[str, set[str]],
    *,
    depth: Mapping[str, int] | None,
    organism: str = "",
    family: str = "",
) -> StructuralStats:
    """Compare the ``direct`` and ``all`` TERM2GENE of one (organism, family).

    Args:
        direct: ``{go_id: set(symbol)}`` built with propagation ``direct``.
        all_: ``{go_id: set(symbol)}`` built with propagation ``all``.
        depth: ``{go_id: depth}`` from :func:`compute_depths` (terms absent from
            it are ``unreachable``); ``None`` marks every term unreachable.
        organism / family: labels carried into the artefacts.
    """
    depth = depth or {}
    genes_direct: set[str] = set()
    genes_all: set[str] = set()
    n_assoc_direct = 0
    n_assoc_all = 0
    n_direct_only_assoc = 0
    for term, genes in direct.items():
        genes_direct |= genes
        n_assoc_direct += len(genes)
        n_direct_only_assoc += len(genes - all_.get(term, set()))
    for genes in all_.values():
        genes_all |= genes
        n_assoc_all += len(genes)

    terms: list[TermSize] = []
    gained_hist: dict[str, int] = {}
    all_hist: dict[str, int] = {}
    n_only_all = n_only_direct = n_changed = n_identical = 0
    n_unreachable = 0
    for term in sorted(set(direct) | set(all_)):
        sd = len(direct.get(term, ()))
        sa = len(all_.get(term, ()))
        d = depth.get(term)
        if term not in direct:
            status = TERM_ONLY_ALL
            n_only_all += 1
            gained_hist[_depth_key(d)] = gained_hist.get(_depth_key(d), 0) + 1
        elif term not in all_:
            status = TERM_ONLY_DIRECT
            n_only_direct += 1
        elif direct[term] == all_[term]:
            status = TERM_BOTH_IDENTICAL
            n_identical += 1
        else:
            status = TERM_BOTH_CHANGED
            n_changed += 1
        if term in all_:
            all_hist[_depth_key(d)] = all_hist.get(_depth_key(d), 0) + 1
            if d is None:
                n_unreachable += 1
        terms.append(TermSize(term=term, size_direct=sd, size_all=sa, depth=d, status=status))

    sizes_direct = [len(g) for g in direct.values()]
    sizes_all = [len(g) for g in all_.values()]
    propagated = n_assoc_all - n_assoc_direct
    return StructuralStats(
        organism=organism,
        family=family,
        n_genes_direct=len(genes_direct),
        n_genes_all=len(genes_all),
        n_terms_direct=len(direct),
        n_terms_all=len(all_),
        n_assoc_direct=n_assoc_direct,
        n_assoc_all=n_assoc_all,
        n_assoc_propagated=propagated,
        pct_assoc_propagated=(100.0 * propagated / n_assoc_all) if n_assoc_all else None,
        n_direct_only_assoc=n_direct_only_assoc,
        term_size_median_direct=_median(sizes_direct),
        term_size_median_all=_median(sizes_all),
        term_size_quantiles_direct=_five_number(sizes_direct),
        term_size_quantiles_all=_five_number(sizes_all),
        n_terms_only_in_all=n_only_all,
        n_terms_only_in_direct=n_only_direct,
        n_terms_both_changed_membership=n_changed,
        n_terms_both_identical=n_identical,
        n_terms_depth_unreachable=n_unreachable,
        depth_distribution_gained_terms=gained_hist,
        depth_distribution_all_terms=all_hist,
        terms=terms,
    )


# --------------------------------------------------------------------------- #
# 2b. Cell tables (one (organism, family, tier) ORA result) + writers/readers
# --------------------------------------------------------------------------- #
@dataclass
class CellTable:
    """All TESTED cells of one (organism, family, tier) under one propagation mode.

    ``p`` / ``padj`` / ``overlap`` are nested ``{mirna: {term: value}}`` maps (the
    nested form shares the miRNA / term string objects and keeps a 40M-cell
    production matrix in memory without a tuple per cell). ``p`` holds EVERY tested
    cell; ``padj`` and ``overlap`` may be sparse (a missing key = unknown).

    Attributes:
        padj_source: :data:`PADJ_AUTHORITATIVE` (from the full long-form cells
            TSV) or :data:`PADJ_RECOMPUTED` (step-up BH over the matrix's tested
            proxy, with authoritative values overlaid where the significant TSV
            has the cell).
        tested_definition: Human-readable statement of what "tested" means here.
        nonauthoritative_padj_floor: For recomputed tables whose overlaid
            authoritative rows are COMPLETE: cells absent from them are KNOWN to
            have an authoritative ``p_adjust`` strictly above this floor (0.25);
            see :meth:`effective_padj` and :func:`production_padj_floor`.
        padj_floor_exempt: miRNAs whose authoritative rows the overlaid files may
            LACK (the export-time species policy moved them to a quarantine that
            was not overlaid): their non-overlaid cells bypass the floor and use
            the recomputed padj as is.
        n_padj_floor_exempt: Number of such (approximate) cells.
        n_authoritative_exempt_rows: Overlaid rows belonging to miRNAs the
            exemption predicate flagged -- when > 0 the files evidently DO cover
            those miRNAs and nothing is exempt.
        significant_paths: The authoritative row files overlaid (significant TSV
            and, when present, the species_excluded quarantine).
    """

    organism: str
    family: str
    tier: str
    p: dict[str, dict[str, float]]
    padj: dict[str, dict[str, float]]
    overlap: dict[str, dict[str, int]]
    padj_source: str
    tested_definition: str
    nonauthoritative_padj_floor: float | None = None
    source_path: str = ""
    n_padj_authoritative: int = 0
    n_padj_recomputed: int = 0
    n_significant_rows_unmatched: int = 0
    significant_paths: tuple[str, ...] = ()
    padj_floor_exempt: frozenset[str] = frozenset()
    n_padj_floor_exempt: int = 0
    n_authoritative_exempt_rows: int = 0

    @property
    def n_tested(self) -> int:
        return sum(len(v) for v in self.p.values())

    def padj_is_authoritative(self, mirna: str, term: str) -> bool:
        """True when the cell's padj came straight from clusterProfiler.

        For a recomputed table the authoritative cells are exactly those overlaid
        from the significant TSV, which are exactly the cells with a known overlap.
        """
        if self.padj_source == PADJ_AUTHORITATIVE:
            return True
        return term in self.overlap.get(mirna, ())

    def effective_padj(self, mirna: str, term: str) -> float | None:
        """padj used for significance calls (``None`` = unknown, never significant).

        Authoritative cells return their padj. A recomputed cell whose authoritative
        padj is known to exceed ``nonauthoritative_padj_floor`` returns
        ``max(recomputed, nextafter(floor))`` so no threshold ``<= floor`` can call
        it significant (the recomputed BH is only an approximation above the floor).
        A miRNA in ``padj_floor_exempt`` is not vouched for by the overlaid rows,
        so its recomputed cells are returned as is.
        """
        val = self.padj.get(mirna, {}).get(term)
        if val is None:
            return None
        if self.padj_source == PADJ_AUTHORITATIVE or term in self.overlap.get(mirna, ()):
            return val
        floor = self.nonauthoritative_padj_floor
        if floor is None or mirna in self.padj_floor_exempt:
            return val
        return max(val, math.nextafter(floor, math.inf))

    def is_significant(self, mirna: str, term: str, fdr: float) -> bool:
        eff = self.effective_padj(mirna, term)
        return eff is not None and eff <= fdr


def write_matrix(
    out_path: str | Path,
    mirnas: list[str],
    terms: list[str],
    pvals: Mapping[str, Mapping[str, tuple[float, float, int]]],
) -> tuple[int, int]:
    """Write a (miRNA x term) p-value matrix CSV in the production ``_write_matrix`` format.

    Header ``miRNA,<term>,...``; one row per miRNA with ``repr(float(p))`` and the
    literal ``1.0`` where the term was not scored for that miRNA (the same lossy
    sentinel the production runners write). Returns ``(n_rows, n_cols)``.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    col_idx = {t: j for j, t in enumerate(terms)}
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["miRNA", *terms])
        for mirna in mirnas:
            row = [1.0] * len(terms)
            for term, (p, _padj, _ov) in pvals.get(mirna, {}).items():
                j = col_idx.get(term)
                if j is not None:
                    row[j] = p
            writer.writerow([mirna, *(repr(float(v)) for v in row)])
    return len(mirnas), len(terms)


class _GzipTextWriter:
    """Text writer over a byte-reproducible gzip container (mtime 0, no filename)."""

    def __init__(self, path: Path) -> None:
        self._raw = path.open("wb")
        self._gz = gzip.GzipFile(filename="", mode="wb", fileobj=self._raw, mtime=0)
        self.text = io.TextIOWrapper(self._gz, encoding="utf-8", newline="")

    def __enter__(self) -> io.TextIOWrapper:
        return self.text

    def __exit__(self, *exc: object) -> None:
        self.text.close()
        self._raw.close()


def _open_text(path: Path, mode: str):
    """Open a possibly-gzipped text file (``.gz`` suffix selects gzip).

    Writes to ``.gz`` go through a fixed-mtime container so identical rows give
    identical bytes.
    """
    if path.suffix == ".gz":
        if mode == "w":
            return _GzipTextWriter(path)
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def write_cells_long_form(
    out_path: str | Path,
    mirnas: Sequence[str],
    pvals: Mapping[str, Mapping[str, tuple[float, float, int]]],
) -> int:
    """Write EVERY tested cell as ``mirna, pathway_id, p_value, p_adjust, overlap``.

    Unlike the production significant subset (``p_adjust <= 0.25``) this keeps all
    scored cells, so a comparison can use exact tested sets and authoritative padj.
    Rows are sorted by (mirna, term); gzip when the path ends in ``.gz``. Returns
    the row count.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with _open_text(out, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(SIGNIFICANT_HEADER)
        for mirna in sorted(mirnas):
            per = pvals.get(mirna, {})
            for term in sorted(per):
                p, padj, ov = per[term]
                writer.writerow([mirna, term, repr(float(p)), repr(float(padj)), int(ov)])
                n += 1
    return n


def load_cell_table_from_long_form(
    path: str | Path, *, organism: str = "", family: str = "", tier: str = ""
) -> CellTable:
    """Load a full long-form cells TSV (see :func:`write_cells_long_form`)."""
    src = Path(path)
    p: dict[str, dict[str, float]] = {}
    padj: dict[str, dict[str, float]] = {}
    overlap: dict[str, dict[str, int]] = {}
    n = 0
    with _open_text(src, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None or tuple(header) != SIGNIFICANT_HEADER:
            raise ValueError(
                f"{src}: unexpected long-form header {header!r}; expected {SIGNIFICANT_HEADER}"
            )
        # csv.reader yields a fresh str per field; interning the term ids shares
        # one object across every miRNA row (tens of millions of cells on LOW).
        interned: dict[str, str] = {}
        for row in reader:
            mirna = row[0]
            term = interned.setdefault(row[1], row[1])
            p.setdefault(mirna, {})[term] = float(row[2])
            padj.setdefault(mirna, {})[term] = float(row[3])
            overlap.setdefault(mirna, {})[term] = int(float(row[4]))
            n += 1
    return CellTable(
        organism=organism,
        family=family,
        tier=tier,
        p=p,
        padj=padj,
        overlap=overlap,
        padj_source=PADJ_AUTHORITATIVE,
        tested_definition=TESTED_DEFINITION_AUTHORITATIVE,
        nonauthoritative_padj_floor=None,
        source_path=str(src),
        n_padj_authoritative=n,
        n_padj_recomputed=0,
    )


def mirna_prefix(name: str) -> str | None:
    """Lowercased organism prefix of a miRNA name (token before the first ``-``).

    ``None`` when there is no ``-`` at all or the leading token is not alphabetic
    (mirrors ``export_frontend_contracts.mirna_prefix``).
    """
    if not name:
        return None
    head, sep, _rest = name.partition("-")
    if not sep:
        return None
    head = head.strip().lower()
    if not head.isalpha():
        return None
    return head


def mirna_retained_by_species_policy(name: str, organism: str) -> bool:
    """True iff the contract export keeps a row with miRNA ``name`` in ``organism``.

    FAIL-CLOSED like ``export_frontend_contracts.retain_row``: an unparseable or
    foreign prefix is excluded; the human-infecting viral prefixes are retained in
    the hsa slice only.
    """
    prefix = mirna_prefix(name)
    if prefix is None:
        return False
    org = organism.strip().lower()
    if prefix == org:
        return True
    return org == VIRAL_HOST_ORGANISM and prefix in VIRAL_PREFIXES


def production_padj_floor(
    *,
    has_significant_tsv: bool,
    quarantine_present: bool,
    n_floor_exempt_mirnas: int = 0,
    n_authoritative_exempt_rows: int = 0,
    inconsistencies: Sequence[str] = (),
) -> tuple[float | None, str]:
    """Decide whether a reused production cell may use the ``> 0.25`` padj floor.

    Inputs are what the loaded :class:`CellTable` and the driver's integrity
    checks observed -- never the emitted contract's ``rows_excluded_species``
    count, which is 0 for every cell after an idempotent re-export and so cannot
    prove that nothing was ever filtered out:

    * no significant TSV -> no floor;
    * any ``inconsistencies`` (row files vs build manifest / matrix / contract
      sha256) -> no floor, the messages are appended to the reason;
    * quarantine overlaid -> floor, rows complete;
    * ``n_floor_exempt_mirnas > 0`` -> floor for the species-policy-retained
      miRNAs only (the exempt ones use the recomputed BH as is);
    * ``n_authoritative_exempt_rows > 0`` -> the files still carry foreign rows,
      i.e. were never filtered -> floor, rows complete;
    * else every matrix miRNA is retained by the policy -> floor, rows complete.

    Returns ``(floor, reason)``; the reason is recorded in the summary / report.
    """
    if not has_significant_tsv:
        return None, FLOOR_REASON_NO_TSV
    if inconsistencies:
        return None, FLOOR_REASON_INCONSISTENT + ": " + "; ".join(inconsistencies)
    if quarantine_present:
        return FDR_CUTOFF, FLOOR_REASON_QUARANTINE
    if n_floor_exempt_mirnas > 0:
        return FDR_CUTOFF, FLOOR_REASON_FOREIGN_EXEMPT
    if n_authoritative_exempt_rows > 0:
        return FDR_CUTOFF, FLOOR_REASON_UNFILTERED
    return FDR_CUTOFF, FLOOR_REASON_NATIVE_ONLY


def load_cell_table_from_production(
    matrix_csv: str | Path,
    significant_tsv: str | Path | Sequence[str | Path] | None,
    *,
    organism: str = "",
    family: str = "",
    tier: str = "",
    padj_floor: float | None = FDR_CUTOFF,
    floor_exempt: Callable[[str], bool] | None = None,
) -> CellTable:
    """Build a :class:`CellTable` from a production matrix (+ authoritative row files).

    The matrix is streamed row by row with :mod:`csv`; only tested cells (proxy
    ``p < 1.0``, see the module docstring) are kept. ``p_adjust`` is recomputed
    per miRNA by :func:`bh_adjust` over its tested terms, then OVERLAID with the
    authoritative ``p_adjust`` + ``overlap`` for every cell present in
    ``significant_tsv`` -- one path or several (the shipped significant TSV plus
    the species_excluded quarantine written by the contract export; all share the
    :data:`SIGNIFICANT_HEADER` layout). Rows that are not tested cells of the
    matrix are counted in ``n_significant_rows_unmatched`` (never silently
    dropped).

    ``padj_floor`` becomes the table's ``nonauthoritative_padj_floor`` when at
    least one authoritative file is supplied. ``floor_exempt(mirna)`` marks the
    matrix miRNAs the supplied files may NOT vouch for (the driver passes the
    negated species policy when the quarantine is absent): their non-overlaid
    cells bypass the floor. The exemption is dropped when the files DO contain
    rows of flagged miRNAs (``n_authoritative_exempt_rows > 0``): the export
    filters a file all-or-nothing, so one such row proves it was never filtered.
    """
    src = Path(matrix_csv)
    if significant_tsv is None:
        sig_paths: tuple[Path, ...] = ()
    elif isinstance(significant_tsv, str | Path):
        sig_paths = (Path(significant_tsv),)
    else:
        sig_paths = tuple(Path(s) for s in significant_tsv)
    p: dict[str, dict[str, float]] = {}
    padj: dict[str, dict[str, float]] = {}
    n_recomputed = 0
    with src.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header or header[0] != "miRNA":
            raise ValueError(f"{src}: expected a 'miRNA,<term>,...' header, got {header!r}")
        terms = header[1:]
        for row in reader:
            if not row:
                continue
            mirna = row[0]
            if len(row) != len(header):
                raise ValueError(f"{src}: row {mirna!r} has {len(row)} cells, header {len(header)}")
            tested_terms: list[str] = []
            tested_p: list[float] = []
            for term, cell in zip(terms, row[1:], strict=True):
                val = float(cell)
                if val < 1.0:
                    tested_terms.append(term)
                    tested_p.append(val)
            if not tested_terms:
                continue
            adj = bh_adjust(tested_p)
            p[mirna] = dict(zip(tested_terms, tested_p, strict=True))
            padj[mirna] = dict(zip(tested_terms, adj, strict=True))
            n_recomputed += len(tested_terms)

    overlap: dict[str, dict[str, int]] = {}
    n_auth = 0
    n_unmatched = 0
    n_exempt_rows = 0
    flagged: dict[str, bool] = {}  # memo of floor_exempt(mirna) over the row files
    for sig in sig_paths:
        with sig.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            if header is None or tuple(header[:5]) != SIGNIFICANT_HEADER:
                raise ValueError(
                    f"{sig}: unexpected header {header!r}; expected {SIGNIFICANT_HEADER}"
                )
            for row in reader:
                mirna, term = row[0], row[1]
                if floor_exempt is not None:
                    flag = flagged.get(mirna)
                    if flag is None:
                        flag = flagged[mirna] = bool(floor_exempt(mirna))
                    if flag:
                        n_exempt_rows += 1
                per = p.get(mirna)
                if per is None or term not in per:
                    n_unmatched += 1
                    continue
                padj[mirna][term] = float(row[3])
                overlap.setdefault(mirna, {})[term] = int(float(row[4]))
                n_auth += 1
    exempt: frozenset[str] = frozenset()
    if floor_exempt is not None and sig_paths and n_exempt_rows == 0:
        exempt = frozenset(m for m in p if floor_exempt(m))
    n_cells_exempt = sum(len(p[m]) - len(overlap.get(m, ())) for m in exempt)
    return CellTable(
        organism=organism,
        family=family,
        tier=tier,
        p=p,
        padj=padj,
        overlap=overlap,
        padj_source=PADJ_RECOMPUTED,
        tested_definition=TESTED_DEFINITION_PRODUCTION,
        nonauthoritative_padj_floor=padj_floor if sig_paths else None,
        source_path=str(src),
        n_padj_authoritative=n_auth,
        n_padj_recomputed=n_recomputed - n_auth,
        n_significant_rows_unmatched=n_unmatched,
        significant_paths=tuple(str(s) for s in sig_paths),
        padj_floor_exempt=exempt,
        n_padj_floor_exempt=n_cells_exempt,
        n_authoritative_exempt_rows=n_exempt_rows,
    )


# --------------------------------------------------------------------------- #
# 2b. Enrichment comparison of one (organism, family, tier)
# --------------------------------------------------------------------------- #
@dataclass
class TermAccumulator:
    """Per-term running totals gathered while comparing one or more cell tables."""

    n_sig_direct: int = 0
    n_sig_all: int = 0
    n_shared_tested: int = 0
    shifts: array = field(default_factory=lambda: array("d"))


@dataclass
class CellComparison:
    """``direct`` vs ``all`` enrichment comparison of one (organism, family, tier)."""

    organism: str
    family: str
    tier: str
    fdr: float
    secondary_fdr_label: str
    padj_source_direct: str
    padj_source_all: str
    padj_floor_direct: float | None
    padj_floor_all: float | None
    n_padj_floor_exempt_direct: int
    n_padj_floor_exempt_all: int
    n_tested_direct: int
    n_tested_all: int
    n_shared_tested: int
    n_only_all_tested: int
    n_only_direct_tested: int
    n_padj_missing_direct: int
    n_padj_missing_all: int
    n_sig_direct: int
    n_sig_all: int
    n_sig_intersection: int
    n_sig_union: int
    jaccard: float | None
    n_sig_only_all: int
    n_sig_only_direct: int
    spearman_p: float | None
    spearman_padj: float | None
    pearson_neglog10_p: float | None
    n_shared_both_padj: int
    n_shared_both_overlap: int
    padj_delta_quantiles: dict[str, float] | None
    overlap_delta_quantiles: dict[str, float] | None
    per_secondary_fdr: dict[str, dict[str, Any]]
    term_acc: dict[str, TermAccumulator] = field(default_factory=dict, repr=False)

    def to_dict(self) -> dict[str, Any]:
        return {f.name: getattr(self, f.name) for f in fields(self) if f.name != "term_acc"}

    def to_row(self) -> dict[str, str]:
        d = self.to_dict()
        pq = d.pop("padj_delta_quantiles") or {}
        oq = d.pop("overlap_delta_quantiles") or {}
        sec = d.pop("per_secondary_fdr")
        row = {k: _fmt(v) for k, v in d.items()}
        for key in ("min", "q05", "q25", "median", "q75", "q95", "max"):
            row[f"padj_delta_{key}"] = _fmt(pq.get(key))
            row[f"overlap_delta_{key}"] = _fmt(oq.get(key))
        for thr in sorted(sec, key=float):
            for k, v in sec[thr].items():
                row[f"secondary_{thr}_{k}"] = _fmt(v)
        return row

    def release_term_data(self) -> None:
        """Drop the per-term accumulators once they have been aggregated."""
        self.term_acc = {}


def _sig_flags(eff: float | None, fdrs: Sequence[float]) -> list[bool]:
    if eff is None:
        return [False] * len(fdrs)
    return [eff <= f for f in fdrs]


def compare_cell_tables(
    direct: CellTable,
    all_: CellTable,
    *,
    fdr: float,
    secondary_fdrs: Sequence[float] = DEFAULT_SECONDARY_FDRS,
    neg_log10_cap: float = NEG_LOG10_CAP,
) -> CellComparison:
    """Compare the tested / significant cells of one (organism, family, tier).

    * tested-set overlap (shared / only-all / only-direct) -- a cell absent from
      a side is NOT tested there, never ``p = 1``;
    * significance sets at ``fdr`` (primary) and each secondary threshold:
      counts, intersection, union, Jaccard, gained (= sig only in ``all``) and
      lost (= sig only in ``direct``) -- calls use :meth:`CellTable.effective_padj`
      (the ``> 0.25`` floor on a reused production side, when valid);
    * Spearman rho on ``-log10(p)`` and on ``-log10(padj)`` over shared tested
      cells (``None`` for < 3 points or a constant side), Pearson r on
      ``-log10(p)``;
    * quantiles of ``padj_all - padj_direct`` (shared cells with both padj known)
      and of ``overlap_all - overlap_direct`` (both overlaps known) -- the padj
      correlations, deltas and per-term shifts use the RAW padj values (never
      the floored ``effective_padj``);
    * per-term accumulators (sig counts, ``-log10(padj)`` shifts) for
      :func:`aggregate_terms`.

    Cells with an unknown padj are excluded from padj statistics and counted in
    ``n_padj_missing_*``. Zero p-values are clamped to ``10**-neg_log10_cap``.
    """
    fdrs = [float(fdr), *(float(f) for f in secondary_fdrs)]
    n_thr = len(fdrs)
    sig_d = [0] * n_thr
    sig_a = [0] * n_thr
    inter = [0] * n_thr
    nlp_d = array("d")
    nlp_a = array("d")
    nlpadj_d = array("d")
    nlpadj_a = array("d")
    padj_delta = array("d")
    ov_delta = array("d")
    term_acc: dict[str, TermAccumulator] = {}
    n_shared = n_only_direct = n_only_all = 0
    n_padj_missing_direct = n_padj_missing_all = 0

    for mirna, terms_d in direct.p.items():
        terms_a = all_.p.get(mirna, {})
        ov_d_row = direct.overlap.get(mirna, {})
        ov_a_row = all_.overlap.get(mirna, {})
        for term, p_d in terms_d.items():
            eff_d = direct.effective_padj(mirna, term)
            if eff_d is None:
                n_padj_missing_direct += 1
            flags_d = _sig_flags(eff_d, fdrs)
            acc = term_acc.get(term)
            if acc is None:
                acc = term_acc[term] = TermAccumulator()
            if flags_d[0]:
                acc.n_sig_direct += 1
            for i, f in enumerate(flags_d):
                if f:
                    sig_d[i] += 1
            p_a = terms_a.get(term)
            if p_a is None:
                n_only_direct += 1
                continue
            n_shared += 1
            acc.n_shared_tested += 1
            nlp_d.append(neg_log10(p_d, cap=neg_log10_cap))
            nlp_a.append(neg_log10(p_a, cap=neg_log10_cap))
            eff_a = all_.effective_padj(mirna, term)
            flags_a = _sig_flags(eff_a, fdrs)
            for i in range(n_thr):
                if flags_d[i] and flags_a[i]:
                    inter[i] += 1
            raw_d = direct.padj.get(mirna, {}).get(term)
            raw_a = all_.padj.get(mirna, {}).get(term)
            if raw_d is not None and raw_a is not None:
                sd = neg_log10(raw_d, cap=neg_log10_cap)
                sa = neg_log10(raw_a, cap=neg_log10_cap)
                nlpadj_d.append(sd)
                nlpadj_a.append(sa)
                padj_delta.append(raw_a - raw_d)
                acc.shifts.append(sa - sd)
            o_d = ov_d_row.get(term)
            o_a = ov_a_row.get(term)
            if o_d is not None and o_a is not None:
                ov_delta.append(float(o_a - o_d))

    for mirna, terms_a in all_.p.items():
        terms_d = direct.p.get(mirna, {})
        for term in terms_a:
            eff_a = all_.effective_padj(mirna, term)
            if eff_a is None:
                n_padj_missing_all += 1
            flags_a = _sig_flags(eff_a, fdrs)
            for i, f in enumerate(flags_a):
                if f:
                    sig_a[i] += 1
            if flags_a[0]:
                acc = term_acc.get(term)
                if acc is None:
                    acc = term_acc[term] = TermAccumulator()
                acc.n_sig_all += 1
            if term not in terms_d:
                n_only_all += 1

    def _set_stats(i: int) -> dict[str, Any]:
        union = sig_d[i] + sig_a[i] - inter[i]
        return {
            "n_sig_direct": sig_d[i],
            "n_sig_all": sig_a[i],
            "n_sig_intersection": inter[i],
            "n_sig_union": union,
            "jaccard": (inter[i] / union) if union else None,
            "n_sig_only_all": sig_a[i] - inter[i],
            "n_sig_only_direct": sig_d[i] - inter[i],
        }

    primary = _set_stats(0)
    secondary = {
        repr(fdrs[i]): {**_set_stats(i), "label": SECONDARY_FDR_LABEL} for i in range(1, n_thr)
    }
    return CellComparison(
        organism=direct.organism or all_.organism,
        family=direct.family or all_.family,
        tier=direct.tier or all_.tier,
        fdr=fdrs[0],
        secondary_fdr_label=SECONDARY_FDR_LABEL,
        padj_source_direct=direct.padj_source,
        padj_source_all=all_.padj_source,
        padj_floor_direct=direct.nonauthoritative_padj_floor,
        padj_floor_all=all_.nonauthoritative_padj_floor,
        n_padj_floor_exempt_direct=direct.n_padj_floor_exempt,
        n_padj_floor_exempt_all=all_.n_padj_floor_exempt,
        n_tested_direct=direct.n_tested,
        n_tested_all=all_.n_tested,
        n_shared_tested=n_shared,
        n_only_all_tested=n_only_all,
        n_only_direct_tested=n_only_direct,
        n_padj_missing_direct=n_padj_missing_direct,
        n_padj_missing_all=n_padj_missing_all,
        spearman_p=_spearman(nlp_d, nlp_a),
        spearman_padj=_spearman(nlpadj_d, nlpadj_a),
        pearson_neglog10_p=_pearson(nlp_d, nlp_a),
        n_shared_both_padj=len(padj_delta),
        n_shared_both_overlap=len(ov_delta),
        padj_delta_quantiles=_seven_number(padj_delta),
        overlap_delta_quantiles=_seven_number(ov_delta),
        per_secondary_fdr=secondary,
        term_acc=term_acc,
        **primary,
    )


# --------------------------------------------------------------------------- #
# 2b. Term-level aggregation + rankings
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TermRecord:
    """Term-level ``direct`` vs ``all`` summary across (tier, mirna) cells."""

    organism: str
    family: str
    term: str
    name: str
    n_sig_direct: int
    n_sig_all: int
    gained: int
    lost: int
    median_neglog10_padj_shift: float | None
    n_shared_tested: int
    size_direct: int | None
    size_all: int | None
    size_delta: int | None
    depth: int | None

    def to_dict(self) -> dict[str, Any]:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    def to_row(self) -> dict[str, str]:
        return {k: _fmt(v) for k, v in self.to_dict().items()}


def aggregate_terms(
    cell_comparisons: Sequence[CellComparison],
    *,
    depth: Mapping[str, int] | None,
    term_names: Mapping[str, str] | None,
    fdr: float,
    term_sizes: Mapping[str, tuple[int, int]] | None = None,
) -> list[TermRecord]:
    """Aggregate per-term significance counts over (organism, family) groups.

    ``gained = n_sig_all - n_sig_direct``, ``lost = max(0, n_sig_direct - n_sig_all)``,
    ``median_neglog10_padj_shift`` = median over shared tested cells of
    ``-log10(padj_all) - (-log10(padj_direct))``. Every comparison must have been
    computed at the same primary ``fdr`` (checked). ``term_sizes`` maps term ->
    ``(size_direct, size_all)`` from the structural comparison. Output is sorted by
    (organism, family, term).
    """
    depth = depth or {}
    term_names = term_names or {}
    term_sizes = term_sizes or {}
    groups: dict[tuple[str, str], dict[str, TermAccumulator]] = {}
    for cmp in cell_comparisons:
        if cmp.fdr != float(fdr):
            raise ValueError(
                f"{cmp.organism}/{cmp.family}/{cmp.tier}: compared at fdr={cmp.fdr}, "
                f"aggregation requested at {fdr}"
            )
        merged = groups.setdefault((cmp.organism, cmp.family), {})
        for term, acc in cmp.term_acc.items():
            tgt = merged.get(term)
            if tgt is None:
                tgt = merged[term] = TermAccumulator()
            tgt.n_sig_direct += acc.n_sig_direct
            tgt.n_sig_all += acc.n_sig_all
            tgt.n_shared_tested += acc.n_shared_tested
            tgt.shifts.extend(acc.shifts)

    out: list[TermRecord] = []
    for (organism, family), merged in sorted(groups.items()):
        for term in sorted(merged):
            acc = merged[term]
            sizes = term_sizes.get(term)
            sd, sa = (sizes if sizes is not None else (None, None))
            out.append(
                TermRecord(
                    organism=organism,
                    family=family,
                    term=term,
                    name=term_names.get(term, ""),
                    n_sig_direct=acc.n_sig_direct,
                    n_sig_all=acc.n_sig_all,
                    gained=acc.n_sig_all - acc.n_sig_direct,
                    lost=max(0, acc.n_sig_direct - acc.n_sig_all),
                    median_neglog10_padj_shift=_median(acc.shifts),
                    n_shared_tested=acc.n_shared_tested,
                    size_direct=sd,
                    size_all=sa,
                    size_delta=(sa - sd) if sizes is not None else None,
                    depth=depth.get(term),
                )
            )
    return out


def _shift_or(rec: TermRecord, default: float) -> float:
    s = rec.median_neglog10_padj_shift
    return default if s is None else s


def rank_gained(records: Sequence[TermRecord]) -> list[TermRecord]:
    """``go_terms_gained`` order: gained desc, median shift desc, term id."""
    return sorted(
        records,
        key=lambda r: (-r.gained, -_shift_or(r, -math.inf), r.term, r.organism, r.family),
    )


def rank_lost(records: Sequence[TermRecord]) -> list[TermRecord]:
    """``go_terms_lost`` order: lost desc, median shift asc, term id."""
    return sorted(
        records,
        key=lambda r: (-r.lost, _shift_or(r, math.inf), r.term, r.organism, r.family),
    )


def top_by_abs_shift(records: Sequence[TermRecord], n: int) -> list[TermRecord]:
    """"largest median shift": top ``n`` by ``abs(median_neglog10_padj_shift)``."""
    with_shift = [r for r in records if r.median_neglog10_padj_shift is not None]
    return sorted(
        with_shift,
        key=lambda r: (-abs(r.median_neglog10_padj_shift or 0.0), r.term, r.organism, r.family),
    )[:n]


def top_by_size_delta(records: Sequence[TermRecord], n: int) -> list[TermRecord]:
    """"largest membership expansion": top ``n`` by ``size_delta``."""
    with_size = [r for r in records if r.size_delta is not None]
    return sorted(with_size, key=lambda r: (-(r.size_delta or 0), r.term, r.organism, r.family))[:n]


def broad_parent_terms(
    records: Sequence[TermRecord], *, max_depth: int = 2, n: int | None = None
) -> list[TermRecord]:
    """"broad parent terms": ``depth <= max_depth`` sorted by ``gained`` desc."""
    broad = [r for r in records if r.depth is not None and r.depth <= max_depth]
    ranked = sorted(broad, key=lambda r: (-r.gained, r.depth or 0, r.term, r.organism, r.family))
    return ranked if n is None else ranked[:n]


def affected_summary(
    cell_comparisons: Sequence[CellComparison],
) -> dict[str, dict[str, dict[str, int]]]:
    """Totals of gained / lost / significant cells by ``(organism, tier)``."""
    out: dict[str, dict[str, dict[str, int]]] = {}
    for cmp in cell_comparisons:
        block = out.setdefault(cmp.organism, {}).setdefault(
            cmp.tier,
            {"gained": 0, "lost": 0, "n_sig_direct": 0, "n_sig_all": 0, "n_families": 0},
        )
        block["gained"] += cmp.n_sig_only_all
        block["lost"] += cmp.n_sig_only_direct
        block["n_sig_direct"] += cmp.n_sig_direct
        block["n_sig_all"] += cmp.n_sig_all
        block["n_families"] += 1
    return out


# --------------------------------------------------------------------------- #
# Threshold resolution
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ThresholdRecord:
    """The resolved significance threshold + where it came from."""

    value: float
    source: str
    production_manifest_value: float | None
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "value": self.value,
            "source": self.source,
            "production_manifest_value": self.production_manifest_value,
            "warnings": list(self.warnings),
        }


def resolve_significance_threshold(
    cfg_raw: Mapping[str, Any] | None,
    cli_value: float | None,
    *,
    provenance_dir: str | Path | None,
) -> ThresholdRecord:
    """Resolve the primary threshold: CLI > config key > code constant.

    Precedence: ``--fdr`` > ``config.yaml`` ``validation.go_propagation.fdr_significant``
    > :data:`FDR_SIGNIFICANT`. When ``<provenance_dir>/production_enrichment_manifest.json``
    exists its ``emitted_contract.fdr_significant_max`` is cross-checked; a mismatch
    is a WARNING recorded in the returned record (not fatal) so the report always
    states which criterion the numbers refer to.
    """
    cfg_raw = cfg_raw or {}
    if cli_value is not None:
        value, source = float(cli_value), THRESHOLD_SOURCE_CLI
    else:
        cfg_val = (
            (cfg_raw.get("validation") or {}).get("go_propagation") or {}
        ).get("fdr_significant")
        if cfg_val is not None:
            value, source = float(cfg_val), THRESHOLD_SOURCE_CONFIG
        else:
            value, source = float(FDR_SIGNIFICANT), THRESHOLD_SOURCE_CODE
    if not (0.0 < value <= 1.0):
        raise ValueError(f"significance threshold must be in (0, 1]; got {value} from {source}")

    warnings: list[str] = []
    manifest_value: float | None = None
    if provenance_dir is not None:
        mpath = Path(provenance_dir) / "production_enrichment_manifest.json"
        if mpath.exists():
            try:
                doc = json.loads(mpath.read_text(encoding="utf-8"))
                raw = (doc.get("emitted_contract") or {}).get("fdr_significant_max")
                manifest_value = float(raw) if raw is not None else None
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                warnings.append(f"could not read {mpath.name}: {exc}")
            if manifest_value is not None and manifest_value != value:
                warnings.append(
                    f"threshold {value} ({source}) differs from the production manifest's "
                    f"emitted_contract.fdr_significant_max = {manifest_value}"
                )
    return ThresholdRecord(
        value=value,
        source=source,
        production_manifest_value=manifest_value,
        warnings=tuple(warnings),
    )


def resolve_secondary_fdrs(values: Sequence[float] | None) -> tuple[float, ...]:
    """Validated, de-duplicated, sorted secondary thresholds (default ``(0.005,)``)."""
    vals = tuple(float(v) for v in (values if values else DEFAULT_SECONDARY_FDRS))
    for v in vals:
        if not (0.0 < v <= 1.0):
            raise ValueError(f"secondary threshold must be in (0, 1]; got {v}")
    return tuple(sorted(set(vals)))


# --------------------------------------------------------------------------- #
# 2c. Run fingerprints (compatibility / reuse logic)
# --------------------------------------------------------------------------- #
def _per_tier_ints(per_tier: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, int]]:
    return {str(t): {str(k): int(v) for k, v in d.items()} for t, d in per_tier.items()}


def _cells_from_manifest_cells(
    cells: Iterable[Mapping[str, Any]], families: Iterable[str]
) -> tuple[tuple[str, str], ...]:
    fams = set(families)
    got = {(str(c.get("family")), str(c.get("tier"))) for c in cells if c.get("family") in fams}
    return tuple(sorted(got))


@dataclass(frozen=True)
class RunFingerprint:
    """The inputs that must match for one side's output to be REUSED.

    Every per-organism field is a ``{organism: ...}`` map so a fingerprint can be
    restricted to one organism (:meth:`subset`) and the reuse decision made per
    organism. Threshold-independent by construction (thresholds only affect the
    comparison, never the ORA outputs).

    ``source_sha256_recorded`` is False for a fingerprint whose source RECORDED
    no target-source digests (the production build manifests): the digest
    comparison is then skipped explicitly rather than passed vacuously, and the
    caller documents what covered the sources instead (the
    ``provenance/manifest.yaml`` sha256 gate + the ``per_tier_raw`` census).
    """

    organisms: tuple[str, ...]
    tiers: tuple[str, ...]
    families: tuple[str, ...]
    go_propagation: dict[str, str]
    org_db_versions: dict[str, str]
    go_db_versions: dict[str, str]
    r_versions: dict[str, str]
    clusterprofiler_versions: dict[str, str]
    source_sha256: dict[str, dict[str, str]]
    per_tier_raw: dict[str, dict[str, dict[str, int]]]
    n_universe: dict[str, dict[str, int]]
    n_terms_total: dict[str, dict[str, int]]
    normalization_mode: dict[str, str]
    cells: dict[str, tuple[tuple[str, str], ...]]
    bh_method: str = "BH"
    source_sha256_recorded: bool = True

    _PER_ORG = (
        "go_propagation",
        "org_db_versions",
        "go_db_versions",
        "r_versions",
        "clusterprofiler_versions",
        "source_sha256",
        "per_tier_raw",
        "n_universe",
        "n_terms_total",
        "normalization_mode",
        "cells",
    )

    def subset(self, organism: str) -> RunFingerprint:
        """Restrict to one organism (absent organism -> empty per-organism maps)."""
        kw: dict[str, Any] = {
            "organisms": (organism,) if organism in self.organisms else (),
            "tiers": self.tiers,
            "families": self.families,
            "bh_method": self.bh_method,
            "source_sha256_recorded": self.source_sha256_recorded,
        }
        for name in self._PER_ORG:
            src = getattr(self, name)
            kw[name] = {organism: src[organism]} if organism in src else {}
        return RunFingerprint(**kw)

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "organisms": list(self.organisms),
            "tiers": list(self.tiers),
            "families": list(self.families),
            "bh_method": self.bh_method,
            "source_sha256_recorded": self.source_sha256_recorded,
        }
        for name in self._PER_ORG:
            out[name] = getattr(self, name)
        out["cells"] = {org: [list(c) for c in cells] for org, cells in self.cells.items()}
        return out

    @classmethod
    def from_dict(cls, doc: Mapping[str, Any]) -> RunFingerprint:
        cells = {
            org: tuple(sorted((str(c[0]), str(c[1])) for c in lst))
            for org, lst in (doc.get("cells") or {}).items()
        }
        return cls(
            organisms=tuple(doc.get("organisms") or ()),
            tiers=tuple(doc.get("tiers") or ()),
            families=tuple(doc.get("families") or ()),
            go_propagation=dict(doc.get("go_propagation") or {}),
            org_db_versions=dict(doc.get("org_db_versions") or {}),
            go_db_versions=dict(doc.get("go_db_versions") or {}),
            r_versions=dict(doc.get("r_versions") or {}),
            clusterprofiler_versions=dict(doc.get("clusterprofiler_versions") or {}),
            source_sha256={k: dict(v) for k, v in (doc.get("source_sha256") or {}).items()},
            per_tier_raw={
                org: {t: dict(v) for t, v in tiers.items()}
                for org, tiers in (doc.get("per_tier_raw") or {}).items()
            },
            n_universe={k: dict(v) for k, v in (doc.get("n_universe") or {}).items()},
            n_terms_total={k: dict(v) for k, v in (doc.get("n_terms_total") or {}).items()},
            normalization_mode=dict(doc.get("normalization_mode") or {}),
            cells=cells,
            bh_method=str(doc.get("bh_method") or "BH"),
            source_sha256_recorded=bool(doc.get("source_sha256_recorded", True)),
        )

    def merged_with(self, other: RunFingerprint) -> RunFingerprint:
        """Union of two fingerprints (``other``'s organisms override on overlap)."""
        kw: dict[str, Any] = {
            "organisms": tuple(sorted(set(self.organisms) | set(other.organisms))),
            "tiers": tuple(sorted(set(self.tiers) | set(other.tiers))),
            "families": tuple(sorted(set(self.families) | set(other.families))),
            "bh_method": other.bh_method,
            "source_sha256_recorded": self.source_sha256_recorded
            and other.source_sha256_recorded,
        }
        for name in self._PER_ORG:
            kw[name] = {**getattr(self, name), **getattr(other, name)}
        return RunFingerprint(**kw)


def compare_fingerprints(expected: RunFingerprint, found: RunFingerprint) -> list[str]:
    """Human-readable mismatches between an expected and a found fingerprint.

    Empty list = compatible. Checks: every expected organism present; every
    per-organism scalar / map equal; every expected ``(family, tier)`` cell present
    (the found side may carry more cells); tiers / families of the expected run
    covered; BH method equal. The target-source digests are compared only when
    the found side RECORDED them (``source_sha256_recorded``); a side that did not
    (the production manifests) is skipped explicitly, never passed vacuously.
    """
    mism: list[str] = []
    if expected.bh_method != found.bh_method:
        mism.append(f"bh_method expected {expected.bh_method!r}, found {found.bh_method!r}")
    missing_t = sorted(set(expected.tiers) - set(found.tiers))
    if missing_t:
        mism.append(f"tiers missing from found run: {missing_t}")
    missing_f = sorted(set(expected.families) - set(found.families))
    if missing_f:
        mism.append(f"families missing from found run: {missing_f}")
    for org in expected.organisms:
        if org not in found.organisms:
            mism.append(f"{org}: absent from found manifest")
            continue
        for name in ("go_propagation", "normalization_mode"):
            e = getattr(expected, name).get(org)
            f = getattr(found, name).get(org)
            if e != f:
                mism.append(f"{org}: {name} expected {e!r}, found {f!r}")
        for name in ("org_db_versions", "go_db_versions", "r_versions", "clusterprofiler_versions"):
            e = getattr(expected, name).get(org)
            f = getattr(found, name).get(org)
            if e != f:
                mism.append(f"{org}: {name[:-1]} expected {e!r}, found {f!r}")
        e_sha = expected.source_sha256.get(org, {})
        f_sha = found.source_sha256.get(org, {})
        if not found.source_sha256_recorded:
            e_sha = {}
        for base in sorted(e_sha):
            if base not in f_sha:
                mism.append(f"{org}: source_sha256 for {base} missing in found run")
            elif e_sha[base] != f_sha[base]:
                mism.append(
                    f"{org}: source_sha256 for {base} expected {e_sha[base][:12]}.., "
                    f"found {f_sha[base][:12]}.."
                )
        e_ptr = expected.per_tier_raw.get(org, {})
        f_ptr = found.per_tier_raw.get(org, {})
        for tier in sorted(e_ptr):
            if e_ptr[tier] != f_ptr.get(tier):
                mism.append(
                    f"{org}: per_tier_raw[{tier}] expected {e_ptr[tier]}, found {f_ptr.get(tier)}"
                )
        for name in ("n_universe", "n_terms_total"):
            e_map = getattr(expected, name).get(org, {})
            f_map = getattr(found, name).get(org, {})
            for fam in sorted(e_map):
                if e_map[fam] != f_map.get(fam):
                    mism.append(
                        f"{org}: {name}[{fam}] expected {e_map[fam]}, found {f_map.get(fam)}"
                    )
        e_cells = set(expected.cells.get(org, ()))
        f_cells = set(found.cells.get(org, ()))
        missing_cells = sorted(e_cells - f_cells)
        if missing_cells:
            mism.append(f"{org}: cells missing from found run: {missing_cells}")
    return mism


def fingerprint_from_production_manifests(
    pathways_manifest: Mapping[str, Any] | None,
    enrichment_manifest: Mapping[str, Any] | None,
    *,
    source_sha: Mapping[str, Mapping[str, str]] | None = None,
    families: Sequence[str],
    tiers: Sequence[str],
) -> RunFingerprint:
    """Fingerprint of the PRODUCTION ``all`` outputs from the two build manifests.

    hsa/mmu come from ``production_pathways_manifest.json`` (organism blocks whose
    ``families[]`` carry GO subfamilies; ``cells[].cp_provenance`` carries the
    clusterProfiler version; target bands re-cased PER FAMILY). rno/mml come from
    ``production_manifest.json`` (``go_provenance`` / ``cp_provenance`` at the
    organism level; bands re-cased ONCE onto the UNION of the three GO universes).

    The production manifests record NO target-source digests (the builds ran
    behind the ``provenance/manifest.yaml`` sha256 gate instead). With
    ``source_sha=None`` (the driver's call) the fingerprint is marked
    ``source_sha256_recorded=False`` so :func:`compare_fingerprints` skips the
    digest check explicitly; a caller that independently verified digests
    (``{organism: {basename: sha256}}``) may supply them.
    """
    fams = tuple(families)
    go_prop: dict[str, str] = {}
    org_db: dict[str, str] = {}
    go_db: dict[str, str] = {}
    r_ver: dict[str, str] = {}
    cp_ver: dict[str, str] = {}
    per_tier: dict[str, dict[str, dict[str, int]]] = {}
    n_uni: dict[str, dict[str, int]] = {}
    n_terms: dict[str, dict[str, int]] = {}
    norm: dict[str, str] = {}
    cells: dict[str, tuple[tuple[str, str], ...]] = {}
    organisms: list[str] = []

    def _first(values: list[str], label: str, org: str) -> str:
        uniq = sorted(set(values))
        if len(uniq) > 1:
            raise ValueError(f"{org}: inconsistent {label} in production manifest: {uniq}")
        return uniq[0] if uniq else ""

    for block in (pathways_manifest or {}).get("organisms", []) or []:
        org = str(block.get("organism"))
        go_blocks = [b for b in block.get("families", []) or [] if b.get("family") in fams]
        if not go_blocks:
            continue
        provs = [b.get("provenance") or {} for b in go_blocks]

        def _prov(key: str, _provs: list = provs, _org: str = org) -> str:
            return _first([str(p.get(key, "")) for p in _provs], key, _org)

        go_prop[org] = _prov("go_propagation")
        org_db[org] = _prov("org_db_version")
        go_db[org] = _prov("go_db_version")
        r_ver[org] = _prov("r_version")
        go_cells = [c for c in block.get("cells", []) or [] if c.get("family") in fams]
        cp_ver[org] = _first(
            [
                str((c.get("cp_provenance") or {}).get("clusterprofiler_version", ""))
                for c in go_cells
            ],
            "clusterprofiler_version",
            org,
        )
        per_tier[org] = _per_tier_ints(block.get("per_tier_raw") or {})
        n_uni[org] = {str(b["family"]): int(b["n_universe"]) for b in go_blocks}
        n_terms[org] = {str(b["family"]): int(b["n_terms_total"]) for b in go_blocks}
        norm[org] = "per_family"
        cells[org] = _cells_from_manifest_cells(go_cells, fams)
        organisms.append(org)

    for block in (enrichment_manifest or {}).get("organisms", []) or []:
        org = str(block.get("organism"))
        if org in organisms:
            continue
        gp = block.get("go_provenance") or {}
        cpv = block.get("cp_provenance") or {}
        go_cells = [c for c in block.get("cells", []) or [] if c.get("family") in fams]
        if not go_cells and not gp:
            continue
        go_prop[org] = str(gp.get("go_propagation", ""))
        org_db[org] = str(gp.get("org_db_version", ""))
        go_db[org] = str(gp.get("go_db_version", ""))
        r_ver[org] = str(gp.get("r_version", ""))
        cp_ver[org] = str(cpv.get("clusterprofiler_version", ""))
        per_tier[org] = _per_tier_ints(block.get("per_tier_raw") or {})
        uni: dict[str, int] = {}
        tot: dict[str, int] = {}
        for c in go_cells:
            fam = str(c["family"])
            uni.setdefault(fam, int(c["n_universe"]))
            tot.setdefault(fam, int(c["n_go_terms_total"]))
        n_uni[org] = uni
        n_terms[org] = tot
        norm[org] = "union"
        cells[org] = _cells_from_manifest_cells(go_cells, fams)
        organisms.append(org)

    return RunFingerprint(
        organisms=tuple(sorted(organisms)),
        tiers=tuple(tiers),
        families=fams,
        go_propagation=go_prop,
        org_db_versions=org_db,
        go_db_versions=go_db,
        r_versions=r_ver,
        clusterprofiler_versions=cp_ver,
        source_sha256={org: dict((source_sha or {}).get(org, {})) for org in organisms},
        per_tier_raw=per_tier,
        n_universe=n_uni,
        n_terms_total=n_terms,
        normalization_mode=norm,
        cells=cells,
        source_sha256_recorded=source_sha is not None,
    )


def fingerprint_from_isolated_manifest(manifest_json: Mapping[str, Any]) -> RunFingerprint:
    """Fingerprint recorded in an isolated side manifest (``<out>/<mode>/manifest.json``)."""
    doc = manifest_json.get("fingerprint")
    if not isinstance(doc, Mapping):
        raise ValueError("isolated manifest has no 'fingerprint' block")
    return RunFingerprint.from_dict(doc)


# --------------------------------------------------------------------------- #
# 2d. Deterministic report writers
# --------------------------------------------------------------------------- #
SUMMARY_JSON = "go_propagation_summary.json"
STRUCTURAL_TSV = "go_structural_comparison.tsv"
STRUCTURAL_TERMS_TSV_GZ = "go_structural_terms.tsv.gz"
CELL_TSV = "go_cell_comparison.tsv"
GAINED_TSV = "go_terms_gained.tsv"
GAINED_FULL_TSV_GZ = "go_terms_gained_full.tsv.gz"
LOST_TSV = "go_terms_lost.tsv"
LOST_FULL_TSV_GZ = "go_terms_lost_full.tsv.gz"
REPORT_MD = "GO_PROPAGATION_SENSITIVITY.md"
MANIFEST_YAML = "manifest.yaml"

DETERMINISTIC_ARTEFACTS = (
    SUMMARY_JSON,
    STRUCTURAL_TSV,
    STRUCTURAL_TERMS_TSV_GZ,
    CELL_TSV,
    GAINED_TSV,
    GAINED_FULL_TSV_GZ,
    LOST_TSV,
    LOST_FULL_TSV_GZ,
    REPORT_MD,
)


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Mapping[str, str]]) -> int:
    """Write rows (dicts keyed by header) as TSV; gzip on ``.gz``. Returns row count."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with _open_text(path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(list(header))
        for row in rows:
            writer.writerow([row.get(h, "") for h in header])
            n += 1
    return n


def write_structural_tables(
    out_dir: str | Path,
    stats: Sequence[StructuralStats],
    *,
    term_names: Mapping[str, Mapping[str, str]] | None = None,
) -> dict[str, int]:
    """Write ``go_structural_comparison.tsv`` + ``go_structural_terms.tsv.gz``."""
    out = Path(out_dir)
    ordered = sorted(stats, key=lambda s: (s.organism, s.family))
    rows = [s.to_row() for s in ordered]
    n_summary = _write_tsv(out / STRUCTURAL_TSV, STRUCTURAL_TSV_HEADER, rows)

    def _term_rows():
        for s in ordered:
            names = (term_names or {}).get(s.organism, {})
            for t in s.terms:
                yield {
                    "organism": s.organism,
                    "family": s.family,
                    "term": t.term,
                    "name": names.get(t.term, ""),
                    "size_direct": _fmt(t.size_direct),
                    "size_all": _fmt(t.size_all),
                    "size_delta": _fmt(t.size_delta),
                    "depth": _fmt(t.depth),
                    "status": t.status,
                }

    n_terms = _write_tsv(out / STRUCTURAL_TERMS_TSV_GZ, STRUCTURAL_TERMS_HEADER, _term_rows())
    return {"structural_rows": n_summary, "structural_term_rows": n_terms}


def _cell_header(comparisons: Sequence[CellComparison]) -> list[str]:
    if comparisons:
        return list(comparisons[0].to_row())
    return []


def write_cell_tables(out_dir: str | Path, comparisons: Sequence[CellComparison]) -> int:
    """Write ``go_cell_comparison.tsv`` (one row per organism x family x tier)."""
    ordered = sorted(comparisons, key=lambda c: (c.organism, c.family, c.tier))
    rows = [c.to_row() for c in ordered]
    header = _cell_header(ordered)
    return _write_tsv(Path(out_dir) / CELL_TSV, header, rows)


def write_term_tables(
    out_dir: str | Path, records: Sequence[TermRecord], *, top_n: int = TOP_N_RANKED
) -> dict[str, int]:
    """Write gained / lost tables (top ``top_n`` + gzipped full versions)."""
    out = Path(out_dir)
    gained = rank_gained(records)
    lost = rank_lost(records)
    return {
        "gained_top": _write_tsv(
            out / GAINED_TSV, TERM_RECORD_HEADER, (r.to_row() for r in gained[:top_n])
        ),
        "gained_full": _write_tsv(
            out / GAINED_FULL_TSV_GZ, TERM_RECORD_HEADER, (r.to_row() for r in gained)
        ),
        "lost_top": _write_tsv(
            out / LOST_TSV, TERM_RECORD_HEADER, (r.to_row() for r in lost[:top_n])
        ),
        "lost_full": _write_tsv(
            out / LOST_FULL_TSV_GZ, TERM_RECORD_HEADER, (r.to_row() for r in lost)
        ),
    }


def write_summary_json(out_dir: str | Path, summary: Mapping[str, Any]) -> Path:
    """Write ``go_propagation_summary.json`` (sorted keys, indent 2, trailing newline)."""
    path = Path(out_dir) / SUMMARY_JSON
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _md_table(header: Sequence[str], rows: Iterable[Sequence[object]]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(v) for v in row) + " |")
    return lines


def _r4(v: object) -> object:
    return round(v, 4) if isinstance(v, float) else v


def render_report_md(summary: Mapping[str, Any]) -> str:
    """Render the human report from the summary dict (deterministic, no timestamps)."""
    thr = summary["thresholds"]
    sides = summary["sides"]
    lines: list[str] = []
    lines.append("# GO propagation sensitivity: `all` (production) vs `direct`")
    lines.append("")
    lines.append(
        "Production GO gene sets use the propagated `GOALL`/`ONTOLOGYALL` org.*.eg.db "
        "columns (`config.yaml` `go_propagation: all`, unchanged). This report quantifies "
        "how the structure of the gene sets and the enrichment calls would differ under "
        "the `direct` (`GO`/`ONTOLOGY`) annotation columns. Nothing under "
        "`data/processed/` was modified; the `direct` side was computed in isolation."
    )
    lines.append("")
    lines.append("## Definitions")
    lines.append("")
    for key, text in sorted(summary["definitions"].items()):
        lines.append(f"* **{key}**: {text}")
    lines.append("")
    lines.append("## Thresholds")
    lines.append("")
    prim = thr["primary"]
    lines.append(
        f"* primary (production SIGNIFICANT criterion): `p_adjust <= {prim['value']}` "
        f"(source: `{prim['source']}`; production manifest value: "
        f"{prim['production_manifest_value']})"
    )
    for w in prim.get("warnings", []):
        lines.append(f"  * WARNING: {w}")
    lines.append(f"* secondary: {thr['secondary']['values']} -- {thr['secondary']['label']}")
    lines.append(f"* multiple testing: {summary['multiple_testing']}")
    lines.append(f"* universe: {summary['universe_definition']}")
    lines.append("")
    lines.append("## Reuse decisions")
    lines.append("")
    for side in ("all", "direct"):
        blk = sides.get(side) or {}
        lines.append(f"### `{side}` side")
        lines.append("")
        lines.append(f"* source: {blk.get('source')}")
        if blk.get("manifest_path"):
            lines.append(
                f"* manifest: `{blk.get('manifest_path')}` (sha256 recorded in manifest.yaml)"
            )
        for org, dec in sorted((blk.get("decisions") or {}).items()):
            lines.append(f"* {org}: {dec.get('source')}")
            for m in dec.get("mismatches_production") or []:
                lines.append(f"  * production output not reusable: {m}")
            for note in dec.get("notes") or []:
                lines.append(f"  * note: {note}")
        prod_cells = blk.get("production_cells") or {}
        if prod_cells:
            lines.append("")
            lines.append(
                "Reused production cells (authoritative rows overlaid on the recomputed "
                "BH; `padj_floor` = 0.25 means cells absent from those rows are known to "
                "have p_adjust > 0.25 -- except for the `floor-exempt` miRNAs, whose "
                "non-overlaid cells use the recomputed BH as is; an empty `padj_floor` "
                "means the recomputed BH was used as is for every cell):"
            )
            lines.append("")
            pc_rows = []
            for org in sorted(prod_cells):
                for key in sorted(prod_cells[org]):
                    pc = prod_cells[org][key]
                    pc_rows.append(
                        [
                            org,
                            key,
                            pc.get("n_padj_authoritative"),
                            pc.get("n_significant_rows_unmatched"),
                            pc.get("species_excluded_tsv") is not None,
                            pc.get("n_floor_exempt_mirnas"),
                            pc.get("n_padj_floor_exempt"),
                            pc.get("significant_tsv_sha256_verified"),
                            pc.get("padj_floor"),
                            pc.get("padj_floor_reason"),
                        ]
                    )
            lines.extend(
                _md_table(
                    [
                        "organism",
                        "family/tier",
                        "authoritative cells",
                        "unmatched rows",
                        "quarantine overlaid",
                        "floor-exempt miRNAs",
                        "floor-exempt cells",
                        "TSV sha256 = contract",
                        "padj_floor",
                        "reason",
                    ],
                    pc_rows,
                )
            )
        lines.append("")
    lines.append("## Structural comparison (per organism x family)")
    lines.append("")
    struct_rows = []
    for org in sorted(summary["structural"]):
        for fam in sorted(summary["structural"][org]):
            s = summary["structural"][org][fam]
            struct_rows.append(
                [
                    org,
                    fam,
                    s["n_terms_direct"],
                    s["n_terms_all"],
                    s["n_assoc_direct"],
                    s["n_assoc_all"],
                    s["n_assoc_propagated"],
                    _r4(s["pct_assoc_propagated"]),
                    s["n_direct_only_assoc"],
                    s["term_size_median_direct"],
                    s["term_size_median_all"],
                    s["n_terms_only_in_all"],
                    s["n_terms_both_changed_membership"],
                    s["n_terms_both_identical"],
                    s["n_terms_depth_unreachable"],
                ]
            )
    lines.extend(
        _md_table(
            [
                "organism",
                "family",
                "terms direct",
                "terms all",
                "assoc direct",
                "assoc all",
                "propagated",
                "% propagated",
                "direct-only assoc (anomaly)",
                "median size direct",
                "median size all",
                "terms only in all",
                "changed membership",
                "identical",
                "depth unreachable",
            ],
            struct_rows,
        )
    )
    lines.append("")
    lines.append("## Enrichment comparison (per organism x family x tier)")
    lines.append("")
    if summary.get("cells"):
        cell_rows = []
        for c in sorted(summary["cells"], key=lambda c: (c["organism"], c["family"], c["tier"])):
            cell_rows.append(
                [
                    c["organism"],
                    c["family"],
                    c["tier"],
                    c["n_tested_direct"],
                    c["n_tested_all"],
                    c["n_shared_tested"],
                    c["n_sig_direct"],
                    c["n_sig_all"],
                    c["n_sig_intersection"],
                    _r4(c["jaccard"]),
                    c["n_sig_only_all"],
                    c["n_sig_only_direct"],
                    _r4(c["spearman_p"]),
                    _r4(c["spearman_padj"]),
                ]
            )
        lines.extend(
            _md_table(
                [
                    "organism",
                    "family",
                    "tier",
                    "tested direct",
                    "tested all",
                    "shared",
                    "sig direct",
                    "sig all",
                    "sig both",
                    "Jaccard",
                    "gained (all only)",
                    "lost (direct only)",
                    "Spearman -log10 p",
                    "Spearman -log10 padj",
                ],
                cell_rows,
            )
        )
        lines.append("")
        lines.append("### Affected cells by organism x tier (totals over families)")
        lines.append("")
        aff_rows = []
        for org in sorted(summary["affected_summary"]):
            for tier in sorted(summary["affected_summary"][org]):
                a = summary["affected_summary"][org][tier]
                aff_rows.append(
                    [org, tier, a["n_sig_direct"], a["n_sig_all"], a["gained"], a["lost"]]
                )
        lines.extend(
            _md_table(["organism", "tier", "sig direct", "sig all", "gained", "lost"], aff_rows)
        )
        lines.append("")
        lines.append("## Ranked terms (formulas below; no manual selection)")
        lines.append("")
        for key, title in (
            ("top_gained", "Top gained terms"),
            ("top_lost", "Top lost terms"),
            ("largest_median_shift", "Largest median -log10(padj) shift"),
            ("largest_membership_expansion", "Largest membership expansion"),
            ("broad_parent_terms", "Broad parent terms (depth <= 2) by gained"),
        ):
            recs = summary["rankings"].get(key) or []
            lines.append(f"### {title}")
            lines.append("")
            if not recs:
                lines.append("(none)")
            else:
                lines.extend(
                    _md_table(
                        [
                            "organism",
                            "family",
                            "term",
                            "name",
                            "sig direct",
                            "sig all",
                            "gained",
                            "lost",
                            "median shift",
                            "size direct",
                            "size all",
                            "depth",
                        ],
                        (
                            [
                                r["organism"],
                                r["family"],
                                r["term"],
                                r["name"],
                                r["n_sig_direct"],
                                r["n_sig_all"],
                                r["gained"],
                                r["lost"],
                                _r4(r["median_neglog10_padj_shift"]),
                                r["size_direct"],
                                r["size_all"],
                                r["depth"],
                            ]
                            for r in recs[:25]
                        ),
                    )
                )
            lines.append("")
    else:
        lines.append("(structural-only run: no enrichment comparison)")
        lines.append("")
    lines.append("## Ranking formulas")
    lines.append("")
    for key, text in sorted(summary["ranking_formulas"].items()):
        lines.append(f"* `{key}`: {text}")
    lines.append("")
    lines.append("## Caveats")
    lines.append("")
    for c in summary["caveats"]:
        lines.append(f"* {c}")
    lines.append("")
    lines.append("## Artefacts")
    lines.append("")
    for name in DETERMINISTIC_ARTEFACTS:
        lines.append(f"* `{name}`")
    lines.append(f"* `{MANIFEST_YAML}` (the only file carrying timestamps / wall times)")
    lines.append("")
    return "\n".join(lines)


def write_report_md(out_dir: str | Path, summary: Mapping[str, Any]) -> Path:
    path = Path(out_dir) / REPORT_MD
    path.write_text(render_report_md(summary), encoding="utf-8")
    return path


CAVEATS: tuple[str, ...] = (
    "The reused production `all` matrices store only the raw p and write 1.0 for untested "
    "cells; 'tested' is therefore the proxy p < 1.0 and p_adjust is recomputed by BH over "
    "that proxy, except for cells in the authoritative significant rows (p_adjust <= 0.25), "
    "which carry the build's own p_adjust + overlap. Those rows are read from BOTH the "
    "shipped enrichment_significant TSV and the species_excluded quarantine written by the "
    "contract export (foreign-species miRNAs are removed from the shipped TSV only, never "
    "from the matrices, and the isolated direct side scores them too). A cell absent from "
    "the overlaid rows is treated as p_adjust > 0.25 (padj_floor = 0.25) only for miRNAs "
    "those rows are known to cover -- every miRNA when the quarantine is overlaid or the "
    "TSV was never filtered, else only the miRNAs the export-time species policy retains "
    "(foreign-species / unparseable miRNAs are floor-exempt and use the recomputed BH as "
    "is; see the reused-production-cells table). Where the floor applies, significance "
    "calls at thresholds <= 0.25 are exact; calls at thresholds > 0.25 on a reused "
    "production side always use the recomputed BH (approximate: the p < 1.0 proxy shrinks "
    "BH's n). padj correlations / deltas / per-term shifts always use the raw recomputed "
    "values.",
    "Terms whose upper-tail p is exactly 1.0 (e.g. very broad propagated terms containing "
    "nearly the whole universe) are invisible under the production tested proxy; they appear "
    "as 'only direct tested' cells when the direct side is exact.",
    "Both sides use the same target sets (production loaders), but each side re-cases the "
    "bands onto its own GO symbol space exactly as the runners do (hsa/mmu per family; "
    "rno/mml onto the union of the three GO family universes).",
    "The GO depth is a graph distance over all GO.db parent relationship types; it is NOT the "
    "annotation-count 'level' used by some tools.",
    "Secondary thresholds are sensitivity checks only; the production criterion is the "
    "primary threshold recorded above.",
    "The production build manifests record no target-source sha256; a reused production "
    "`all` side is therefore accepted on the provenance/manifest.yaml sha256 gate plus the "
    "per_tier_raw target-set census, universes, term counts and R/package versions -- not on "
    "per-run source digests (see manifest.yaml provenance_gate and the side notes).",
)


__all__ = [
    "CAVEATS",
    "CellComparison",
    "CellTable",
    "DEFAULT_SECONDARY_FDRS",
    "DETERMINISTIC_ARTEFACTS",
    "GO_DEPTH_DEFINITION",
    "GO_ROOTS",
    "MULTIPLE_TESTING",
    "NEG_LOG10_CAP",
    "PADJ_AUTHORITATIVE",
    "PADJ_RECOMPUTED",
    "RANKING_FORMULAS",
    "RunFingerprint",
    "SECONDARY_FDR_LABEL",
    "StructuralStats",
    "TermRecord",
    "TermSize",
    "ThresholdRecord",
    "UNIVERSE_DEFINITION",
    "affected_summary",
    "aggregate_terms",
    "bh_adjust",
    "broad_parent_terms",
    "compare_cell_tables",
    "compare_fingerprints",
    "compute_depths",
    "depths_by_family",
    "fingerprint_from_isolated_manifest",
    "fingerprint_from_production_manifests",
    "load_cell_table_from_long_form",
    "load_cell_table_from_production",
    "mirna_prefix",
    "mirna_retained_by_species_policy",
    "neg_log10",
    "production_padj_floor",
    "quantiles",
    "rank_gained",
    "rank_lost",
    "read_go_parent_edges",
    "render_report_md",
    "resolve_secondary_fdrs",
    "resolve_significance_threshold",
    "structural_comparison",
    "top_by_abs_shift",
    "top_by_size_delta",
    "write_cell_tables",
    "write_cells_long_form",
    "write_matrix",
    "write_report_md",
    "write_structural_tables",
    "write_summary_json",
    "write_term_tables",
]
