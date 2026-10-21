"""Cell-for-cell validation of the controlled v2 reproduction.

This is the **engine-isolating** comparison: our clusterProfiler matrices vs the
frozen miRPathDB2 GeneTrail2 p-value matrices, on identical (miRNA, pathway)
frames built from v2's own gene sets. Because both sides share v2's exact
gene-set membership (we reconstructed TERM2GENE from the GeneTrail2 ``Info``
column) and v2's exact miRNA-row / pathway-column order, any divergence is
attributable to the **enrichment engine** (clusterProfiler's hypergeometric ORA
vs GeneTrail2's), NOT to gene-set drift or predictor-version differences.

The validator emits, per comparison cell (organism x variant x family):

* **Exact concordance** of the comparison frame — pathway set and miRNA set must
  be byte-identical to v2 (HARD FAIL otherwise), plus the matched-miRNA /
  matched-pathway counts and any rows/cols that could not be matched.
* **Tested-coverage** asymmetry — how many cells each side actually scored
  (``p < 1.0``; v2's sentinel ``1.0`` == not tested), the both-tested overlap,
  and the new-only / v2-only tested counts.
* **Rank correlation** — Spearman rho on ``-log10(p)`` over the both-tested
  cells.
* **Directional significance** — confusion of significant / not-significant
  calls at raw ``p <= 0.05`` AND BH-FDR ``<= 0.05`` (both reported); agreement
  %, both-sig / v2-only / new-only / both-not counts.
* **Numeric divergence** — the distribution of ``Delta = (-log10 p_new) -
  (-log10 p_v2)`` over both-tested cells: mean, median, IQR, 5th/95th pct, and
  worst-case.

Thresholds (explicit): significance is evaluated at raw ``p <= 0.05`` and,
separately, at ``BH-FDR <= 0.05``. **BH context** — v2's matrices store raw
ORA p-values, so for the FDR axis we apply Benjamini-Hochberg *within each
matrix cell* over that cell's tested p-values (``p < 1.0``), identically to both
sides, so the comparison is apples-to-apples. This is a comparison-time BH for
the directional-agreement axis only; it is NOT the production per-(miRNA, tier,
db) BH family. The not-tested sentinel (``1.0``) is excluded from BH so it does
not inflate the denominator.

All functions here are pure (NumPy / SciPy / statsmodels); figure generation and
file IO live in the driver script.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# v2's "not tested" sentinel (a literal 1.0 in the matrix).
NOT_TESTED = 1.0
RAW_ALPHA = 0.05
FDR_ALPHA = 0.05
# -log10 floor for p == 0 (avoid +inf); v2 stores doubles, true 0 is rare but
# can occur from underflow. Mirror it on both sides so it cancels in Delta.
_NEG_LOG10_CAP = 350.0


def _neg_log10(p: np.ndarray) -> np.ndarray:
    """``-log10(p)`` with a finite cap for ``p == 0`` (applied to both sides)."""
    with np.errstate(divide="ignore"):
        out = -np.log10(p)
    return np.minimum(out, _NEG_LOG10_CAP)


@dataclass
class FrameConcordance:
    """Exact-concordance result for one cell's comparison frame.

    Attributes:
        organism / variant / family: Cell identity.
        n_pathways_new / n_pathways_v2: Pathway (column) counts.
        n_mirnas_new / n_mirnas_v2: miRNA (row) counts.
        pathways_identical / mirnas_identical: Whether the ordered label lists
            are byte-identical (the HARD concordance requirement).
        matched_pathways / matched_mirnas: Count of labels present in BOTH (in
            the identical case, == the full count).
        v2_only_pathways / new_only_pathways: Labels present on only one side.
        v2_only_mirnas / new_only_mirnas: Same, for miRNAs.
        passed: True iff both label lists are byte-identical.
        reasons: Human-readable explanations for any unmatched labels / failure.
    """

    organism: str
    variant: str
    family: str
    n_pathways_new: int
    n_pathways_v2: int
    n_mirnas_new: int
    n_mirnas_v2: int
    pathways_identical: bool
    mirnas_identical: bool
    matched_pathways: int
    matched_mirnas: int
    v2_only_pathways: list[str]
    new_only_pathways: list[str]
    v2_only_mirnas: list[str]
    new_only_mirnas: list[str]
    passed: bool
    reasons: list[str] = field(default_factory=list)


def check_frame_concordance(
    *,
    organism: str,
    variant: str,
    family: str,
    new_pathways: Sequence[str],
    new_mirnas: Sequence[str],
    v2_pathways: Sequence[str],
    v2_mirnas: Sequence[str],
) -> FrameConcordance:
    """Verify our cell's frame is byte-identical to v2's (exact concordance).

    Because we reuse v2's GMT + v2's matrix row/col order by construction, the
    label lists MUST match exactly; this is the gate-0 concordance check. Any
    drift is a hard fail and is reported with the offending labels.
    """
    new_p, v2_p = list(new_pathways), list(v2_pathways)
    new_m, v2_m = list(new_mirnas), list(v2_mirnas)

    pathways_identical = new_p == v2_p
    mirnas_identical = new_m == v2_m

    new_pset, v2_pset = set(new_p), set(v2_p)
    new_mset, v2_mset = set(new_m), set(v2_m)
    v2_only_pathways = sorted(v2_pset - new_pset)
    new_only_pathways = sorted(new_pset - v2_pset)
    v2_only_mirnas = sorted(v2_mset - new_mset)
    new_only_mirnas = sorted(new_mset - v2_mset)

    reasons: list[str] = []
    if not pathways_identical:
        if v2_only_pathways or new_only_pathways:
            reasons.append(
                f"pathway label sets differ: {len(v2_only_pathways)} v2-only, "
                f"{len(new_only_pathways)} new-only"
            )
        else:
            reasons.append("pathway labels identical as sets but ORDER differs")
    if not mirnas_identical:
        if v2_only_mirnas or new_only_mirnas:
            reasons.append(
                f"miRNA label sets differ: {len(v2_only_mirnas)} v2-only, "
                f"{len(new_only_mirnas)} new-only"
            )
        else:
            reasons.append("miRNA labels identical as sets but ORDER differs")

    passed = pathways_identical and mirnas_identical
    return FrameConcordance(
        organism=organism,
        variant=variant,
        family=family,
        n_pathways_new=len(new_p),
        n_pathways_v2=len(v2_p),
        n_mirnas_new=len(new_m),
        n_mirnas_v2=len(v2_m),
        pathways_identical=pathways_identical,
        mirnas_identical=mirnas_identical,
        matched_pathways=len(new_pset & v2_pset),
        matched_mirnas=len(new_mset & v2_mset),
        v2_only_pathways=v2_only_pathways,
        new_only_pathways=new_only_pathways,
        v2_only_mirnas=v2_only_mirnas,
        new_only_mirnas=new_only_mirnas,
        passed=passed,
        reasons=reasons,
    )


def _bh_significant(
    pmatrix: np.ndarray, tested: np.ndarray, alpha: float
) -> np.ndarray:
    """BH-significant boolean mask over a p-value matrix.

    BH is applied over the *tested* cells only (``p < 1.0``); not-tested cells
    (the ``1.0`` sentinel) are excluded so they cannot inflate the BH
    denominator. Returns a full-shape boolean mask (False at not-tested cells).
    """
    sig = np.zeros(pmatrix.shape, dtype=bool)
    flat_tested = tested.ravel()
    if not flat_tested.any():
        return sig
    pflat = pmatrix.ravel()[flat_tested]
    reject, _, _, _ = multipletests(pflat, alpha=alpha, method="fdr_bh")
    sig_flat = np.zeros(pmatrix.size, dtype=bool)
    sig_flat[np.flatnonzero(flat_tested)] = reject
    return sig_flat.reshape(pmatrix.shape)


@dataclass
class Confusion:
    """A significant/not 2x2 confusion + agreement for one threshold axis.

    Counts are over the comparison domain (cells tested by at least one side,
    unless ``domain == 'both'``). ``both_sig`` = sig in new AND v2; ``v2_only`` =
    sig in v2 not new; ``new_only`` = sig in new not v2; ``both_not`` = sig in
    neither.
    """

    threshold: str  # "raw_p<=0.05" | "bh_fdr<=0.05"
    domain: str  # "either_tested" | "both_tested"
    n: int
    both_sig: int
    v2_only: int
    new_only: int
    both_not: int
    agreement: float  # fraction (both_sig + both_not) / n


def _confusion(
    sig_new: np.ndarray,
    sig_v2: np.ndarray,
    mask: np.ndarray,
    threshold: str,
    domain: str,
) -> Confusion:
    """Build a confusion matrix from two significance masks over ``mask`` cells."""
    sn = sig_new[mask]
    sv = sig_v2[mask]
    both_sig = int(np.sum(sn & sv))
    v2_only = int(np.sum(~sn & sv))
    new_only = int(np.sum(sn & ~sv))
    both_not = int(np.sum(~sn & ~sv))
    n = int(mask.sum())
    agree = float((both_sig + both_not) / n) if n else float("nan")
    return Confusion(
        threshold=threshold,
        domain=domain,
        n=n,
        both_sig=both_sig,
        v2_only=v2_only,
        new_only=new_only,
        both_not=both_not,
        agreement=agree,
    )


@dataclass
class CellComparison:
    """Full per-cell comparison metrics (the validator's atomic unit)."""

    organism: str
    variant: str
    family: str
    # frame
    concordance: FrameConcordance
    # coverage
    n_tested_new: int
    n_tested_v2: int
    n_both_tested: int
    n_new_only_tested: int
    n_v2_only_tested: int
    n_either_tested: int
    # rank correlation (on both-tested -log10 p)
    spearman_rho: float
    spearman_p: float
    # directional agreement
    confusions: list[Confusion]
    # numeric divergence (Delta = -log10 p_new - -log10 p_v2 over both-tested)
    delta_mean: float
    delta_median: float
    delta_iqr: float
    delta_p05: float
    delta_p95: float
    delta_min: float
    delta_max: float
    # paired arrays for figures (both-tested -log10 p)
    nlp_new: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))
    nlp_v2: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))

    def metrics_dict(self) -> dict:
        """JSON-serializable metrics view (drops the large paired arrays)."""
        c = self.concordance
        return {
            "organism": self.organism,
            "variant": self.variant,
            "family": self.family,
            "concordance": {
                "passed": c.passed,
                "pathways_identical": c.pathways_identical,
                "mirnas_identical": c.mirnas_identical,
                "n_pathways_new": c.n_pathways_new,
                "n_pathways_v2": c.n_pathways_v2,
                "n_mirnas_new": c.n_mirnas_new,
                "n_mirnas_v2": c.n_mirnas_v2,
                "matched_pathways": c.matched_pathways,
                "matched_mirnas": c.matched_mirnas,
                "v2_only_pathways": c.v2_only_pathways,
                "new_only_pathways": c.new_only_pathways,
                "v2_only_mirnas": c.v2_only_mirnas,
                "new_only_mirnas": c.new_only_mirnas,
                "reasons": c.reasons,
            },
            "coverage": {
                "n_tested_new": self.n_tested_new,
                "n_tested_v2": self.n_tested_v2,
                "n_both_tested": self.n_both_tested,
                "n_new_only_tested": self.n_new_only_tested,
                "n_v2_only_tested": self.n_v2_only_tested,
                "n_either_tested": self.n_either_tested,
            },
            "spearman_rho": self.spearman_rho,
            "spearman_p": self.spearman_p,
            "confusions": [
                {
                    "threshold": cf.threshold,
                    "domain": cf.domain,
                    "n": cf.n,
                    "both_sig": cf.both_sig,
                    "v2_only": cf.v2_only,
                    "new_only": cf.new_only,
                    "both_not": cf.both_not,
                    "agreement": cf.agreement,
                }
                for cf in self.confusions
            ],
            "delta": {
                "mean": self.delta_mean,
                "median": self.delta_median,
                "iqr": self.delta_iqr,
                "p05": self.delta_p05,
                "p95": self.delta_p95,
                "min": self.delta_min,
                "max": self.delta_max,
            },
        }


def compare_cell(
    *,
    organism: str,
    variant: str,
    family: str,
    new_pathways: Sequence[str],
    new_mirnas: Sequence[str],
    new_values: np.ndarray,
    v2_pathways: Sequence[str],
    v2_mirnas: Sequence[str],
    v2_values: np.ndarray,
) -> CellComparison:
    """Compute the full comparison for one (organism, variant, family) cell.

    Args:
        new_*: Our matrix labels + ``(n_mirnas, n_pathways)`` p-value array.
        v2_*: The frozen v2 matrix labels + p-value array.

    Both matrices MUST share the same orientation (miRNA rows, pathway cols) and,
    by construction, the same label order — the frame concordance check verifies
    this and the metrics below assume the arrays are aligned.

    Returns:
        A :class:`CellComparison` with concordance, coverage, Spearman rho,
        raw + BH-FDR directional confusions, and Delta distribution stats.
    """
    concordance = check_frame_concordance(
        organism=organism,
        variant=variant,
        family=family,
        new_pathways=new_pathways,
        new_mirnas=new_mirnas,
        v2_pathways=v2_pathways,
        v2_mirnas=v2_mirnas,
    )

    new = np.asarray(new_values, dtype=np.float64)
    v2 = np.asarray(v2_values, dtype=np.float64)
    if new.shape != v2.shape:
        raise ValueError(
            f"{organism}/{family}/{variant}: shape mismatch new {new.shape} vs v2 {v2.shape}"
        )

    tested_new = new < NOT_TESTED
    tested_v2 = v2 < NOT_TESTED
    both = tested_new & tested_v2
    either = tested_new | tested_v2

    n_both = int(both.sum())

    # --- rank correlation on both-tested -log10 p ---
    nlp_new_both = _neg_log10(new[both])
    nlp_v2_both = _neg_log10(v2[both])
    if n_both >= 2 and np.ptp(nlp_new_both) > 0 and np.ptp(nlp_v2_both) > 0:
        rho, sp = spearmanr(nlp_new_both, nlp_v2_both)
        rho, sp = float(rho), float(sp)
    else:
        rho, sp = float("nan"), float("nan")

    # --- directional significance (raw + BH-FDR) ---
    sig_new_raw = tested_new & (new <= RAW_ALPHA)
    sig_v2_raw = tested_v2 & (v2 <= RAW_ALPHA)
    sig_new_bh = _bh_significant(new, tested_new, FDR_ALPHA)
    sig_v2_bh = _bh_significant(v2, tested_v2, FDR_ALPHA)

    confusions = [
        _confusion(sig_new_raw, sig_v2_raw, either, "raw_p<=0.05", "either_tested"),
        _confusion(sig_new_bh, sig_v2_bh, either, "bh_fdr<=0.05", "either_tested"),
    ]

    # --- numeric divergence (Delta over both-tested) ---
    if n_both:
        delta = nlp_new_both - nlp_v2_both
        q25, q75 = np.percentile(delta, [25, 75])
        d = {
            "mean": float(np.mean(delta)),
            "median": float(np.median(delta)),
            "iqr": float(q75 - q25),
            "p05": float(np.percentile(delta, 5)),
            "p95": float(np.percentile(delta, 95)),
            "min": float(np.min(delta)),
            "max": float(np.max(delta)),
        }
    else:
        d = dict.fromkeys(("mean", "median", "iqr", "p05", "p95", "min", "max"), float("nan"))

    return CellComparison(
        organism=organism,
        variant=variant,
        family=family,
        concordance=concordance,
        n_tested_new=int(tested_new.sum()),
        n_tested_v2=int(tested_v2.sum()),
        n_both_tested=n_both,
        n_new_only_tested=int((tested_new & ~tested_v2).sum()),
        n_v2_only_tested=int((tested_v2 & ~tested_new).sum()),
        n_either_tested=int(either.sum()),
        spearman_rho=rho,
        spearman_p=sp,
        confusions=confusions,
        delta_mean=d["mean"],
        delta_median=d["median"],
        delta_iqr=d["iqr"],
        delta_p05=d["p05"],
        delta_p95=d["p95"],
        delta_min=d["min"],
        delta_max=d["max"],
        nlp_new=nlp_new_both,
        nlp_v2=nlp_v2_both,
    )


def pooled_spearman(comparisons: Sequence[CellComparison]) -> tuple[float, float, int]:
    """Spearman rho on the pooled both-tested ``-log10 p`` across cells.

    Returns ``(rho, p_value, n_pairs)``. Used for the per-organism and
    per-family roll-ups (a single rho over all matched cells in the group).
    """
    a = np.concatenate([c.nlp_new for c in comparisons]) if comparisons else np.empty(0)
    b = np.concatenate([c.nlp_v2 for c in comparisons]) if comparisons else np.empty(0)
    if a.size >= 2 and np.ptp(a) > 0 and np.ptp(b) > 0:
        rho, sp = spearmanr(a, b)
        return float(rho), float(sp), int(a.size)
    return float("nan"), float("nan"), int(a.size)


def pooled_confusion(
    comparisons: Sequence[CellComparison], threshold: str
) -> Confusion:
    """Sum the confusion counts across cells for one threshold axis."""
    both_sig = v2_only = new_only = both_not = 0
    domain = "either_tested"
    for c in comparisons:
        for cf in c.confusions:
            if cf.threshold == threshold:
                both_sig += cf.both_sig
                v2_only += cf.v2_only
                new_only += cf.new_only
                both_not += cf.both_not
                domain = cf.domain
    n = both_sig + v2_only + new_only + both_not
    agree = float((both_sig + both_not) / n) if n else float("nan")
    return Confusion(
        threshold=threshold,
        domain=domain,
        n=n,
        both_sig=both_sig,
        v2_only=v2_only,
        new_only=new_only,
        both_not=both_not,
        agreement=agree,
    )


def pooled_delta(comparisons: Sequence[CellComparison]) -> dict:
    """Pooled Delta distribution stats over all both-tested cells in a group."""
    deltas = [c.nlp_new - c.nlp_v2 for c in comparisons if c.nlp_new.size]
    if not deltas:
        keys = ("mean", "median", "iqr", "p05", "p95", "min", "max", "n")
        return dict.fromkeys(keys, float("nan"))
    delta = np.concatenate(deltas)
    q25, q75 = np.percentile(delta, [25, 75])
    return {
        "mean": float(np.mean(delta)),
        "median": float(np.median(delta)),
        "iqr": float(q75 - q25),
        "p05": float(np.percentile(delta, 5)),
        "p95": float(np.percentile(delta, 95)),
        "min": float(np.min(delta)),
        "max": float(np.max(delta)),
        "n": int(delta.size),
    }
