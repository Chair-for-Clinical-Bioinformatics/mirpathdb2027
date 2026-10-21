"""Unchanged production kernel adapter + agreement audit against :mod:`.kernel`.

:func:`run_production_kernel` calls
:func:`mirpathdb_prod.enrichment.clusterprofiler.run_ora_batch_with_provenance`
-- the exact function and parameters every production runner used -- on a
:class:`FamilySpace` and a set of queries. The memberships handed to R are the
space's eligible memberships and the universe is the space's background, so a
restricted background (target-union) is passed to R as restricted memberships
too (``enricher`` counts ``n`` over annotated query genes, which coincides with
"query genes inside the background" only when every annotated gene is in the
background).

:func:`agreement_report` compares the R output with the independent kernel on
the SAME queries: scored term sets must coincide, overlaps must match exactly,
raw and adjusted P-values must agree to a relative tolerance.
"""

from __future__ import annotations

import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from mirpathdb_prod.enrichment import clusterprofiler as cp
from mirpathdb_prod.validation.robustness.inputs import FamilySpace
from mirpathdb_prod.validation.robustness.kernel import OraBatch


class ProductionKernelError(RuntimeError):
    """The pinned R kernel failed (surfaced, never swallowed into empty results)."""


def run_production_kernel(
    space: FamilySpace,
    queries: Mapping[str, np.ndarray],
    *,
    rscript: str | None = None,
    max_cores: int | None = None,
) -> tuple[dict[str, dict[str, tuple[float, float, int]]], dict[str, str]]:
    """Run the unchanged production kernel; returns ``({name: {term: (p, padj, k)}}, provenance)``.

    ``max_cores`` sets ``$MPD_MAX_CORES`` for the call (the R driver still clamps
    to ``detectCores() - 4``). Raises :class:`ProductionKernelError` on any
    engine failure.
    """
    gene_sets = {
        name: [space.genes[i] for i in np.asarray(idx, dtype=np.int64)]
        for name, idx in queries.items()
    }
    term2gene = space.term2gene()
    universe = list(space.genes)
    prev = os.environ.get(cp.MAX_CORES_ENV)
    try:
        if max_cores is not None:
            os.environ[cp.MAX_CORES_ENV] = str(int(max_cores))
        result, prov = cp.run_ora_batch_with_provenance(
            gene_sets, term2gene, universe, term2name=None, rscript=rscript
        )
    except RuntimeError as exc:  # the bridge surfaces Rscript failures as RuntimeError
        raise ProductionKernelError(str(exc)) from exc
    finally:
        if max_cores is not None:
            if prev is None:
                os.environ.pop(cp.MAX_CORES_ENV, None)
            else:
                os.environ[cp.MAX_CORES_ENV] = prev
    try:
        n_err = int(prov.get("n_worker_errors", "0") or "0")
    except ValueError:
        n_err = 0
    if n_err:
        raise ProductionKernelError(
            f"production kernel: {n_err} gene set(s) failed in workers: "
            f"{prov.get('worker_error_examples', '')}"
        )
    return result, dict(prov)


@dataclass
class AgreementReport:
    n_queries: int
    n_term_rows: int
    term_set_mismatches: int
    overlap_mismatches: int
    max_rel_diff_p: float
    max_rel_diff_padj: float
    max_abs_diff_p: float
    max_abs_diff_padj: float
    tolerance: float
    passed: bool
    examples: list[dict[str, object]]

    def as_dict(self) -> dict[str, object]:
        return {
            "n_queries": self.n_queries,
            "n_term_rows": self.n_term_rows,
            "term_set_mismatches": self.term_set_mismatches,
            "overlap_mismatches": self.overlap_mismatches,
            "max_rel_diff_p": self.max_rel_diff_p,
            "max_rel_diff_padj": self.max_rel_diff_padj,
            "max_abs_diff_p": self.max_abs_diff_p,
            "max_abs_diff_padj": self.max_abs_diff_padj,
            "tolerance_rel": self.tolerance,
            "passed": self.passed,
            "examples": self.examples,
        }


def _rel(a: float, b: float) -> float:
    d = abs(a - b)
    scale = max(abs(a), abs(b))
    return d / scale if scale > 0 else 0.0


def agreement_report(
    space: FamilySpace,
    names: Sequence[str],
    r_result: Mapping[str, Mapping[str, tuple[float, float, int]]],
    py: OraBatch,
    *,
    tolerance: float = 1e-9,
    max_examples: int = 10,
) -> AgreementReport:
    """Compare production-kernel output with the independent kernel, query by query."""
    terms = space.terms
    n_rows = 0
    set_mm = 0
    ov_mm = 0
    max_rp = 0.0
    max_rq = 0.0
    max_ap = 0.0
    max_aq = 0.0
    examples: list[dict[str, object]] = []
    for i, name in enumerate(names):
        r_terms = r_result.get(name, {})
        py_terms = {
            terms[t]: (float(p), float(q), int(k))
            for t, p, q, k in zip(py.hits[i], py.p[i], py.padj_overlap[i], py.k[i], strict=True)
        }
        if set(r_terms) != set(py_terms):
            set_mm += 1
            if len(examples) < max_examples:
                only_r = sorted(set(r_terms) - set(py_terms))[:3]
                only_py = sorted(set(py_terms) - set(r_terms))[:3]
                examples.append(
                    {"query": name, "kind": "term_set", "only_r": only_r, "only_python": only_py}
                )
        for term, (rp, rq, rk) in r_terms.items():
            if term not in py_terms:
                continue
            pp, pq, pk = py_terms[term]
            n_rows += 1
            if rk != pk:
                ov_mm += 1
                if len(examples) < max_examples:
                    examples.append(
                        {"query": name, "term": term, "kind": "overlap", "r": rk, "python": pk}
                    )
            rel_p = _rel(rp, pp)
            rel_q = _rel(rq, pq)
            max_rp = max(max_rp, rel_p)
            max_rq = max(max_rq, rel_q)
            max_ap = max(max_ap, abs(rp - pp))
            max_aq = max(max_aq, abs(rq - pq))
            if (rel_p > tolerance or rel_q > tolerance) and len(examples) < max_examples:
                examples.append(
                    {
                        "query": name,
                        "term": term,
                        "kind": "pvalue",
                        "r_p": rp,
                        "python_p": pp,
                        "r_padj": rq,
                        "python_padj": pq,
                    }
                )
    passed = set_mm == 0 and ov_mm == 0 and max_rp <= tolerance and max_rq <= tolerance
    return AgreementReport(
        n_queries=len(names),
        n_term_rows=n_rows,
        term_set_mismatches=set_mm,
        overlap_mismatches=ov_mm,
        max_rel_diff_p=max_rp,
        max_rel_diff_padj=max_rq,
        max_abs_diff_p=max_ap,
        max_abs_diff_padj=max_aq,
        tolerance=tolerance,
        passed=passed,
        examples=examples,
    )
