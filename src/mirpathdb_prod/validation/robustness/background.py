"""Stage 3: background sensitivity on the REAL queries (paired 2 x 2 design).

For one stratum (organism x family x cumulative tier) the alternative background
is the union of the normalised targets of ALL eligible miRNAs of that tier,
intersected with the annotation background (it is a subset of it by
construction, since normalisation already restricted every target set to the
family's symbol space). Pathway memberships are restricted to that background
and every statistic is recomputed for every eligible miRNA under

    {annotation background, target-union background} x {overlap-only BH, fixed-family BH}

so background effects and correction effects can be told apart. Each query is
asserted to be preserved under the target-union background (``n`` unchanged);
because ``k`` is unchanged as well, the overlap-only BH family is identical
under both backgrounds and the background acts only through ``N`` and ``K``.

A background-dependent change is reported as such; it does not quantify or
isolate 3'UTR bias.
"""

from __future__ import annotations

import json
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from mirpathdb_prod.validation.robustness.inputs import FamilySpace, Stratum
from mirpathdb_prod.validation.robustness.kernel import OraBatch, fold_enrichment, ora_batch
from mirpathdb_prod.validation.robustness.r_kernel import agreement_report, run_production_kernel

CONDITIONS = (
    "annotation_overlap",
    "annotation_fixed",
    "target_union_overlap",
    "target_union_fixed",
)

PAIRS = {
    "background_effect_overlap": ("annotation_overlap", "target_union_overlap"),
    "background_effect_fixed": ("annotation_fixed", "target_union_fixed"),
    "correction_effect_annotation": ("annotation_overlap", "annotation_fixed"),
    "correction_effect_target_union": ("target_union_overlap", "target_union_fixed"),
}


def _sig(hits: np.ndarray, padj: np.ndarray, alpha: float) -> set[int]:
    return set(hits[padj <= alpha].tolist())


def compare_sets(a: set[int], b: set[int]) -> dict[str, float]:
    both = len(a & b)
    union = len(a | b)
    return {
        "n_a": len(a),
        "n_b": len(b),
        "n_both": both,
        "gains": len(b - a),
        "losses": len(a - b),
        "jaccard": both / union if union else 1.0,
    }


def _condition_sets(py_u: OraBatch, py_a: OraBatch, i: int, alpha: float) -> dict[str, set[int]]:
    return {
        "annotation_overlap": _sig(py_u.hits[i], py_u.padj_overlap[i], alpha),
        "annotation_fixed": _sig(py_u.hits[i], py_u.padj_fixed[i], alpha),
        "target_union_overlap": _sig(py_a.hits[i], py_a.padj_overlap[i], alpha),
        "target_union_fixed": _sig(py_a.hits[i], py_a.padj_fixed[i], alpha),
    }


def build_target_union_space(
    stratum: Stratum,
) -> tuple[FamilySpace, dict[str, np.ndarray], list[str]]:
    """Restrict the family space to the tier-wide target union; remap the queries.

    Returns ``(alt_space, alt_queries, violations)`` where ``violations`` lists
    miRNAs whose query is NOT fully preserved (must be empty by construction).
    """
    space = stratum.space
    union = stratum.target_union()
    alt = space.restrict_background([space.genes[i] for i in union], "target_union")
    alt_index = alt.gene_index()
    alt_queries: dict[str, np.ndarray] = {}
    violations: list[str] = []
    for m, idx in stratum.queries.items():
        genes = [space.genes[i] for i in idx]
        mapped = [alt_index[g] for g in genes if g in alt_index]
        if len(mapped) != len(genes):
            violations.append(m)
        alt_queries[m] = np.array(sorted(mapped), dtype=np.int64)
    return alt, alt_queries, violations


def run_background_stratum(
    stratum: Stratum,
    *,
    alphas: Sequence[float],
    out_dir: Path,
    figure_mirnas: Mapping[str, str],
    figure_terms: Sequence[str],
    cohort: set[str] | None = None,
    r_kernel: Mapping[str, object] | None = None,
    force: bool = False,
    log=print,
) -> dict[str, object]:
    """Run (or resume) the 2 x 2 background/correction comparison for one stratum."""
    key = stratum.key
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / f"{key}_summary.json"
    if summary_path.exists() and not force:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("status") == "complete":
            log(f"[background] {key}: complete (resume)")
            return summary
    space = stratum.space
    names = sorted(stratum.queries)
    if not names:
        summary = {"key": key, "status": "empty", "note": "no eligible miRNA"}
        summary_path.write_text(json.dumps(summary, indent=1), encoding="utf-8")
        return summary
    t0 = time.perf_counter()
    alt, alt_queries, violations = build_target_union_space(stratum)
    q_u = [stratum.queries[m] for m in names]
    q_a = [alt_queries[m] for m in names]
    py_u = ora_batch(space, q_u)
    py_a = ora_batch(alt, q_a)
    if not np.array_equal(py_u.n, py_a.n):
        violations = sorted(set(violations) | {names[i] for i in np.flatnonzero(py_u.n != py_a.n)})

    # ---- per-query rows ----------------------------------------------------
    n_sizes = py_u.n
    quart_edges = np.percentile(n_sizes, [25, 50, 75]) if n_sizes.size else np.zeros(3)
    rows = []
    for i, m in enumerate(names):
        row: dict[str, object] = {
            "mirna": m,
            "figure_example": figure_mirnas.get(m, ""),
            "in_shared_cohort": bool(cohort and m in cohort),
            "n": int(py_u.n[i]),
            "n_target_union": int(py_a.n[i]),
            "size_quartile": int(np.searchsorted(quart_edges, py_u.n[i], side="right") + 1),
            "N_annotation": int(space.N),
            "N_target_union": int(alt.N),
            "m_overlap": int(py_u.hits[i].size),
            "m_overlap_target_union": int(py_a.hits[i].size),
            "m_full_annotation": int(space.n_terms),
            "m_full_target_union": int(alt.n_terms),
            "query_preserved": m not in violations,
        }
        # p-value shift over the (identical) overlapping terms
        if py_u.hits[i].size and np.array_equal(py_u.hits[i], _map_hits(alt, space, py_a.hits[i])):
            with np.errstate(divide="ignore"):
                d = np.log10(np.maximum(py_a.p[i], 1e-300)) - np.log10(
                    np.maximum(py_u.p[i], 1e-300)
                )
            row["median_delta_log10p_target_union_minus_annotation"] = float(np.median(d))
            row["fraction_terms_p_increased"] = float((d > 0).mean())
            row["fraction_terms_p_decreased"] = float((d < 0).mean())
        else:
            row["median_delta_log10p_target_union_minus_annotation"] = float("nan")
            row["fraction_terms_p_increased"] = float("nan")
            row["fraction_terms_p_decreased"] = float("nan")
        for a in alphas:
            sets = _condition_sets(py_u, py_a, i, a)
            for cond in CONDITIONS:
                row[f"n_sig_{cond}_{a:g}"] = len(sets[cond])
            for pair, (c1, c2) in PAIRS.items():
                cmp_ = compare_sets(sets[c1], sets[c2])
                for stat in ("n_both", "gains", "losses", "jaccard"):
                    row[f"{pair}_{stat}_{a:g}"] = cmp_[stat]
        rows.append(row)
    per_query = pd.DataFrame(rows)
    per_query.to_csv(out_dir / f"{key}_queries.tsv", sep="\t", index=False)

    # ---- per-term rows for the figure miRNAs -------------------------------
    fig_rows = []
    fig_present = [m for m in names if m in figure_mirnas]
    fterms = set(figure_terms)
    for m in fig_present:
        i = names.index(m)
        fig_rows.extend(_term_rows(space, alt, py_u, py_a, i, m, figure_mirnas[m], fterms))
    fig_df = pd.DataFrame(fig_rows)
    fig_df.to_csv(out_dir / f"{key}_figure_terms.tsv", sep="\t", index=False)

    # ---- optional production-kernel run on the figure miRNAs (both backgrounds)
    r_block: dict[str, object] = {}
    if (r_kernel or {}).get("enabled") and fig_present:
        r_block = _production_on_figure(
            space,
            alt,
            stratum,
            alt_queries,
            fig_present,
            py_u,
            py_a,
            names,
            r_kernel,
            out_dir,
            key,
            log,
        )  # type: ignore[arg-type]

    # ---- summary -------------------------------------------------------------
    summary: dict[str, object] = {
        "key": key,
        "organism": stratum.organism,
        "family": stratum.family,
        "tier": stratum.tier,
        "status": "complete",
        "wall_s": round(time.perf_counter() - t0, 2),
        "n_queries": len(names),
        "n_queries_in_shared_cohort": int(per_query["in_shared_cohort"].sum()),
        "N_annotation": int(space.N),
        "N_target_union": int(alt.N),
        "m_full_annotation": int(space.n_terms),
        "m_full_target_union": int(alt.n_terms),
        "terms_unavailable_under_target_union": len(alt.dropped_terms),
        "terms_unavailable_examples": list(alt.dropped_terms[:10]),
        "query_preservation_violations": violations,
        "figure_mirnas_present": fig_present,
        "figure_terms_requested": list(figure_terms),
        "figure_terms_present": [t for t in figure_terms if t in space.term_index()],
        "aggregate": _aggregate(per_query, alphas),
        "aggregate_by_size_quartile": {
            str(q): _aggregate(per_query[per_query["size_quartile"] == q], alphas)
            for q in sorted(per_query["size_quartile"].unique())
        },
        "aggregate_shared_cohort": _aggregate(per_query[per_query["in_shared_cohort"]], alphas)
        if cohort
        else {},
        "production_kernel_figure_check": r_block,
    }
    summary_path.write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    log(
        f"[background] {key}: {len(names)} queries, N {space.N}->{alt.N}, "
        f"m_full {space.n_terms}->{alt.n_terms}, "
        f"unavailable {len(alt.dropped_terms)}, violations {len(violations)}, {summary['wall_s']}s"
    )
    return summary


def _map_hits(alt: FamilySpace, space: FamilySpace, alt_hits: np.ndarray) -> np.ndarray:
    tindex = space.term_index()
    return np.array(sorted(tindex[alt.terms[j]] for j in alt_hits), dtype=np.int64)


def _term_rows(
    space: FamilySpace,
    alt: FamilySpace,
    py_u: OraBatch,
    py_a: OraBatch,
    i: int,
    mirna: str,
    example: str,
    figure_terms: set[str],
) -> list[dict[str, object]]:
    alt_by_term = {alt.terms[j]: pos for pos, j in enumerate(py_a.hits[i])}
    n = int(py_u.n[i])
    rows = []
    fe_u = fold_enrichment(py_u.k[i], n, py_u.K[i], space.N)
    for pos, j in enumerate(py_u.hits[i]):
        term = space.terms[j]
        apos = alt_by_term.get(term)
        row: dict[str, object] = {
            "mirna": mirna,
            "figure_example": example,
            "term": term,
            "is_figure_term": term in figure_terms,
            "k": int(py_u.k[i][pos]),
            "n": n,
            "K_annotation": int(py_u.K[i][pos]),
            "N_annotation": int(space.N),
            "fold_enrichment_annotation": float(fe_u[pos]),
            "p_annotation": float(py_u.p[i][pos]),
            "padj_annotation_overlap": float(py_u.padj_overlap[i][pos]),
            "padj_annotation_fixed": float(py_u.padj_fixed[i][pos]),
        }
        if apos is not None:
            K_a = int(py_a.K[i][apos])
            row.update(
                {
                    "K_target_union": K_a,
                    "N_target_union": int(alt.N),
                    "fold_enrichment_target_union": float(
                        fold_enrichment(np.array([py_a.k[i][apos]]), n, np.array([K_a]), alt.N)[0]
                    ),
                    "p_target_union": float(py_a.p[i][apos]),
                    "padj_target_union_overlap": float(py_a.padj_overlap[i][apos]),
                    "padj_target_union_fixed": float(py_a.padj_fixed[i][apos]),
                }
            )
        rows.append(row)
    return rows


def _aggregate(df: pd.DataFrame, alphas: Sequence[float]) -> dict[str, object]:
    out: dict[str, object] = {"n_queries": int(len(df))}
    if df.empty:
        return out
    out["median_delta_log10p"] = float(
        df["median_delta_log10p_target_union_minus_annotation"].median()
    )
    out["mean_fraction_terms_p_increased"] = float(df["fraction_terms_p_increased"].mean())
    for a in alphas:
        block: dict[str, object] = {}
        for cond in CONDITIONS:
            col = df[f"n_sig_{cond}_{a:g}"]
            block[f"mean_n_sig_{cond}"] = float(col.mean())
            block[f"median_n_sig_{cond}"] = float(col.median())
            block[f"fraction_queries_any_sig_{cond}"] = float((col > 0).mean())
        for pair in PAIRS:
            j = df[f"{pair}_jaccard_{a:g}"]
            block[f"{pair}_median_jaccard"] = float(j.median())
            block[f"{pair}_mean_jaccard"] = float(j.mean())
            block[f"{pair}_fraction_identical"] = float((j >= 1.0).mean())
            block[f"{pair}_mean_gains"] = float(df[f"{pair}_gains_{a:g}"].mean())
            block[f"{pair}_mean_losses"] = float(df[f"{pair}_losses_{a:g}"].mean())
        out[f"alpha_{a:g}"] = block
    return out


def _production_on_figure(
    space, alt, stratum, alt_queries, fig_present, py_u, py_a, names, r_kernel, out_dir, key, log
) -> dict[str, object]:
    """Run the unchanged production kernel on the figure miRNAs under both backgrounds."""
    out: dict[str, object] = {}
    idx = [names.index(m) for m in fig_present]
    for label, sp_, queries, py in (
        ("annotation", space, {m: stratum.queries[m] for m in fig_present}, py_u),
        ("target_union", alt, {m: alt_queries[m] for m in fig_present}, py_a),
    ):
        path = out_dir / f"{key}_figure_{label}_r_kernel.tsv"
        try:
            result, prov = run_production_kernel(
                sp_,
                queries,
                rscript=r_kernel.get("rscript") or None,
                max_cores=int(r_kernel.get("max_cores") or 0) or None,
            )
        except Exception as exc:  # surfaced in the summary, never silently dropped
            out[label] = {"status": "failed", "error": str(exc)[:2000]}
            log(f"[background] {key}: production kernel FAILED on figure miRNAs ({label}): {exc}")
            continue
        sub = OraBatch(
            N=py.N,
            m_full=py.m_full,
            n=py.n[idx],
            hits=[py.hits[i] for i in idx],
            k=[py.k[i] for i in idx],
            K=[py.K[i] for i in idx],
            p=[py.p[i] for i in idx],
            padj_overlap=[py.padj_overlap[i] for i in idx],
            padj_fixed=[py.padj_fixed[i] for i in idx],
        )
        agree = agreement_report(
            sp_, fig_present, result, sub, tolerance=float(r_kernel.get("tolerance_rel", 1e-9))
        )
        with path.open("w", encoding="utf-8") as fh:
            fh.write("mirna\tterm\tp_value\tp_adjust\toverlap\n")
            for m in fig_present:
                for term, (p, q, k) in sorted(result.get(m, {}).items()):
                    fh.write(f"{m}\t{term}\t{p!r}\t{q!r}\t{k}\n")
        out[label] = {
            "status": "complete",
            "provenance": prov,
            "agreement": agree.as_dict(),
            "file": path.name,
        }
        log(
            f"[background] {key}: production kernel on figure miRNAs ({label}) "
            f"agreement passed={agree.passed}"
        )
    return out
