"""Stage 4: target-frequency-aware (degree-binned) sensitivity check.

Nulls: (1) size-matched resampling within predefined gene-degree bins, preserving the
focal query's count per bin (degree-matched); (2) uniform draws of the same size from the
whole background; (3) uniform draws of the same size from the target-union pool only
(genes with degree >= 1) -- the control that separates pool restriction from degree
matching. Details of (1): size-matched resampling within predefined gene-degree bins, preserving the
focal query's count per bin. Degree = number of eligible miRNAs of the tier
targeting the gene in the complete organism/tier network restricted to the
family background; bins are degree quantiles defined WITHOUT pathway labels;
focal-miRNA edges DO contribute to the degrees (documented). Sampling is
without replacement; genes with degree 0 are never drawn.

This answers a different question from uniform sampling ("does a random set of
genes with the same target-frequency profile produce the same rejections?"). It
does not establish independence from sequence features, assay biases or
genuine shared regulation; seed/family structure is NOT preserved (each focal
miRNA is resampled alone, siblings are not co-resampled); tiers are randomised
independently (no joint null preserving cumulative nesting). Monte-Carlo tail
probabilities from R draws have a floor of 1/(R+1) and are NOT adjusted P-values.
"""

from __future__ import annotations

import json
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from mirpathdb_prod.validation.robustness.calibration import clopper_pearson
from mirpathdb_prod.validation.robustness.inputs import Stratum
from mirpathdb_prod.validation.robustness.kernel import (
    OraBatch,
    fold_enrichment,
    ora_batch,
    ora_single,
)
from mirpathdb_prod.validation.robustness.sampling import (
    DegreeBins,
    achieved_degree_balance,
    config_rng,
    degree_matched_samples,
    network_degrees,
    uniform_pool_samples,
    uniform_samples,
)


def _any_rejection(batch: OraBatch, alphas: Sequence[float]) -> dict[str, object]:
    out: dict[str, object] = {}
    reps = batch.size
    for fam in ("overlap", "fixed"):
        for a in alphas:
            counts = batch.n_rejections(a, fam)
            x = int((counts > 0).sum())
            lo, hi = clopper_pearson(x, reps)
            out[f"{fam}_{a:g}"] = {
                "any_rejection_fraction": x / reps if reps else float("nan"),
                "any_rejection_ci95": [lo, hi],
                "mean_rejections": float(counts.mean()) if reps else float("nan"),
                "rejection_count_quantiles": {
                    q: float(np.percentile(counts, q)) for q in (50, 90, 99, 100)
                }
                if reps
                else {},
            }
    return out


def _term_null(batch: OraBatch, term_idx: int) -> np.ndarray:
    """Per-replicate overlap k for one term (0 when the term was not hit)."""
    ks = np.zeros(batch.size, dtype=np.int64)
    for i in range(batch.size):
        pos = np.searchsorted(batch.hits[i], term_idx)
        if pos < batch.hits[i].size and batch.hits[i][pos] == term_idx:
            ks[i] = batch.k[i][pos]
    return ks


def run_degree_stratum(
    stratum: Stratum,
    *,
    seed: int,
    alphas: Sequence[float],
    out_dir: Path,
    figure_mirnas: Mapping[str, str],
    figure_terms: Sequence[str],
    n_bins: int = 10,
    reps: int = 500,
    random_focal: int = 20,
    force: bool = False,
    log=print,
) -> dict[str, object]:
    """Run (or resume) the degree-matched null for the focal miRNAs of one stratum."""
    key = stratum.key
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / f"{key}_summary.json"
    if summary_path.exists() and not force:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("status") == "complete":
            log(f"[degree] {key}: complete (resume)")
            return summary
    space = stratum.space
    if not stratum.queries:
        summary = {"key": key, "status": "empty", "note": "no eligible miRNA"}
        summary_path.write_text(json.dumps(summary, indent=1), encoding="utf-8")
        return summary
    t0 = time.perf_counter()
    degrees = network_degrees(space, stratum.queries)
    bins = DegreeBins.from_degrees(degrees, n_bins=n_bins)
    pool = np.flatnonzero(degrees > 0)  # the tier's target-union pool inside the background
    fig_present = [m for m in sorted(stratum.queries) if m in figure_mirnas]
    others = [m for m in sorted(stratum.queries) if m not in figure_mirnas]
    rng_f = config_rng(seed, "degree_focal", key)
    random_pick = (
        sorted(rng_f.choice(others, size=min(random_focal, len(others)), replace=False).tolist())
        if others and random_focal > 0
        else []
    )
    focal = fig_present + random_pick
    tindex = space.term_index()
    fig_term_idx = {t: tindex[t] for t in figure_terms if t in tindex}

    focal_rows = []
    fig_term_rows = []
    for m in focal:
        q = stratum.queries[m]
        n = int(q.size)
        diag = bins.diagnostics(q)
        rng_d = config_rng(seed, "degree", key, m)
        samples_d = degree_matched_samples(rng_d, bins, q, reps)
        rng_u = config_rng(seed, "degree_uniform", key, m)
        samples_u = uniform_samples(rng_u, space.N, n, reps)
        # Control null: uniform draws restricted to the target-union pool (degree >= 1)
        # without degree matching, to separate pool restriction from degree matching.
        rng_p = config_rng(seed, "degree_pool_uniform", key, m)
        samples_p = uniform_pool_samples(rng_p, pool, n, reps)
        obs = ora_single(space, q)
        py_d = ora_batch(space, samples_d)
        py_u = ora_batch(space, samples_u)
        py_p = ora_batch(space, samples_p)
        balance = achieved_degree_balance(bins, q, samples_d)
        rec: dict[str, object] = {
            "mirna": m,
            "figure_example": figure_mirnas.get(m, ""),
            "n": n,
            "observed": {
                f"n_sig_{fam}_{a:g}": int(obs.n_rejections(a, fam)[0])
                for fam in ("overlap", "fixed")
                for a in alphas
            },
            "degree_matched": _any_rejection(py_d, alphas),
            "uniform": _any_rejection(py_u, alphas),
            "pool_uniform": _any_rejection(py_p, alphas),
            "diagnostics": diag,
            "balance": balance,
        }
        focal_rows.append(rec)
        # Null distributions for the prespecified figure terms.
        obs_terms = obs.as_dict(0, space.terms)
        for term, j in fig_term_idx.items():
            k_obs = obs_terms.get(term, (1.0, 1.0, 0))[2]
            K = int(space.term_sizes[j])
            kd = _term_null(py_d, j)
            ku = _term_null(py_u, j)
            kp = _term_null(py_p, j)
            fe_obs = (
                float(fold_enrichment(np.array([k_obs]), n, np.array([K]), space.N)[0])
                if n
                else float("nan")
            )
            fig_term_rows.append(
                {
                    "mirna": m,
                    "figure_example": figure_mirnas.get(m, ""),
                    "term": term,
                    "K": K,
                    "n": n,
                    "k_observed": int(k_obs),
                    "fold_enrichment_observed": fe_obs,
                    "p_observed": obs_terms.get(term, (1.0, 1.0, 0))[0],
                    "padj_overlap_observed": obs_terms.get(term, (1.0, 1.0, 0))[1],
                    "k_null_degree_mean": float(kd.mean()),
                    "k_null_degree_q50": float(np.percentile(kd, 50)),
                    "k_null_degree_q95": float(np.percentile(kd, 95)),
                    "k_null_degree_max": int(kd.max()),
                    "mc_tail_degree": float((1 + (kd >= k_obs).sum()) / (reps + 1)),
                    "k_null_uniform_mean": float(ku.mean()),
                    "k_null_uniform_q50": float(np.percentile(ku, 50)),
                    "k_null_uniform_q95": float(np.percentile(ku, 95)),
                    "k_null_uniform_max": int(ku.max()),
                    "mc_tail_uniform": float((1 + (ku >= k_obs).sum()) / (reps + 1)),
                    "k_null_pool_mean": float(kp.mean()),
                    "k_null_pool_q50": float(np.percentile(kp, 50)),
                    "k_null_pool_q95": float(np.percentile(kp, 95)),
                    "k_null_pool_max": int(kp.max()),
                    "mc_tail_pool": float((1 + (kp >= k_obs).sum()) / (reps + 1)),
                    "mc_floor": 1.0 / (reps + 1),
                    "fold_enrichment_null_degree_q95": float(
                        np.percentile(fold_enrichment(kd, n, np.full(reps, K), space.N), 95)
                    )
                    if n
                    else float("nan"),
                }
            )
        (out_dir / f"{key}_{_safe(m)}.json").write_text(
            json.dumps(rec, indent=1, default=float), encoding="utf-8"
        )
        log(
            f"[degree] {key} {m}: n={n} bins={diag['n_bins']} "
            f"max_occ={diag['max_bin_occupancy']:.2f} any-rej(0.05) "
            f"degree={rec['degree_matched']['overlap_0.05']['any_rejection_fraction']:.3f} "
            f"uniform={rec['uniform']['overlap_0.05']['any_rejection_fraction']:.3f} "
            f"pool={rec['pool_uniform']['overlap_0.05']['any_rejection_fraction']:.3f}"
        )
    pd.DataFrame(fig_term_rows).to_csv(out_dir / f"{key}_figure_terms.tsv", sep="\t", index=False)
    flat = []
    for rec in focal_rows:
        row: dict[str, object] = {
            "mirna": rec["mirna"],
            "figure_example": rec["figure_example"],
            "n": rec["n"],
        }
        row.update({f"observed_{k}": v for k, v in rec["observed"].items()})  # type: ignore[union-attr]
        for null in ("degree_matched", "uniform", "pool_uniform"):
            for fk, block in rec[null].items():  # type: ignore[union-attr]
                row[f"{null}_any_rej_{fk}"] = block["any_rejection_fraction"]
                row[f"{null}_mean_rej_{fk}"] = block["mean_rejections"]
        d = rec["diagnostics"]
        row["max_bin_occupancy"] = d["max_bin_occupancy"]  # type: ignore[index]
        row["n_bins_nearly_fixed"] = d["n_bins_nearly_fixed"]  # type: ignore[index]
        row["expected_self_overlap_fraction"] = d["expected_self_overlap_fraction"]  # type: ignore[index]
        row["distinct_genes_over_query_size"] = rec["balance"]["distinct_genes_over_query_size"]  # type: ignore[index]
        flat.append(row)
    pd.DataFrame(flat).to_csv(out_dir / f"{key}_focal.tsv", sep="\t", index=False)

    rand_rows = [r for r in focal_rows if not r["figure_example"]]
    summary: dict[str, object] = {
        "key": key,
        "organism": stratum.organism,
        "family": stratum.family,
        "tier": stratum.tier,
        "status": "complete",
        "wall_s": round(time.perf_counter() - t0, 2),
        "replicates": reps,
        "n_bins_requested": n_bins,
        "n_bins_effective": bins.n_bins,
        "bin_degree_ranges": bins.bin_ranges(),
        "bin_populations": [int(m.size) for m in bins.members],
        "n_genes_degree0_never_sampled": int((degrees == 0).sum()),
        "pool_size": int(pool.size),
        "N": int(space.N),
        "network_mirnas": len(stratum.queries),
        "focal_edges_included_in_degree": True,
        "seed_family_structure_preserved": False,
        "tiers_randomised_independently": True,
        "figure_mirnas_present": fig_present,
        "random_focal": random_pick,
        "figure_terms_present": sorted(fig_term_idx),
        "aggregate_random_focal": {
            fk: {
                "mean_any_rejection_degree_matched": float(
                    np.mean([r["degree_matched"][fk]["any_rejection_fraction"] for r in rand_rows])
                ),  # type: ignore[index]
                "mean_any_rejection_uniform": float(
                    np.mean([r["uniform"][fk]["any_rejection_fraction"] for r in rand_rows])
                ),  # type: ignore[index]
                "mean_any_rejection_pool_uniform": float(
                    np.mean([r["pool_uniform"][fk]["any_rejection_fraction"] for r in rand_rows])
                ),  # type: ignore[index]
            }
            for fk in (f"{fam}_{a:g}" for fam in ("overlap", "fixed") for a in alphas)
        }
        if rand_rows
        else {},
        "n_focal_nearly_fixed": int(
            sum(1 for r in focal_rows if r["diagnostics"]["n_bins_nearly_fixed"] > 0)
        ),  # type: ignore[index]
    }
    summary_path.write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    log(f"[degree] {key}: {len(focal)} focal miRNAs in {summary['wall_s']}s")
    return summary


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)
