"""Stage 2: uniform-null calibration of overlap-only BH vs fixed-family BH.

For one (organism, family, tier) stratum and one query size ``n`` (a percentile
of the actual per-miRNA query sizes), draw ``reps`` uniform samples of ``n``
background genes without replacement, score every sample with the independent
kernel and -- when enabled -- with the unchanged production kernel on the SAME
samples, and summarise:

* fraction of replicates with ANY rejection at each adjusted-P threshold, with
  exact (Clopper-Pearson) binomial confidence intervals -- under this complete
  null that fraction IS the FDR of the procedure (every rejection is false);
* the distribution of the number of rejections per replicate;
* the mean fraction of eligible terms rejected (a different statistic, reported
  separately and never conflated with the FDR);
* per-term rejection frequency across replicates;
* retention at the storage bound (0.25), reported separately.

Nothing here measures biological false discoveries among real associations.
"""

from __future__ import annotations

import json
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

from mirpathdb_prod.validation.robustness.inputs import FamilySpace, Stratum
from mirpathdb_prod.validation.robustness.kernel import OraBatch, bh_fixed_family, ora_batch
from mirpathdb_prod.validation.robustness.r_kernel import agreement_report, run_production_kernel
from mirpathdb_prod.validation.robustness.sampling import config_rng, uniform_samples


def clopper_pearson(x: int, n: int, level: float = 0.95) -> tuple[float, float]:
    """Exact binomial confidence interval for ``x`` successes in ``n`` trials."""
    if n <= 0:
        return (float("nan"), float("nan"))
    a = (1.0 - level) / 2.0
    lo = 0.0 if x == 0 else float(beta.ppf(a, x, n - x + 1))
    hi = 1.0 if x == n else float(beta.ppf(1.0 - a, x + 1, n - x))
    return (lo, hi)


def replicate_table(
    batch: OraBatch, alphas: Sequence[float], retention_alpha: float
) -> pd.DataFrame:
    """Per-replicate rejection counts for both correction families."""
    rows = []
    for i in range(batch.size):
        p = batch.p[i]
        po = batch.padj_overlap[i]
        pf = batch.padj_fixed[i]
        row: dict[str, object] = {
            "rep": i,
            "n": int(batch.n[i]),
            "m_overlap": int(batch.hits[i].size),
            "min_p": float(p.min()) if p.size else 1.0,
            "min_padj_overlap": float(po.min()) if po.size else 1.0,
            "min_padj_fixed": float(pf.min()) if pf.size else 1.0,
        }
        for a in alphas:
            row[f"rej_overlap_{a:g}"] = int((po <= a).sum())
            row[f"rej_fixed_{a:g}"] = int((pf <= a).sum())
        row[f"ret_overlap_{retention_alpha:g}"] = int((po <= retention_alpha).sum())
        row[f"ret_fixed_{retention_alpha:g}"] = int((pf <= retention_alpha).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def summarise_replicates(
    table: pd.DataFrame, m_full: int, alphas: Sequence[float], retention_alpha: float
) -> dict[str, object]:
    reps = int(len(table))
    out: dict[str, object] = {"replicates": reps, "m_full": int(m_full)}
    if reps == 0:
        return out
    out["m_overlap_mean"] = float(table["m_overlap"].mean())
    out["m_overlap_min"] = int(table["m_overlap"].min())
    out["m_overlap_max"] = int(table["m_overlap"].max())
    for fam in ("overlap", "fixed"):
        block: dict[str, object] = {}
        for a in alphas:
            col = table[f"rej_{fam}_{a:g}"]
            x = int((col > 0).sum())
            lo, hi = clopper_pearson(x, reps)
            block[f"alpha_{a:g}"] = {
                "any_rejection_fraction": x / reps,
                "any_rejection_ci95": [lo, hi],
                "replicates_with_rejection": x,
                "mean_rejections": float(col.mean()),
                "mean_fraction_eligible_terms_rejected": float(col.mean() / m_full)
                if m_full
                else float("nan"),
                "rejection_count_quantiles": {
                    q: float(np.percentile(col, q)) for q in (50, 90, 95, 99, 100)
                },
                "rejection_count_histogram": {
                    str(int(k)): int(v) for k, v in col.value_counts().sort_index().items()
                },
            }
        rc = table[f"ret_{fam}_{retention_alpha:g}"]
        xr = int((rc > 0).sum())
        lo, hi = clopper_pearson(xr, reps)
        block[f"retention_{retention_alpha:g}"] = {
            "any_retained_fraction": xr / reps,
            "any_retained_ci95": [lo, hi],
            "mean_retained": float(rc.mean()),
        }
        out[fam] = block
    return out


def term_rejection_frequency(
    space: FamilySpace, batch: OraBatch, alphas: Sequence[float]
) -> pd.DataFrame:
    """Per-term frequency of being hit (k >= 1) and rejected, across replicates."""
    reps = batch.size
    m_full = space.n_terms
    hit = np.zeros(m_full, dtype=np.int64)
    rej: dict[str, np.ndarray] = {}
    for fam in ("overlap", "fixed"):
        for a in alphas:
            rej[f"{fam}_{a:g}"] = np.zeros(m_full, dtype=np.int64)
    for i in range(reps):
        hit[batch.hits[i]] += 1
        for a in alphas:
            rej[f"overlap_{a:g}"][batch.rejected_terms(i, a, "overlap")] += 1
            rej[f"fixed_{a:g}"][batch.rejected_terms(i, a, "fixed")] += 1
    any_rej = np.zeros(m_full, dtype=bool)
    for arr in rej.values():
        any_rej |= arr > 0
    idx = np.flatnonzero(any_rej)
    df = pd.DataFrame(
        {
            "term": [space.terms[j] for j in idx],
            "K": space.term_sizes[idx],
            "hit_frequency": hit[idx] / max(reps, 1),
        }
    )
    for key, arr in rej.items():
        df[f"reject_freq_{key}"] = arr[idx] / max(reps, 1)
    return df.sort_values(
        ["reject_freq_overlap_" + f"{alphas[-1]:g}", "term"], ascending=[False, True]
    )


def batch_from_production(
    space: FamilySpace,
    names: Sequence[str],
    result: Mapping[str, Mapping[str, tuple[float, float, int]]],
    n: Sequence[int],
) -> OraBatch:
    """Wrap production-kernel output in an :class:`OraBatch` (adds fixed-family BH from raw P)."""
    tindex = space.term_index()
    hits, ks, Ks, ps, po, pf = [], [], [], [], [], []
    for name in names:
        rows = result.get(name, {})
        terms = sorted(rows, key=lambda t: tindex[t])
        h = np.array([tindex[t] for t in terms], dtype=np.int64)
        p = np.array([rows[t][0] for t in terms], dtype=np.float64)
        q = np.array([rows[t][1] for t in terms], dtype=np.float64)
        k = np.array([rows[t][2] for t in terms], dtype=np.int64)
        hits.append(h)
        ks.append(k)
        Ks.append(space.term_sizes[h] if h.size else np.zeros(0, dtype=np.int64))
        ps.append(p)
        po.append(q)
        pf.append(bh_fixed_family(p, space.n_terms))
    return OraBatch(
        N=space.N,
        m_full=space.n_terms,
        n=np.asarray(n, dtype=np.int64),
        hits=hits,
        k=ks,
        K=Ks,
        p=ps,
        padj_overlap=po,
        padj_fixed=pf,
    )


def run_calibration_config(
    stratum: Stratum,
    *,
    size_label: str,
    n: int,
    reps: int,
    seed: int,
    alphas: Sequence[float],
    retention_alpha: float,
    out_dir: Path,
    r_kernel: Mapping[str, object] | None = None,
    force: bool = False,
    force_python: bool = False,
    log=print,
) -> dict[str, object]:
    """Run (or resume) one calibration configuration; returns its summary dict.

    ``force`` recomputes everything; ``force_python`` recomputes only the
    independent kernel, its replicate tables and the agreement audit while
    re-using cached production-kernel outputs (the expensive part).
    """
    space = stratum.space
    key = f"{stratum.key}_{size_label}_n{n}"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / f"{key}.json"
    if summary_path.exists() and not force and not force_python:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("status") == "complete" and (
            not (r_kernel or {}).get("enabled")
            or summary.get("production_kernel", {}).get("status") == "complete"
        ):
            log(f"[calibration] {key}: complete (resume)")
            return summary
    else:
        summary = {}

    rng = config_rng(seed, "calibration", stratum.organism, stratum.family, stratum.tier, n)
    samples = uniform_samples(rng, space.N, n, reps)
    names = [f"rep_{i:04d}" for i in range(reps)]

    t0 = time.perf_counter()
    py = ora_batch(space, samples)
    t_py = time.perf_counter() - t0
    table = replicate_table(py, alphas, retention_alpha)
    table.to_csv(out_dir / f"{key}_replicates.tsv", sep="\t", index=False)
    term_rejection_frequency(space, py, alphas).to_csv(
        out_dir / f"{key}_terms.tsv", sep="\t", index=False
    )
    summary.update(
        {
            "key": key,
            "organism": stratum.organism,
            "family": stratum.family,
            "tier": stratum.tier,
            "size_label": size_label,
            "n": int(n),
            "N": int(space.N),
            "m_full": int(space.n_terms),
            "replicates": int(reps),
            "seed": int(seed),
            "sampling": "uniform without replacement from the annotation background",
            "independent_kernel": {
                "status": "complete",
                "wall_s": round(t_py, 3),
                **summarise_replicates(table, space.n_terms, alphas, retention_alpha),
            },
            "status": "complete",
        }
    )
    log(f"[calibration] {key}: independent kernel {reps} reps in {t_py:.1f}s")

    if (r_kernel or {}).get("enabled"):
        r_dir = out_dir / "r_kernel"
        r_dir.mkdir(exist_ok=True)
        r_path = r_dir / f"{key}.tsv"
        prov_path = r_dir / f"{key}_provenance.json"
        if r_path.exists() and prov_path.exists() and not force:
            result = _read_r_long(r_path)
            prov = json.loads(prov_path.read_text(encoding="utf-8"))
            t_r = float(prov.get("_wall_s", float("nan")))
            log(f"[calibration] {key}: production kernel loaded (resume)")
        else:
            t0 = time.perf_counter()
            result, prov = run_production_kernel(
                space,
                dict(zip(names, samples, strict=True)),
                rscript=r_kernel.get("rscript") or None,  # type: ignore[union-attr]
                max_cores=int(r_kernel.get("max_cores") or 0) or None,  # type: ignore[union-attr]
            )
            t_r = time.perf_counter() - t0
            prov = dict(prov, _wall_s=round(t_r, 3))
            _write_r_long(r_path, names, result)
            prov_path.write_text(json.dumps(prov, indent=1), encoding="utf-8")
            log(f"[calibration] {key}: production kernel {reps} reps in {t_r:.1f}s")
        rb = batch_from_production(space, names, result, [len(s) for s in samples])
        r_table = replicate_table(rb, alphas, retention_alpha)
        r_table.to_csv(r_dir / f"{key}_replicates.tsv", sep="\t", index=False)
        agree = agreement_report(
            space, names, result, py, tolerance=float(r_kernel.get("tolerance_rel", 1e-9))
        )  # type: ignore[union-attr]
        (r_dir / f"{key}_agreement.json").write_text(
            json.dumps(agree.as_dict(), indent=1), encoding="utf-8"
        )
        summary["production_kernel"] = {
            "status": "complete",
            "wall_s": round(t_r, 3),
            "provenance": prov,
            **summarise_replicates(r_table, space.n_terms, alphas, retention_alpha),
        }
        summary["kernel_agreement"] = agree.as_dict()
        log(
            f"[calibration] {key}: agreement passed={agree.passed} "
            f"max_rel_p={agree.max_rel_diff_p:.2e} "
            f"max_rel_padj={agree.max_rel_diff_padj:.2e} rows={agree.n_term_rows}"
        )
    summary_path.write_text(json.dumps(summary, indent=1), encoding="utf-8")
    return summary


def _write_r_long(
    path: Path, names: Sequence[str], result: Mapping[str, Mapping[str, tuple[float, float, int]]]
) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("name\tterm\tp_value\tp_adjust\toverlap\n")
        for name in names:
            for term, (p, q, k) in sorted(result.get(name, {}).items()):
                fh.write(f"{name}\t{term}\t{p!r}\t{q!r}\t{k}\n")


def _read_r_long(path: Path) -> dict[str, dict[str, tuple[float, float, int]]]:
    df = pd.read_csv(path, sep="\t", dtype={"name": str, "term": str})
    out: dict[str, dict[str, tuple[float, float, int]]] = {}
    for name, term, p, q, k in df.itertuples(index=False):
        out.setdefault(str(name), {})[str(term)] = (float(p), float(q), int(k))
    return out
