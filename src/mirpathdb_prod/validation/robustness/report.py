"""Stage 5: machine-readable summary tables and compact plots from the stage outputs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def _read_tsv(path: Path) -> pd.DataFrame:
    """Read a stage TSV; an empty file (a stratum with nothing to list) is an empty frame."""
    if path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _load_json_dir(d: Path, suffix: str) -> list[dict]:
    out = []
    for p in sorted(d.glob(f"*{suffix}")):
        try:
            out.append(json.loads(p.read_text(encoding="utf-8")))
        except json.JSONDecodeError:
            continue
    return out


def calibration_table(out_dir: Path) -> pd.DataFrame:
    rows = []
    for s in _load_json_dir(out_dir / "calibration", ".json"):
        if s.get("status") != "complete":
            continue
        for kernel in ("independent_kernel", "production_kernel"):
            blk = s.get(kernel)
            if not blk or blk.get("status") != "complete":
                continue
            for fam in ("overlap", "fixed"):
                for akey, v in blk.get(fam, {}).items():
                    if not akey.startswith("alpha_"):
                        continue
                    rows.append(
                        {
                            "key": s["key"],
                            "organism": s["organism"],
                            "family": s["family"],
                            "tier": s["tier"],
                            "size_label": s["size_label"],
                            "n": s["n"],
                            "N": s["N"],
                            "m_full": s["m_full"],
                            "replicates": s["replicates"],
                            "kernel": kernel,
                            "correction": fam,
                            "alpha": float(akey.split("_", 1)[1]),
                            "any_rejection_fraction": v["any_rejection_fraction"],
                            "ci95_low": v["any_rejection_ci95"][0],
                            "ci95_high": v["any_rejection_ci95"][1],
                            "mean_rejections": v["mean_rejections"],
                            "mean_fraction_eligible_terms_rejected": v[
                                "mean_fraction_eligible_terms_rejected"
                            ],
                            "rej_q50": v["rejection_count_quantiles"].get("50"),
                            "rej_q95": v["rejection_count_quantiles"].get("95"),
                            "rej_max": v["rejection_count_quantiles"].get("100"),
                            "m_overlap_mean": blk.get("m_overlap_mean"),
                            "wall_s": blk.get("wall_s"),
                            "kernel_agreement_passed": s.get("kernel_agreement", {}).get("passed"),
                        }
                    )
    return pd.DataFrame(rows)


def retention_table(out_dir: Path) -> pd.DataFrame:
    rows = []
    for s in _load_json_dir(out_dir / "calibration", ".json"):
        if s.get("status") != "complete":
            continue
        for kernel in ("independent_kernel", "production_kernel"):
            blk = s.get(kernel)
            if not blk or blk.get("status") != "complete":
                continue
            for fam in ("overlap", "fixed"):
                for rkey, v in blk.get(fam, {}).items():
                    if not rkey.startswith("retention_"):
                        continue
                    rows.append(
                        {
                            "key": s["key"],
                            "organism": s["organism"],
                            "family": s["family"],
                            "tier": s["tier"],
                            "size_label": s["size_label"],
                            "n": s["n"],
                            "kernel": kernel,
                            "correction": fam,
                            "retention_alpha": float(rkey.split("_", 1)[1]),
                            "any_retained_fraction": v["any_retained_fraction"],
                            "ci95_low": v["any_retained_ci95"][0],
                            "ci95_high": v["any_retained_ci95"][1],
                            "mean_retained": v["mean_retained"],
                        }
                    )
    return pd.DataFrame(rows)


def agreement_table(out_dir: Path) -> pd.DataFrame:
    rows = []
    for s in _load_json_dir(out_dir / "calibration", ".json"):
        a = s.get("kernel_agreement")
        if not a:
            continue
        rows.append(
            {
                "key": s["key"],
                "organism": s["organism"],
                "family": s["family"],
                "tier": s["tier"],
                "n": s["n"],
                **{k: v for k, v in a.items() if k != "examples"},
            }
        )
    for s in _load_json_dir(out_dir / "background", "_summary.json"):
        for label, blk in (s.get("production_kernel_figure_check") or {}).items():
            a = blk.get("agreement") if isinstance(blk, dict) else None
            if a:
                rows.append(
                    {
                        "key": f"{s['key']}_figure_{label}",
                        "organism": s["organism"],
                        "family": s["family"],
                        "tier": s["tier"],
                        "n": None,
                        **{k: v for k, v in a.items() if k != "examples"},
                    }
                )
    return pd.DataFrame(rows)


def background_table(out_dir: Path) -> pd.DataFrame:
    rows = []
    for s in _load_json_dir(out_dir / "background", "_summary.json"):
        if s.get("status") != "complete":
            continue
        base = {
            "key": s["key"],
            "organism": s["organism"],
            "family": s["family"],
            "tier": s["tier"],
            "n_queries": s["n_queries"],
            "n_queries_shared_cohort": s.get("n_queries_in_shared_cohort"),
            "N_annotation": s["N_annotation"],
            "N_target_union": s["N_target_union"],
            "m_full_annotation": s["m_full_annotation"],
            "m_full_target_union": s["m_full_target_union"],
            "terms_unavailable_under_target_union": s["terms_unavailable_under_target_union"],
            "n_preservation_violations": len(s.get("query_preservation_violations", [])),
            "median_delta_log10p": s["aggregate"].get("median_delta_log10p"),
        }
        for scope_name, agg in (
            ("all", s["aggregate"]),
            ("shared_cohort", s.get("aggregate_shared_cohort") or {}),
        ):
            for akey, blk in agg.items():
                if not akey.startswith("alpha_"):
                    continue
                rows.append(
                    {**base, "scope": scope_name, "alpha": float(akey.split("_", 1)[1]), **blk}
                )
        for q, agg in (s.get("aggregate_by_size_quartile") or {}).items():
            for akey, blk in agg.items():
                if not akey.startswith("alpha_"):
                    continue
                rows.append(
                    {
                        **base,
                        "scope": f"size_quartile_{q}",
                        "n_queries": agg.get("n_queries"),
                        "alpha": float(akey.split("_", 1)[1]),
                        **blk,
                    }
                )
    return pd.DataFrame(rows)


def degree_table(out_dir: Path) -> pd.DataFrame:
    frames = []
    for p in sorted((out_dir / "degree").glob("*_focal.tsv")):
        df = _read_tsv(p)
        key = p.name[: -len("_focal.tsv")]
        org, fam, tier = (
            key.split("_", 1)[0],
            key.split("_", 1)[1].rsplit("_", 1)[0],
            key.rsplit("_", 1)[1],
        )
        df.insert(0, "tier", tier)
        df.insert(0, "family", fam)
        df.insert(0, "organism", org)
        df.insert(0, "key", key)
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def figure_terms_table(out_dir: Path) -> pd.DataFrame:
    frames = []
    for p in sorted((out_dir / "degree").glob("*_figure_terms.tsv")):
        df = _read_tsv(p)
        if df.empty:
            continue
        key = p.name[: -len("_figure_terms.tsv")]
        df.insert(0, "key", key)
        frames.append(df)
    deg = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    frames = []
    for p in sorted((out_dir / "background").glob("*_figure_terms.tsv")):
        df = _read_tsv(p)
        if df.empty:
            continue
        df = df[df["is_figure_term"]]
        key = p.name[: -len("_figure_terms.tsv")]
        df.insert(0, "key", key)
        frames.append(df)
    bg = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if deg.empty or bg.empty:
        return deg if not deg.empty else bg
    return deg.merge(bg, on=["key", "mirna", "term"], how="outer", suffixes=("", "_bg"))


def paired_table(out_dir: Path) -> pd.DataFrame:
    """Old-versus-candidate on identical draws, per configuration and kernel.

    Every replicate table scores the same uniform draw under the production overlap-only BH
    (``rej_overlap_*``) and under BH over the full eligible family (``rej_fixed_*``; the
    candidate convention, numerator ``m_full``).  The paired quantity of interest is the
    fraction of draws rejecting under the production convention but not under the candidate,
    with a Clopper-Pearson interval; the reverse fraction is reported and must be zero because
    ``q_full >= q_overlap`` on every draw.  Rows with ``measure == "retention"`` do the same
    for the 0.25 storage bound.
    """
    from mirpathdb_prod.validation.robustness.calibration import clopper_pearson

    meta = {
        s["key"]: s
        for s in _load_json_dir(out_dir / "calibration", ".json")
        if s.get("status") == "complete"
    }
    rows = []
    for kernel, d in (
        ("independent_kernel", out_dir / "calibration"),
        ("production_kernel", out_dir / "calibration" / "r_kernel"),
    ):
        if not d.is_dir():
            continue
        for p in sorted(d.glob("*_replicates.tsv")):
            key = p.name[: -len("_replicates.tsv")]
            s = meta.get(key)
            if s is None:
                continue
            df = _read_tsv(p)
            if df.empty:
                continue
            reps = int(len(df))
            pairs = [
                ("rejection", a, f"rej_overlap_{a}", f"rej_fixed_{a}") for a in ("0.005", "0.05")
            ]
            pairs.append(("retention", "0.25", "ret_overlap_0.25", "ret_fixed_0.25"))
            for measure, alpha, co, cf in pairs:
                if co not in df.columns or cf not in df.columns:
                    continue
                prod = df[co].to_numpy() > 0
                cand = df[cf].to_numpy() > 0
                excess = int((prod & ~cand).sum())
                reverse = int((cand & ~prod).sum())
                lo, hi = clopper_pearson(excess, reps)
                rows.append(
                    {
                        "key": key,
                        "organism": s["organism"],
                        "family": s["family"],
                        "tier": s["tier"],
                        "size_label": s["size_label"],
                        "n": s["n"],
                        "m_full": s.get("m_full"),
                        "kernel": kernel,
                        "measure": measure,
                        "alpha": float(alpha),
                        "replicates": reps,
                        "production_any": float(prod.mean()),
                        "candidate_any": float(cand.mean()),
                        "excess_any": excess / reps,
                        "excess_ci95_low": lo,
                        "excess_ci95_high": hi,
                        "reverse_any": reverse / reps,
                        "production_mean_count": float(df[co].mean()),
                        "candidate_mean_count": float(df[cf].mean()),
                        "mean_paired_difference": float((df[co] - df[cf]).mean()),
                        "candidate_never_exceeds_production": bool((df[cf] <= df[co]).all()),
                    }
                )
    return pd.DataFrame(rows)


def write_summaries(out_dir: Path) -> dict[str, Path]:
    sdir = out_dir / "summary"
    sdir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for name, fn in (
        ("calibration", calibration_table),
        ("retention", retention_table),
        ("paired_old_vs_candidate", paired_table),
        ("kernel_agreement", agreement_table),
        ("background", background_table),
        ("degree", degree_table),
        ("figure_terms", figure_terms_table),
    ):
        df = fn(out_dir)
        path = sdir / f"{name}.tsv"
        df.to_csv(path, sep="\t", index=False)
        written[name] = path
    return written


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def _mpl():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def plot_calibration(out_dir: Path, plots_dir: Path) -> list[Path]:
    plt = _mpl()
    df = calibration_table(out_dir)
    paths: list[Path] = []
    if df.empty:
        return paths
    kernel = (
        "production_kernel" if (df["kernel"] == "production_kernel").any() else "independent_kernel"
    )
    d = df[df["kernel"] == kernel]
    families = sorted(d["family"].unique())
    alphas = sorted(d["alpha"].unique())
    fig, axes = plt.subplots(
        len(alphas), len(families), figsize=(6 * len(families), 3.6 * len(alphas)), squeeze=False
    )
    for r, a in enumerate(alphas):
        for c, fam in enumerate(families):
            ax = axes[r][c]
            sub = d[(d["alpha"] == a) & (d["family"] == fam)].copy()
            sub["x"] = (
                sub["organism"]
                + "/"
                + sub["tier"]
                + "/"
                + sub["size_label"]
                + " (n="
                + sub["n"].astype(str)
                + ")"
            )
            xs = sorted(sub["x"].unique())
            pos = {x: i for i, x in enumerate(xs)}
            for corr, colour, off in (("overlap", "#c0392b", -0.15), ("fixed", "#2471a3", 0.15)):
                ss = sub[sub["correction"] == corr]
                x = np.array([pos[v] for v in ss["x"]]) + off
                y = ss["any_rejection_fraction"].to_numpy()
                lo = y - ss["ci95_low"].to_numpy()
                hi = ss["ci95_high"].to_numpy() - y
                ax.errorbar(
                    x,
                    y,
                    yerr=[lo, hi],
                    fmt="o",
                    ms=4,
                    color=colour,
                    capsize=2,
                    label=f"{corr}-family BH",
                )
            ax.axhline(a, color="grey", ls="--", lw=1)
            ax.set_xticks(range(len(xs)))
            ax.set_xticklabels(xs, rotation=90, fontsize=6)
            ax.set_ylim(0, 1)
            ax.set_title(f"{fam}, threshold {a:g} ({kernel.replace('_', ' ')})", fontsize=9)
            ax.set_ylabel("P(any rejection) under uniform null")
            if r == 0 and c == 0:
                ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = plots_dir / f"calibration_any_rejection.{ext}"
        fig.savefig(p, dpi=150)
        paths.append(p)
    plt.close(fig)
    return paths


def plot_calibration_vs_n(out_dir: Path, plots_dir: Path) -> list[Path]:
    """Any-rejection fraction as a function of the query size (log axis), per family."""
    plt = _mpl()
    df = calibration_table(out_dir)
    paths: list[Path] = []
    if df.empty:
        return paths
    kernel = (
        "production_kernel" if (df["kernel"] == "production_kernel").any() else "independent_kernel"
    )
    d = df[df["kernel"] == kernel]
    families = sorted(d["family"].unique())
    alphas = sorted(d["alpha"].unique())
    fig, axes = plt.subplots(
        len(alphas), len(families), figsize=(6 * len(families), 3.8 * len(alphas)), squeeze=False
    )
    markers = {"hsa": "o", "mmu": "s", "rno": "^", "mml": "D"}
    for r, a in enumerate(alphas):
        for c, fam in enumerate(families):
            ax = axes[r][c]
            sub = d[(d["alpha"] == a) & (d["family"] == fam)]
            for org in sorted(sub["organism"].unique()):
                for corr, colour in (("overlap", "#c0392b"), ("fixed", "#2471a3")):
                    ss = sub[(sub["organism"] == org) & (sub["correction"] == corr)].sort_values(
                        "n"
                    )
                    if ss.empty:
                        continue
                    y = ss["any_rejection_fraction"].to_numpy()
                    lo = y - ss["ci95_low"].to_numpy()
                    hi = ss["ci95_high"].to_numpy() - y
                    ax.errorbar(
                        ss["n"],
                        y,
                        yerr=[lo, hi],
                        fmt=markers.get(org, "o") + "-",
                        ms=4,
                        lw=1,
                        color=colour,
                        alpha=0.85,
                        capsize=2,
                        label=f"{org} {corr}-family BH",
                    )
            ax.axhline(a, color="grey", ls="--", lw=1)
            ax.set_xscale("log")
            ax.set_ylim(-0.02, 1.0)
            ax.set_xlabel("query size n (genes in background)")
            ax.set_ylabel(f"P(any rejection) at FDR {a:g}")
            ax.set_title(f"{fam}: uniform null, {kernel.replace('_', ' ')}", fontsize=9)
            ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = plots_dir / f"calibration_vs_query_size.{ext}"
        fig.savefig(p, dpi=150)
        paths.append(p)
    plt.close(fig)
    return paths


def plot_background(out_dir: Path, plots_dir: Path, alpha: float = 0.05) -> list[Path]:
    plt = _mpl()
    paths: list[Path] = []
    frames = []
    for p in sorted((out_dir / "background").glob("*_queries.tsv")):
        df = _read_tsv(p)
        if df.empty:
            continue
        key = p.name[: -len("_queries.tsv")]
        org, rest = key.split("_", 1)
        fam, tier = rest.rsplit("_", 1)
        df["organism"], df["family"], df["tier"] = org, fam, tier
        frames.append(df)
    if not frames:
        return paths
    d = pd.concat(frames, ignore_index=True)
    col_bg = f"background_effect_overlap_jaccard_{alpha:g}"
    col_corr = f"correction_effect_annotation_jaccard_{alpha:g}"
    if col_bg not in d.columns:
        return paths
    families = sorted(d["family"].unique())
    fig, axes = plt.subplots(2, len(families), figsize=(6 * len(families), 7), squeeze=False)
    for c, fam in enumerate(families):
        sub = d[d["family"] == fam]
        groups = sorted(sub.groupby(["organism", "tier"]).groups.keys())
        for r, (col, title) in enumerate(
            (
                (col_bg, "background effect (annotation vs target-union, overlap-only BH)"),
                (
                    col_corr,
                    "correction effect (overlap-only vs fixed-family BH, annotation background)",
                ),
            )
        ):
            ax = axes[r][c]
            data = [
                sub[(sub["organism"] == o) & (sub["tier"] == t)][col].dropna().to_numpy()
                for o, t in groups
            ]
            ax.boxplot(data, showfliers=False)
            ax.set_xticks(range(1, len(groups) + 1))
            ax.set_xticklabels([f"{o}/{t}" for o, t in groups], rotation=90, fontsize=7)
            ax.set_ylim(-0.02, 1.02)
            ax.set_ylabel(f"Jaccard of significant sets (FDR {alpha:g})")
            ax.set_title(f"{fam}: {title}", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = plots_dir / f"background_sensitivity_jaccard.{ext}"
        fig.savefig(p, dpi=150)
        paths.append(p)
    plt.close(fig)
    return paths


def plot_degree(out_dir: Path, plots_dir: Path, alpha: float = 0.05) -> list[Path]:
    plt = _mpl()
    paths: list[Path] = []
    d = degree_table(out_dir)
    if d.empty:
        return paths
    col_d = f"degree_matched_any_rej_overlap_{alpha:g}"
    col_u = f"uniform_any_rej_overlap_{alpha:g}"
    if col_d not in d.columns:
        return paths
    col_p = f"pool_uniform_any_rej_overlap_{alpha:g}"
    panels = [(col_u, "uniform null (whole background)")]
    if col_p in d.columns:
        panels.append((col_p, "uniform null (target-union pool)"))
    fig, axes = plt.subplots(1, len(panels), figsize=(6 * len(panels), 6), squeeze=False)
    markers = {"GO_BP": "o", "REACTOME": "s", "GO_CC": "^", "GO_MF": "v", "WIKIPATHWAYS": "D"}
    for ax, (col_x, xlabel) in zip(axes[0], panels, strict=True):
        for fam in sorted(d["family"].unique()):
            sub = d[d["family"] == fam]
            marker = markers.get(fam, "o")
            fig_mask = sub["figure_example"].fillna("").astype(str) != ""
            ax.scatter(
                sub.loc[~fig_mask, col_x],
                sub.loc[~fig_mask, col_d],
                marker=marker,
                s=18,
                alpha=0.6,
                label=f"{fam} random focal",
            )
            ax.scatter(
                sub.loc[fig_mask, col_x],
                sub.loc[fig_mask, col_d],
                marker=marker,
                s=60,
                edgecolor="k",
                label=f"{fam} figure miRNA",
            )
        ax.plot([0, 1], [0, 1], color="grey", ls="--", lw=1)
        ax.set_xlabel(f"P(any rejection) under {xlabel} (FDR {alpha:g})")
        ax.set_ylabel(f"P(any rejection) under degree-matched null (FDR {alpha:g})")
        ax.set_title("Focal queries (overlap-only BH)", fontsize=9)
        ax.legend(fontsize=6)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = plots_dir / f"degree_null_any_rejection.{ext}"
        fig.savefig(p, dpi=150)
        paths.append(p)
    plt.close(fig)
    # Figure-term null distributions for the figure miRNAs.
    ft = figure_terms_table(out_dir)
    if not ft.empty and "k_null_degree_q95" in ft.columns:
        sub = ft[ft["figure_example"].fillna("").astype(str) != ""].copy()
        if not sub.empty:
            sub["label"] = sub["key"] + " " + sub["mirna"] + " " + sub["term"]
            sub = sub.sort_values("label")
            fig, ax = plt.subplots(figsize=(8, max(3, 0.22 * len(sub))))
            y = np.arange(len(sub))
            ax.barh(
                y,
                sub["k_null_degree_q95"],
                color="#f5b7b1",
                label="degree-matched null k (95th pct)",
            )
            ax.barh(
                y,
                sub["k_null_uniform_q95"],
                color="#aed6f1",
                alpha=0.6,
                label="uniform null k (95th pct)",
            )
            ax.scatter(sub["k_observed"], y, color="k", s=12, zorder=3, label="observed k")
            ax.set_yticks(y)
            ax.set_yticklabels(sub["label"], fontsize=5)
            ax.set_xlabel("overlap k")
            ax.legend(fontsize=7)
            ax.set_title(
                "Prespecified figure terms: observed overlap vs null 95th percentiles", fontsize=9
            )
            fig.tight_layout()
            for ext in ("png", "pdf"):
                p = plots_dir / f"figure_terms_null_overlap.{ext}"
                fig.savefig(p, dpi=150)
                paths.append(p)
            plt.close(fig)
    return paths


def write_plots(out_dir: Path) -> list[Path]:
    plots_dir = out_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    paths += plot_calibration(out_dir, plots_dir)
    paths += plot_calibration_vs_n(out_dir, plots_dir)
    paths += plot_background(out_dir, plots_dir)
    paths += plot_degree(out_dir, plots_dir)
    return paths
