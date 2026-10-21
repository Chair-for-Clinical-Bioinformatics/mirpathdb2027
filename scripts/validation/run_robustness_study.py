#!/usr/bin/env python
"""Single entry point of the statistical robustness study (validation only).

Stages (``--stage``; ``all`` runs them in order, every stage is resumable):

  scope        enumerate strata, shared cohorts, actual query sizes; write manifest
  calibration  uniform-null calibration: overlap-only BH vs fixed-family BH
               (independent kernel + unchanged production kernel on the same samples)
  background   target-union background vs annotation background on real queries (2x2)
  degree       degree-binned (target-frequency-aware) null for focal queries
  report       summary tables + plots

Usage (from a checkout, with the ``validation`` extra installed)::

    PYTHONPATH=src python scripts/validation/run_robustness_study.py \\
        --config config/robustness_study.yaml \\
        --inputs-dir /path/to/mirpathdb-data-2026-08-31 \\
        --out /path/to/results --stage all

``--dry-run`` writes the scope + manifest and prints the planned work without
computing anything. ``--r-kernel off`` disables the production-kernel runs
(results are then labelled "independent implementation" only).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from mirpathdb_prod.validation.robustness.background import run_background_stratum
from mirpathdb_prod.validation.robustness.calibration import run_calibration_config
from mirpathdb_prod.validation.robustness.degree import run_degree_stratum
from mirpathdb_prod.validation.robustness.inputs import ReleaseInputs
from mirpathdb_prod.validation.robustness.report import write_plots, write_summaries
from mirpathdb_prod.validation.robustness.study import (
    StudyConfig,
    build_scope,
    update_manifest,
    write_manifest,
    write_scope,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
STAGES = ("scope", "calibration", "background", "degree", "report", "all")


def _log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--config", type=Path, default=REPO_ROOT / "config" / "robustness_study.yaml")
    ap.add_argument(
        "--inputs-dir",
        type=Path,
        required=True,
        help="release payload dir (enrichment_inputs/, target_pairs/)",
    )
    ap.add_argument("--out", type=Path, required=True, help="output directory (created; resumable)")
    ap.add_argument("--stage", choices=STAGES, default="all")
    ap.add_argument("--seed", type=int, default=None, help="override the study seed")
    ap.add_argument(
        "--replicates", type=int, default=None, help="override replicate count (all stages)"
    )
    ap.add_argument(
        "--strata", nargs="*", default=None, help="restrict to stratum keys <org>_<FAMILY>_<TIER>"
    )
    ap.add_argument(
        "--no-grid",
        action="store_true",
        help="skip the extra_query_sizes grid (percentile configurations only)",
    )
    ap.add_argument("--r-kernel", choices=("config", "on", "off"), default="config")
    ap.add_argument("--rscript", default=None, help="Rscript of the pinned clusterProfiler env")
    ap.add_argument("--max-cores", type=int, default=None, help="MPD_MAX_CORES for the R driver")
    ap.add_argument("--provenance-dir", type=Path, default=REPO_ROOT / "provenance")
    ap.add_argument("--dry-run", action="store_true", help="scope + manifest only; print the plan")
    ap.add_argument("--force", action="store_true", help="recompute even if outputs exist")
    ap.add_argument(
        "--force-python",
        action="store_true",
        help="calibration: recompute the independent kernel + agreement, keep cached R outputs",
    )
    ap.add_argument("--no-hashes", action="store_true", help="skip SHA-256 of inputs (tests)")
    args = ap.parse_args(argv)

    cfg = StudyConfig.from_yaml(args.config)
    if args.seed is not None:
        cfg.seed = int(args.seed)
        cfg.raw["seed"] = cfg.seed
    if args.replicates is not None:
        cfg.replicates = int(args.replicates)
        cfg.raw["replicates"] = cfg.replicates
        cfg.degree_null["replicates"] = cfg.replicates
    if args.r_kernel == "on":
        cfg.r_kernel["enabled"] = True
    elif args.r_kernel == "off":
        cfg.r_kernel["enabled"] = False
    if args.rscript:
        cfg.r_kernel["rscript"] = args.rscript
    if args.max_cores:
        cfg.r_kernel["max_cores"] = int(args.max_cores)

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    inputs = ReleaseInputs(args.inputs_dir, with_hashes=not args.no_hashes)
    limitations: list[str] = [str(x) for x in (cfg.raw.get("limitations") or [])]
    if not cfg.r_kernel.get("enabled"):
        limitations.append(
            "production R kernel disabled for this run: numbers come from the independent "
            "numpy/scipy implementation and are not production-kernel runs"
        )

    _log(f"scope: inputs={inputs.root} release={inputs.release_info().get('release_id', '?')}")
    strata, configs, cohorts = build_scope(
        inputs,
        cfg,
        provenance_dir=args.provenance_dir,
        strata_filter=args.strata,
        include_grid=not args.no_grid,
        log=_log,
    )
    for s in strata:
        if s.n_manifest is not None and s.n_eligible > s.n_manifest:
            limitations.append(
                f"{s.key}: reconstructed cohort {s.n_eligible} > manifest {s.n_manifest}"
            )
        if s.tier != "*" and not s.available:
            limitations.append(f"{s.key}: empty stratum (no eligible miRNA)")
    if args.strata or args.no_grid:
        # A restricted run must never overwrite the canonical scope/manifest of the
        # full study; it only reports what it is about to run.
        _log(
            f"scope: restricted run ({len(configs)} configurations) -- "
            "canonical scope/manifest left untouched"
        )
    else:
        write_scope(out, strata, configs, cohorts, inputs, cfg)
        write_manifest(
            out,
            cfg=cfg,
            inputs=inputs,
            strata=strata,
            configs=configs,
            repo_root=REPO_ROOT,
            limitations=limitations,
        )
    avail = [s for s in strata if s.available and s.tier != "*"]
    _log(
        f"scope: {len(avail)} strata, {len(configs)} calibration configurations, "
        f"replicates={cfg.replicates}"
    )
    if args.dry_run or args.stage == "scope":
        for c in configs:
            _log(f"  plan calibration {c.key}")
        return 0

    stage = args.stage
    figure_terms = {
        (s.organism, s.family): cfg.figure_terms_for(s.organism, s.family) for s in avail
    }

    if stage in ("calibration", "all"):
        t0 = time.perf_counter()
        r_prov: dict[str, object] = {}
        for c in configs:
            st = inputs.stratum(c.organism, c.family, c.tier)
            summary = run_calibration_config(
                st,
                size_label=c.label,
                n=c.n,
                reps=cfg.replicates,
                seed=cfg.seed,
                alphas=cfg.alphas,
                retention_alpha=cfg.retention_alpha,
                out_dir=out / "calibration",
                r_kernel=cfg.r_kernel,
                force=args.force,
                force_python=args.force_python,
                log=_log,
            )
            pk = summary.get("production_kernel", {})
            if pk.get("provenance") and not r_prov:
                r_prov = {
                    k: v
                    for k, v in pk["provenance"].items()
                    if k in ("clusterprofiler_version", "r_version", "detected_cores")
                }
        update_manifest(
            out,
            production_kernel_provenance=r_prov,
            calibration_wall_s=round(time.perf_counter() - t0, 1),
        )
        _log(f"calibration done in {time.perf_counter() - t0:.0f}s")

    if stage in ("background", "all"):
        t0 = time.perf_counter()
        for s in avail:
            st = inputs.stratum(s.organism, s.family, s.tier)
            run_background_stratum(
                st,
                alphas=cfg.alphas,
                out_dir=out / "background",
                figure_mirnas=cfg.figure_mirnas(s.organism),
                figure_terms=figure_terms[(s.organism, s.family)],
                cohort=cohorts.get(f"{s.organism}_{s.family}"),
                r_kernel=cfg.r_kernel,
                force=args.force,
                log=_log,
            )
        update_manifest(out, background_wall_s=round(time.perf_counter() - t0, 1))
        _log(f"background done in {time.perf_counter() - t0:.0f}s")

    if stage in ("degree", "all"):
        t0 = time.perf_counter()
        dn = cfg.degree_null
        for s in avail:
            st = inputs.stratum(s.organism, s.family, s.tier)
            run_degree_stratum(
                st,
                seed=cfg.seed,
                alphas=cfg.alphas,
                out_dir=out / "degree",
                figure_mirnas=cfg.figure_mirnas(s.organism),
                figure_terms=figure_terms[(s.organism, s.family)],
                n_bins=int(dn.get("n_bins", 10)),
                reps=int(dn.get("replicates", cfg.replicates)),
                random_focal=int(dn.get("random_focal_per_stratum", 20)),
                force=args.force,
                log=_log,
            )
        update_manifest(out, degree_wall_s=round(time.perf_counter() - t0, 1))
        _log(f"degree done in {time.perf_counter() - t0:.0f}s")

    if stage in ("report", "all"):
        written = write_summaries(out)
        plots = write_plots(out)
        _log(
            "report: " + ", ".join(str(p.relative_to(out)) for p in list(written.values()) + plots)
        )
        update_manifest(
            out,
            report={
                "summaries": [str(p.relative_to(out)) for p in written.values()],
                "plots": [str(p.relative_to(out)) for p in plots],
            },
        )
    _log("done")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print(json.dumps({"status": "interrupted"}), file=sys.stderr)
        sys.exit(130)
