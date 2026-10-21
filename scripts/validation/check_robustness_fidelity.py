#!/usr/bin/env python
"""Release-fidelity audit of the robustness study's reconstruction chain.

For every requested stratum, rebuild each eligible miRNA's query from the
released ``target_pairs`` (production normalisation), score it with the
independent kernel against the shipped ``enrichment_inputs``, and compare with
EVERY row of the shipped ``enrichment_significant/<org>/<FAMILY>_<TIER>.tsv``
(the retained rows, BH-FDR <= 0.25): raw P, adjusted P and overlap must agree.

This is a consistency audit of the study's inputs and kernel against the
release, not a reconstruction of the study from the retained rows (queries and
BH families are always taken from the complete inputs). It also reports
reconstructed rows at or below the storage bound that are absent from the
store, so a silently truncated store would be visible.

Usage::

    PYTHONPATH=src python scripts/validation/check_robustness_fidelity.py \\
        --inputs-dir /path/to/mirpathdb-data-2026-08-31 --out /path/to/results \\
        --organisms hsa mmu rno mml --families GO_BP REACTOME --tiers HIGH MEDIUM LOW
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from mirpathdb_prod.validation.robustness.inputs import ReleaseInputs
from mirpathdb_prod.validation.robustness.kernel import ora_batch

STORE_BOUND = 0.25


def audit_stratum(
    inputs: ReleaseInputs,
    organism: str,
    family: str,
    tier: str,
    *,
    tolerance: float = 1e-9,
    chunk: int = 256,
    log=print,
) -> dict[str, object]:
    st = inputs.stratum(organism, family, tier)
    space = st.space
    store_path = inputs.root / "enrichment_significant" / organism / f"{family}_{tier}.tsv"
    if not store_path.exists():
        return {"key": st.key, "status": "no_store_file", "path": str(store_path)}
    t0 = time.perf_counter()
    store = pd.read_csv(
        store_path,
        sep="\t",
        dtype={
            "mirna": str,
            "pathway_id": str,
            "p_value": float,
            "p_adjust": float,
            "overlap": int,
        },
    )
    tindex = space.term_index()
    store["term_idx"] = store["pathway_id"].map(tindex)
    n_store = int(len(store))
    unknown_term = int(store["term_idx"].isna().sum())
    store_mirnas = set(store["mirna"].unique())
    eligible = set(st.queries)
    store_not_eligible = sorted(store_mirnas - eligible)
    groups = {m: g for m, g in store.groupby("mirna", sort=False)}
    names = sorted(eligible)
    matched = 0
    overlap_mm = 0
    p_bad = 0
    padj_bad = 0
    max_rp = 0.0
    max_rq = 0.0
    missing_in_recon = 0  # store rows whose term the reconstruction did not score
    recon_retained_not_in_store = 0
    examples: list[dict[str, object]] = []
    for start in range(0, len(names), chunk):
        block = names[start : start + chunk]
        py = ora_batch(space, [st.queries[m] for m in block])
        for i, m in enumerate(block):
            hits = py.hits[i]
            g = groups.get(m)
            retained_mask = py.padj_overlap[i] <= STORE_BOUND
            recon_retained = set(hits[retained_mask].tolist())
            if g is None:
                recon_retained_not_in_store += len(recon_retained)
                continue
            g = g[g["term_idx"].notna()]
            idx = g["term_idx"].astype(np.int64).to_numpy()
            pos = np.searchsorted(hits, idx)
            ok = (
                (pos < hits.size) & (hits[np.minimum(pos, max(hits.size - 1, 0))] == idx)
                if hits.size
                else np.zeros(len(idx), bool)
            )
            missing_in_recon += int((~ok).sum())
            pos = pos[ok]
            sp = g["p_value"].to_numpy()[ok]
            sq = g["p_adjust"].to_numpy()[ok]
            sk = g["overlap"].to_numpy()[ok]
            rp = np.abs(py.p[i][pos] - sp) / np.maximum(
                np.maximum(np.abs(sp), np.abs(py.p[i][pos])), 1e-300
            )
            rq = np.abs(py.padj_overlap[i][pos] - sq) / np.maximum(
                np.maximum(np.abs(sq), np.abs(py.padj_overlap[i][pos])), 1e-300
            )
            matched += int(ok.sum())
            overlap_mm += int((py.k[i][pos] != sk).sum())
            p_bad += int((rp > tolerance).sum())
            padj_bad += int((rq > tolerance).sum())
            if rp.size:
                max_rp = max(max_rp, float(rp.max()))
                max_rq = max(max_rq, float(rq.max()))
            store_terms = set(idx[ok].tolist())
            recon_retained_not_in_store += len(recon_retained - store_terms)
            if len(examples) < 10:
                bad = np.flatnonzero((rp > tolerance) | (rq > tolerance) | (py.k[i][pos] != sk))
                for b in bad[:3]:
                    examples.append(
                        {
                            "mirna": m,
                            "term": space.terms[hits[pos[b]]],
                            "store_p": float(sp[b]),
                            "recon_p": float(py.p[i][pos[b]]),
                            "store_padj": float(sq[b]),
                            "recon_padj": float(py.padj_overlap[i][pos[b]]),
                            "store_overlap": int(sk[b]),
                            "recon_overlap": int(py.k[i][pos[b]]),
                        }
                    )
    out = {
        "key": st.key,
        "status": "complete",
        "store_file": str(store_path),
        "store_rows": n_store,
        "store_rows_unknown_term": unknown_term,
        "store_mirnas": len(store_mirnas),
        "store_mirnas_not_eligible": store_not_eligible[:20],
        "n_store_mirnas_not_eligible": len(store_not_eligible),
        "eligible_mirnas": len(eligible),
        "rows_matched": matched,
        "rows_missing_in_reconstruction": missing_in_recon,
        "overlap_mismatches": overlap_mm,
        "rows_p_beyond_tolerance": p_bad,
        "rows_padj_beyond_tolerance": padj_bad,
        "max_rel_diff_p": max_rp,
        "max_rel_diff_padj": max_rq,
        "reconstructed_retained_rows_absent_from_store": recon_retained_not_in_store,
        "tolerance_rel": tolerance,
        "passed": (
            unknown_term == 0
            and not store_not_eligible
            and missing_in_recon == 0
            and overlap_mm == 0
            and p_bad == 0
            and padj_bad == 0
            and recon_retained_not_in_store == 0
        ),
        "examples": examples,
        "wall_s": round(time.perf_counter() - t0, 1),
    }
    log(
        f"[fidelity] {st.key}: store rows {n_store}, matched {matched}, "
        f"overlap mismatches {overlap_mm}, p>tol {p_bad}, padj>tol {padj_bad}, "
        f"max rel p {max_rp:.2e} padj {max_rq:.2e}, "
        f"absent-from-store {recon_retained_not_in_store}, passed={out['passed']} "
        f"({out['wall_s']}s)"
    )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--inputs-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--organisms", nargs="+", default=["hsa", "mmu", "rno", "mml"])
    ap.add_argument("--families", nargs="+", default=["GO_BP", "REACTOME"])
    ap.add_argument("--tiers", nargs="+", default=["HIGH", "MEDIUM", "LOW"])
    ap.add_argument("--tolerance", type=float, default=1e-9)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    inputs = ReleaseInputs(args.inputs_dir, with_hashes=False)
    out_dir = args.out / "fidelity"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for org in args.organisms:
        if org not in inputs.organisms():
            continue
        for fam in args.families:
            if fam not in inputs.available_families(org):
                continue
            for tier in args.tiers:
                path = out_dir / f"{org}_{fam}_{tier}.json"
                if path.exists() and not args.force:
                    rec = json.loads(path.read_text(encoding="utf-8"))
                else:
                    rec = audit_stratum(inputs, org, fam, tier, tolerance=args.tolerance)
                    path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
                rows.append(
                    {
                        k: v
                        for k, v in rec.items()
                        if k not in ("examples", "store_mirnas_not_eligible")
                    }
                )
                inputs.drop_cached_stratum(org, fam, tier)
    pd.DataFrame(rows).to_csv(out_dir / "summary.tsv", sep="\t", index=False)
    print(f"[fidelity] summary -> {out_dir / 'summary.tsv'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
