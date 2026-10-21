#!/usr/bin/env python
"""Migrate a released payload to the contract-v2 full-family adjusted store.

For every (organism, family, tier) cell the script rebuilds every eligible
miRNA's query from the released ``target_pairs`` (production normalisation),
scores it with the independent kernel against the shipped ``enrichment_inputs``
and derives the full-family BH value ``q_full`` (BH over the ``M`` eligible terms
of the cell, ``M = space.n_terms``). Nothing is derived from the retained rows:
``m`` (the complete overlap-tested count of a query) and ``M`` both come from the
complete release inputs, never from the FDR-filtered store or from ``P < 1``.

Per cell, in this order:

1. ``M`` must equal the number of distinct ``pathway_id`` of the term2gene file
   (production cells: every listed term is eligible);
2. hard identity gate against the shipped
   ``enrichment_significant/<org>/<FAMILY>_<TIER>.tsv``: every shipped row is
   found in the reconstruction with the same overlap and with ``p_value`` /
   ``p_adjust`` within ``--tolerance`` (relative), and every reconstructed row
   with ``q_overlap <= 0.25`` (the store bound) exists in the store;
3. migration identity ``q_full == min(1, q_overlap * M / m)`` against the direct
   full-family BH for every miRNA (relative 1e-12);
4. only when 1-3 pass: the v2 cell file is written.

Output (``--out``)::

    enrichment_significant/<org>/<FAMILY>_<TIER>.tsv
        mirna, pathway_id, p_value, p_adjust, overlap,
        p_adjust_overlap, m_overlap, m_family
    pathway_labels.tsv      copied unchanged
    contract_version        "2" (stamped last, only when every cell passed)
    MIGRATION_REPORT.json   per-cell audit + totals
    COMPARE_REPORT.json     with --compare-with only

Row order of a v2 cell file: miRNAs in sorted order (the production writer's
``sorted(gene_sets)``), and within a miRNA the terms in term-index order, i.e.
sorted ``pathway_id`` (``FamilySpace.terms``). ``p_adjust`` is the direct
full-family BH value; ``p_adjust_overlap`` is the legacy overlap-only value the
shipped store carried. The first five columns keep the v1 name/order/type.

Usage::

    PYTHONPATH=src python scripts/migrate_full_family_adjust.py \\
        --payload /path/to/mirpathdb-data-2026-08-31 --out /path/to/v2 \\
        [--organisms hsa mmu] [--families GO_BP] [--tiers LOW] \\
        [--tolerance 1e-9] [--retention 0.25] [--compare-with /path/to/candidate]
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from mirpathdb_prod.enrichment import family_bh as fb
from mirpathdb_prod.reporting.significant import FDR_CUTOFF, FDR_SIGNIFICANT
from mirpathdb_prod.validation.robustness.inputs import TIERS, ReleaseInputs
from mirpathdb_prod.validation.robustness.kernel import ora_batch

CONTRACT_VERSION = "2"
STORE_BOUND = FDR_CUTOFF  # retention bound of the shipped v1 store
STRINGENT_ALPHA = 0.005
MIGRATION_RTOL = 1e-12
# --compare-with compares this script's scipy-based values with a candidate produced by the
# production R kernel; the two hypergeometric tails agree to ~6e-12 relative on the released
# cells (the fidelity audit's tolerance is the same 1e-9), so the cross-engine comparison
# cannot use the within-engine 1e-12 of the migration identity.
COMPARE_RTOL = 1e-9

V2_HEADER = (
    "mirna",
    "pathway_id",
    "p_value",
    "p_adjust",
    "overlap",
    "p_adjust_overlap",
    "m_overlap",
    "m_family",
)
V1_HEADER = ("mirna", "pathway_id", "p_value", "p_adjust", "overlap")

M_SOURCE = "reconstructed_from_release_inputs"
M_RULE = (
    "m (overlap-tested count per miRNA) = number of eligible terms with overlap >= 1 in the "
    "query rebuilt from target_pairs + enrichment_inputs, including terms whose raw P is "
    "exactly 1; M (family size) = eligible terms of the cell from the complete term2gene "
    "within the universe. m is never taken from retained rows of the shipped store and "
    "never from the count of rows with P < 1."
)

ROW_ORDER = (
    "mirnas sorted; within a miRNA, terms in term-index order (sorted pathway_id); "
    "rows with p_adjust (q_full) <= retention"
)


def _rel_diff(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    return np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)


def read_store(path: Path, header: tuple[str, ...] = V1_HEADER) -> pd.DataFrame:
    dtype = {
        "mirna": str,
        "pathway_id": str,
        "p_value": float,
        "p_adjust": float,
        "overlap": np.int64,
        "p_adjust_overlap": float,
        "m_overlap": np.int64,
        "m_family": np.int64,
    }
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={k: v for k, v in dtype.items() if k in header},
        float_precision="round_trip",
    )
    missing = [c for c in header if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    return df


def cell_path(root: Path, organism: str, family: str, tier: str) -> Path:
    return root / "enrichment_significant" / organism / f"{family}_{tier}.tsv"


def migrate_cell(
    inputs: ReleaseInputs,
    out_root: Path,
    organism: str,
    family: str,
    tier: str,
    *,
    tolerance: float = 1e-9,
    retention: float = FDR_CUTOFF,
    chunk: int = 256,
    log=print,
) -> dict[str, object]:
    """Audit one cell and, when every gate passes, write its v2 file.

    Returns the per-cell report record (``passed`` + counters). No output file
    exists for a cell whose record has ``passed == False``.
    """
    t0 = time.perf_counter()
    st = inputs.stratum(organism, family, tier)
    space = st.space
    store_path = cell_path(inputs.root, organism, family, tier)
    out_path = cell_path(out_root, organism, family, tier)
    m_family = int(space.n_terms)
    n_terms_in_file = int(space.provenance["n_terms_in_file"])
    rec: dict[str, object] = {
        "key": st.key,
        "organism": organism,
        "family": family,
        "tier": tier,
        "store_file": str(store_path),
        "output_file": str(out_path),
        "m_family": m_family,
        "n_terms_in_file": n_terms_in_file,
        "universe_equals_annotated_union": bool(
            space.provenance.get("universe_equals_annotated_union", True)
        ),
        "mirnas_eligible": len(st.queries),
        "tolerance_rel": tolerance,
        "migration_rtol": MIGRATION_RTOL,
        "retention": retention,
        "store_bound": STORE_BOUND,
    }
    failures: list[str] = []
    if m_family != n_terms_in_file:
        failures.append(
            f"family_size_mismatch: M={m_family} eligible terms but the term2gene file lists "
            f"{n_terms_in_file} distinct pathway_id (ineligible: {list(space.dropped_terms)[:5]})"
        )
    if not store_path.exists():
        failures.append(f"no_store_file: {store_path}")
    if failures:
        rec.update(status="aborted", passed=False, failures=failures)
        log(f"[migrate] {st.key}: ABORTED before scoring: {failures[0]}")
        return rec

    store = read_store(store_path)
    tindex = space.term_index()
    store["term_idx"] = store["pathway_id"].map(tindex)
    unknown_term = int(store["term_idx"].isna().sum())
    store_mirnas = set(store["mirna"].unique())
    store_not_eligible = sorted(store_mirnas - set(st.queries))
    groups = {m: g for m, g in store.groupby("mirna", sort=False)}
    names = sorted(st.queries)

    matched = overlap_mm = p_bad = padj_bad = 0
    max_rp = max_rq = 0.0
    missing_in_recon = 0
    recon_retained_absent = 0
    migration_max = 0.0
    migration_bad = 0
    rows_out = rows_sig = rows_str = 0
    examples: list[dict[str, object]] = []

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(".tsv.tmp")
    with tmp_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(V2_HEADER)
        for start in range(0, len(names), chunk):
            block = names[start : start + chunk]
            py = ora_batch(space, [st.queries[m] for m in block])
            for i, m in enumerate(block):
                hits = py.hits[i]
                m_overlap = int(hits.size)
                q_over = py.padj_overlap[i]
                q_full = py.padj_fixed[i]
                # ---- (a) identity gate against the shipped rows ----------------
                recon_retained = set(hits[q_over <= STORE_BOUND].tolist())
                g = groups.get(m)
                if g is None:
                    recon_retained_absent += len(recon_retained)
                else:
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
                    rp = _rel_diff(py.p[i][pos], sp)
                    rq = _rel_diff(q_over[pos], sq)
                    kmm = py.k[i][pos] != sk
                    matched += int(ok.sum())
                    overlap_mm += int(kmm.sum())
                    p_bad += int((rp > tolerance).sum())
                    padj_bad += int((rq > tolerance).sum())
                    if rp.size:
                        max_rp = max(max_rp, float(rp.max()))
                        max_rq = max(max_rq, float(rq.max()))
                    recon_retained_absent += len(recon_retained - set(idx[ok].tolist()))
                    if len(examples) < 10:
                        for b in np.flatnonzero((rp > tolerance) | (rq > tolerance) | kmm)[:3]:
                            examples.append(
                                {
                                    "mirna": m,
                                    "term": space.terms[hits[pos[b]]],
                                    "store_p": float(sp[b]),
                                    "recon_p": float(py.p[i][pos[b]]),
                                    "store_padj": float(sq[b]),
                                    "recon_padj": float(q_over[pos[b]]),
                                    "store_overlap": int(sk[b]),
                                    "recon_overlap": int(py.k[i][pos[b]]),
                                }
                            )
                # ---- (b) migration identity ---------------------------------------
                if m_overlap:
                    migrated = fb.migrate_overlap_adjusted(q_over, m_overlap, m_family)
                    rel = _rel_diff(migrated, q_full)
                    migration_max = max(migration_max, float(rel.max()))
                    migration_bad += int((rel > MIGRATION_RTOL).sum())
                # ---- (c) v2 rows ---------------------------------------------------
                keep = np.flatnonzero(q_full <= retention)
                for j in keep:
                    writer.writerow(
                        [
                            m,
                            space.terms[hits[j]],
                            repr(float(py.p[i][j])),
                            repr(float(q_full[j])),
                            int(py.k[i][j]),
                            repr(float(q_over[j])),
                            m_overlap,
                            m_family,
                        ]
                    )
                rows_out += int(keep.size)
                rows_sig += int((q_full <= FDR_SIGNIFICANT).sum())
                rows_str += int((q_full <= STRINGENT_ALPHA).sum())

    if unknown_term:
        failures.append(f"shipped_rows_unknown_term={unknown_term}")
    if store_not_eligible:
        failures.append(f"shipped_mirnas_not_eligible={len(store_not_eligible)}")
    if missing_in_recon:
        failures.append(f"shipped_rows_missing_in_reconstruction={missing_in_recon}")
    if overlap_mm:
        failures.append(f"overlap_mismatches={overlap_mm}")
    if p_bad:
        failures.append(f"rows_p_beyond_tolerance={p_bad}")
    if padj_bad:
        failures.append(f"rows_padj_beyond_tolerance={padj_bad}")
    if recon_retained_absent:
        failures.append(f"reconstructed_retained_rows_absent_from_store={recon_retained_absent}")
    if migration_bad:
        failures.append(f"migration_identity_violations={migration_bad}")
    passed = not failures
    if passed:
        tmp_path.replace(out_path)
    else:
        tmp_path.unlink(missing_ok=True)
    rec.update(
        {
            "status": "complete" if passed else "aborted",
            "passed": passed,
            "failures": failures,
            "shipped_rows": int(len(store)),
            "shipped_rows_unknown_term": unknown_term,
            "shipped_mirnas": len(store_mirnas),
            "shipped_mirnas_not_eligible": store_not_eligible[:20],
            "n_shipped_mirnas_not_eligible": len(store_not_eligible),
            "rows_matched": matched,
            "shipped_rows_missing_in_reconstruction": missing_in_recon,
            "overlap_mismatches": overlap_mm,
            "rows_p_beyond_tolerance": p_bad,
            "rows_padj_beyond_tolerance": padj_bad,
            "max_rel_diff_p": max_rp,
            "max_rel_diff_p_adjust": max_rq,
            "reconstructed_retained_rows_absent_from_store": recon_retained_absent,
            "migration_max_rel_diff": migration_max,
            "migration_identity_violations": migration_bad,
            "rows_out": rows_out if passed else 0,
            "rows_significant": rows_sig if passed else 0,
            "rows_stringent": rows_str if passed else 0,
            "examples": examples,
            "wall_s": round(time.perf_counter() - t0, 1),
        }
    )
    log(
        f"[migrate] {st.key}: M={m_family}, miRNAs {len(st.queries)}, shipped rows {len(store)}, "
        f"matched {matched}, overlap mismatches {overlap_mm}, p>tol {p_bad}, padj>tol {padj_bad}, "
        f"max rel p {max_rp:.2e} padj {max_rq:.2e}, absent-from-store {recon_retained_absent}, "
        f"migration max rel {migration_max:.2e}, rows_out {rec['rows_out']} "
        f"(sig {rec['rows_significant']}, stringent {rec['rows_stringent']}), "
        f"passed={passed} ({rec['wall_s']}s)"
    )
    return rec


def compare_cell(written: Path, candidate: Path) -> dict[str, object]:
    """Row-by-row comparison of a written v2 cell with the candidate's file."""
    rec: dict[str, object] = {"written": str(written), "candidate": str(candidate)}
    if not candidate.exists():
        rec.update(passed=False, failures=["candidate_file_missing"])
        return rec
    a = read_store(written, V2_HEADER)
    b = read_store(candidate, V2_HEADER)
    key = ["mirna", "pathway_id"]
    ka = list(zip(a["mirna"], a["pathway_id"], strict=True))
    kb = list(zip(b["mirna"], b["pathway_id"], strict=True))
    sa, sb = set(ka), set(kb)
    failures: list[str] = []
    if len(ka) != len(sa) or len(kb) != len(sb):
        failures.append("duplicate_keys")
    only_a = sorted(sa - sb)
    only_b = sorted(sb - sa)
    if only_a:
        failures.append(f"rows_only_in_written={len(only_a)}")
    if only_b:
        failures.append(f"rows_only_in_candidate={len(only_b)}")
    merged = a.merge(b, on=key, suffixes=("_w", "_c"))
    # Raw P is compared with the same relative tolerance as the adjusted values: the production
    # runner writes R's 15-significant-digit float strings, this script Python's shortest repr.
    counts = {
        "p_value_beyond_rtol": int(
            (_rel_diff(merged["p_value_w"], merged["p_value_c"]) > COMPARE_RTOL).sum()
        ),
        "overlap_differs": int((merged["overlap_w"] != merged["overlap_c"]).sum()),
        "p_adjust_beyond_rtol": int(
            (_rel_diff(merged["p_adjust_w"], merged["p_adjust_c"]) > COMPARE_RTOL).sum()
        ),
        "p_adjust_overlap_beyond_rtol": int(
            (
                _rel_diff(merged["p_adjust_overlap_w"], merged["p_adjust_overlap_c"]) > COMPARE_RTOL
            ).sum()
        ),
        "m_overlap_differs": int((merged["m_overlap_w"] != merged["m_overlap_c"]).sum()),
        "m_family_differs": int((merged["m_family_w"] != merged["m_family_c"]).sum()),
    }
    failures.extend(f"{k}={v}" for k, v in counts.items() if v)
    rec.update(
        {
            "rows_written": int(len(a)),
            "rows_candidate": int(len(b)),
            "rows_common": int(len(merged)),
            "rows_only_in_written": only_a[:10],
            "rows_only_in_candidate": only_b[:10],
            "row_order_identical": ka == kb,
            "max_rel_diff_p_value": (
                float(_rel_diff(merged["p_value_w"], merged["p_value_c"]).max())
                if len(merged)
                else 0.0
            ),
            "max_rel_diff_p_adjust": (
                float(_rel_diff(merged["p_adjust_w"], merged["p_adjust_c"]).max())
                if len(merged)
                else 0.0
            ),
            "max_rel_diff_p_adjust_overlap": (
                float(_rel_diff(merged["p_adjust_overlap_w"], merged["p_adjust_overlap_c"]).max())
                if len(merged)
                else 0.0
            ),
            **counts,
            "compare_rtol": COMPARE_RTOL,
            "passed": not failures,
            "failures": failures,
        }
    )
    return rec


def select_cells(
    inputs: ReleaseInputs,
    organisms: list[str] | None,
    families: list[str] | None,
    tiers: list[str] | None,
) -> list[tuple[str, str, str]]:
    cells = []
    for org in organisms or inputs.organisms():
        if org not in inputs.organisms():
            continue
        avail = inputs.available_families(org)
        for fam in families or avail:
            if fam not in avail:
                continue
            for tier in tiers or list(TIERS):
                cells.append((org, fam, tier))
    return cells


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--payload", type=Path, required=True, help="released payload directory")
    ap.add_argument("--out", type=Path, required=True, help="v2 output directory")
    ap.add_argument("--organisms", nargs="+", default=None)
    ap.add_argument("--families", nargs="+", default=None)
    ap.add_argument("--tiers", nargs="+", default=None)
    ap.add_argument("--tolerance", type=float, default=1e-9)
    ap.add_argument("--retention", type=float, default=FDR_CUTOFF)
    ap.add_argument("--compare-with", type=Path, default=None)
    ap.add_argument("--chunk", type=int, default=256)
    args = ap.parse_args(argv)

    payload = args.payload.resolve()
    out_root = args.out.resolve()
    if out_root == payload:
        raise SystemExit("--out must differ from --payload (inputs are never overwritten)")
    inputs = ReleaseInputs(payload, with_hashes=False)
    out_root.mkdir(parents=True, exist_ok=True)
    cells = select_cells(inputs, args.organisms, args.families, args.tiers)
    if not cells:
        raise SystemExit("no cell matches the requested organisms/families/tiers")

    labels_src = payload / "pathway_labels.tsv"
    labels_copied = labels_src.exists()
    if labels_copied:
        shutil.copyfile(labels_src, out_root / "pathway_labels.tsv")
    else:
        print(f"[migrate] WARNING: {labels_src} not found; labels not copied", flush=True)

    records: list[dict[str, object]] = []
    for org, fam, tier in cells:
        records.append(
            migrate_cell(
                inputs,
                out_root,
                org,
                fam,
                tier,
                tolerance=args.tolerance,
                retention=args.retention,
                chunk=args.chunk,
            )
        )
        inputs.drop_cached_stratum(org, fam, tier)

    all_passed = all(r["passed"] for r in records)
    stamp = out_root / "contract_version"
    if all_passed:
        stamp.write_text(CONTRACT_VERSION + "\n", encoding="utf-8")
    else:
        stamp.unlink(missing_ok=True)

    def total(key: str) -> int:
        return int(sum(int(r.get(key, 0)) for r in records))

    report = {
        "generator": "scripts/migrate_full_family_adjust.py",
        "payload": str(payload),
        "out": str(out_root),
        "release_info": inputs.release_info(),
        "contract_version": CONTRACT_VERSION if all_passed else None,
        "contract_version_stamped": all_passed,
        "pathway_labels_copied": labels_copied,
        "columns": list(V2_HEADER),
        "p_adjust_definition": "full-family BH over m_family eligible terms (direct computation)",
        "p_adjust_overlap_definition": "legacy overlap-only BH over m_overlap terms (v1 p_adjust)",
        "row_order": ROW_ORDER,
        "m_source": M_SOURCE,
        "m_rule": M_RULE,
        "retention": args.retention,
        "significant_alpha": FDR_SIGNIFICANT,
        "stringent_alpha": STRINGENT_ALPHA,
        "tolerance_rel": args.tolerance,
        "migration_rtol": MIGRATION_RTOL,
        "cells": records,
        "totals": {
            "cells": len(records),
            "cells_passed": sum(1 for r in records if r["passed"]),
            "cells_aborted": sum(1 for r in records if not r["passed"]),
            "mirnas_eligible": total("mirnas_eligible"),
            "shipped_rows": total("shipped_rows"),
            "rows_matched": total("rows_matched"),
            "overlap_mismatches": total("overlap_mismatches"),
            "rows_p_beyond_tolerance": total("rows_p_beyond_tolerance"),
            "rows_padj_beyond_tolerance": total("rows_padj_beyond_tolerance"),
            "reconstructed_retained_rows_absent_from_store": total(
                "reconstructed_retained_rows_absent_from_store"
            ),
            "max_rel_diff_p": max(
                (float(r.get("max_rel_diff_p", 0.0)) for r in records), default=0.0
            ),
            "max_rel_diff_p_adjust": max(
                (float(r.get("max_rel_diff_p_adjust", 0.0)) for r in records), default=0.0
            ),
            "migration_max_rel_diff": max(
                (float(r.get("migration_max_rel_diff", 0.0)) for r in records), default=0.0
            ),
            "rows_out": total("rows_out"),
            "rows_significant": total("rows_significant"),
            "rows_stringent": total("rows_stringent"),
        },
        "passed": all_passed,
    }
    (out_root / "MIGRATION_REPORT.json").write_text(
        json.dumps(report, indent=1) + "\n", encoding="utf-8"
    )
    print(
        f"[migrate] {report['totals']['cells_passed']}/{len(records)} cells passed, "
        f"rows_out {report['totals']['rows_out']} -> {out_root / 'MIGRATION_REPORT.json'}",
        flush=True,
    )
    rc = 0 if all_passed else 1

    if args.compare_with is not None:
        cand = args.compare_with.resolve()
        comparisons = []
        for r in records:
            if not r["passed"]:
                comparisons.append({"key": r["key"], "passed": False, "failures": ["cell_aborted"]})
                continue
            c = compare_cell(
                Path(str(r["output_file"])),
                cell_path(cand, str(r["organism"]), str(r["family"]), str(r["tier"])),
            )
            c["key"] = r["key"]
            comparisons.append(c)
            print(f"[compare] {r['key']}: passed={c['passed']} {c['failures']}", flush=True)
        cmp_passed = bool(comparisons) and all(c["passed"] for c in comparisons)
        (out_root / "COMPARE_REPORT.json").write_text(
            json.dumps(
                {
                    "written": str(out_root),
                    "candidate": str(cand),
                    "compare_rtol": COMPARE_RTOL,
                    "cells": comparisons,
                    "passed": cmp_passed,
                },
                indent=1,
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"[compare] passed={cmp_passed} -> {out_root / 'COMPARE_REPORT.json'}", flush=True)
        if not cmp_passed:
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
