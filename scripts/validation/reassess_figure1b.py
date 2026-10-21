#!/usr/bin/env python
"""Reassess a cross-species DicePlot figure (manuscript Figure 1B) under full-family BH.

For every MirGeneDB family-arm group x GO term x organism x evidence tier
("position"), resolve every mapped member miRNA of the group in that organism,
rebuild its query from the released inputs (production normalisation), score it
with the independent kernel (raw P, overlap-only BH = the shipped convention,
full-family BH = the candidate convention with M eligible terms), verify the
shipped row where one is retained, and derive the organism seat state under both
conventions from the smallest adjusted P among members with retained results.

States: ``stringent`` (adjusted P <= stringent), ``significant`` (<= alpha),
``none`` (member(s) tested, nothing at or below alpha), ``absent`` (no mapped member
of the group is eligible in that organism/tier -- no query exists), ``unavailable``
(term not in the organism's collection). ``none`` never means "tested negative"
for an individual member; the member table carries the per-member numbers.

Usage::

    PYTHONPATH=src python scripts/validation/reassess_figure1b.py \\
        --inputs-dir /path/to/mirpathdb-data-2026-08-31 --out /path/to/out \\
        --groups "MIR-124|3p" "MIR-9|5p" "MIR-132|3p" "MIR-137|3p" \\
        --terms GO:0030182 GO:0016358 GO:0007411 GO:0050808 \\
                GO:0007268 GO:0048167 GO:0043523 GO:0007611
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from mirpathdb_prod.validation.robustness.inputs import ReleaseInputs
from mirpathdb_prod.validation.robustness.kernel import ora_batch

ORGS = ("hsa", "mmu", "rno", "mml")


def load_members(mapping_path: Path, groups: list[str]) -> dict[str, dict[str, list[dict]]]:
    """``{group: {organism: [member rows]}}`` from the MirGeneDB mapping TSV."""
    out: dict[str, dict[str, list[dict]]] = {g: {o: [] for o in ORGS} for g in groups}
    with mapping_path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for r in reader:
            g = r.get("group_id")
            o = r.get("species")
            if g in out and o in ORGS:
                out[g][o].append(dict(r))
    return out


def state_from_q(q: float | None, alpha: float, stringent: float) -> str:
    if q is None or not np.isfinite(q):
        return "none"
    if q <= stringent:
        return "stringent"
    if q <= alpha:
        return "significant"
    return "none"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--inputs-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--groups", nargs="+", required=True)
    ap.add_argument("--terms", nargs="+", required=True)
    ap.add_argument("--family", default="GO_BP")
    ap.add_argument("--tiers", nargs="+", default=["HIGH", "LOW"])
    ap.add_argument("--organisms", nargs="+", default=list(ORGS))
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--stringent", type=float, default=0.005)
    ap.add_argument("--retention", type=float, default=0.25)
    ap.add_argument(
        "--mapping",
        type=Path,
        default=None,
        help="MirGeneDB mapping TSV (default: <inputs>/mirgenedb_v3_miRNA_mapping.tsv)",
    )
    ap.add_argument(
        "--in-atlas-only", action="store_true", help="count only in_atlas == yes members as mapped"
    )
    args = ap.parse_args(argv)

    inputs = ReleaseInputs(args.inputs_dir, with_hashes=False)
    mapping = args.mapping or (inputs.root / "mirgenedb_v3_miRNA_mapping.tsv")
    members = load_members(mapping, args.groups)
    args.out.mkdir(parents=True, exist_ok=True)

    member_rows: list[dict] = []
    positions: list[dict] = []
    for org in args.organisms:
        for tier in args.tiers:
            if args.family not in inputs.available_families(org):
                for g in args.groups:
                    for t in args.terms:
                        positions.append(
                            _position(
                                g,
                                t,
                                org,
                                tier,
                                "unavailable",
                                "unavailable",
                                [],
                                [],
                                "family not shipped",
                            )
                        )
                continue
            st = inputs.stratum(org, args.family, tier)
            space = st.space
            tindex = space.term_index()
            m_full = space.n_terms
            store = _load_store(
                inputs.root / "enrichment_significant" / org / f"{args.family}_{tier}.tsv"
            )
            # Resolve members and score the eligible ones in one batch.
            group_members: dict[str, list[dict]] = {}
            names: list[str] = []
            for g in args.groups:
                rows = members[g][org]
                if args.in_atlas_only:
                    rows = [r for r in rows if r.get("in_atlas") == "yes"]
                seen: set[str] = set()
                uniq = []
                for r in rows:
                    nm = r["mirna_name"]
                    if nm in seen:
                        continue
                    seen.add(nm)
                    uniq.append(r)
                group_members[g] = uniq
                names.extend(r["mirna_name"] for r in uniq)
            eligible_names = [nm for nm in dict.fromkeys(names) if nm in st.queries]
            batch = (
                ora_batch(space, [st.queries[nm] for nm in eligible_names])
                if eligible_names
                else None
            )
            per_member: dict[str, dict] = {}
            for i, nm in enumerate(eligible_names):
                hits = batch.hits[i]
                d = {
                    "n": int(batch.n[i]),
                    "m_overlap": int(hits.size),
                    "terms": {},
                }
                for t in args.terms:
                    j = tindex.get(t)
                    if j is None:
                        d["terms"][t] = None
                        continue
                    pos = np.searchsorted(hits, j)
                    if pos < hits.size and hits[pos] == j:
                        d["terms"][t] = {
                            "k": int(batch.k[i][pos]),
                            "K": int(batch.K[i][pos]),
                            "p": float(batch.p[i][pos]),
                            "q_overlap": float(batch.padj_overlap[i][pos]),
                            "q_full": float(batch.padj_fixed[i][pos]),
                        }
                    else:
                        d["terms"][t] = {
                            "k": 0,
                            "K": int(space.term_sizes[j]),
                            "p": 1.0,
                            "q_overlap": None,
                            "q_full": 1.0,
                        }
                per_member[nm] = d
            for g in args.groups:
                for t in args.terms:
                    term_present = t in tindex
                    orig_candidates: list[tuple[float, str]] = []
                    corr_candidates: list[tuple[float, str]] = []
                    n_eligible = 0
                    member_details = []
                    for r in group_members[g]:
                        nm = r["mirna_name"]
                        info = per_member.get(nm)
                        raw_n = st.n_raw.get(nm, 0)
                        shipped = store.get((nm, t))
                        row = {
                            "group": g,
                            "term": t,
                            "organism": org,
                            "tier": tier,
                            "member": nm,
                            "in_atlas": r.get("in_atlas"),
                            "mimat": r.get("mimat") or r.get("accession") or "",
                            "member_has_targets_at_tier": raw_n > 0,
                            "n_raw_cumulative": raw_n,
                            "member_eligible": info is not None,
                            "n": info["n"] if info else 0,
                            "m_overlap": info["m_overlap"] if info else 0,
                            "M_full": m_full,
                            "N": space.N,
                            "term_in_collection": term_present,
                        }
                        if info and term_present and info["terms"].get(t):
                            tt = info["terms"][t]
                            row.update(
                                {
                                    "k": tt["k"],
                                    "K": tt["K"],
                                    "p_reconstructed": tt["p"],
                                    "q_overlap_reconstructed": tt["q_overlap"],
                                    "q_full": tt["q_full"],
                                    "shipped_p_value": shipped["p_value"] if shipped else None,
                                    "shipped_p_adjust": shipped["p_adjust"] if shipped else None,
                                    "shipped_overlap": shipped["overlap"] if shipped else None,
                                    "shipped_row_retained": shipped is not None,
                                    "reconstruction_matches_shipped": (
                                        abs(shipped["p_adjust"] - tt["q_overlap"])
                                        <= 1e-9 * max(shipped["p_adjust"], 1e-300)
                                        and shipped["overlap"] == tt["k"]
                                    )
                                    if (shipped and tt["q_overlap"] is not None)
                                    else None,
                                    "state_original_member": state_from_q(
                                        shipped["p_adjust"], args.alpha, args.stringent
                                    )
                                    if shipped
                                    else (
                                        "none"
                                        if tt["q_overlap"] is not None
                                        and tt["q_overlap"] > args.retention
                                        else ("none" if tt["k"] == 0 else "none")
                                    ),
                                    "state_corrected_member": state_from_q(
                                        tt["q_full"], args.alpha, args.stringent
                                    ),
                                    "retained_original": shipped is not None,
                                    "retained_corrected": tt["q_full"] <= args.retention,
                                }
                            )
                            n_eligible += 1
                            if shipped:
                                orig_candidates.append((shipped["p_adjust"], nm))
                            if tt["q_full"] <= args.retention:
                                corr_candidates.append((tt["q_full"], nm))
                        elif info and not term_present:
                            row.update({"note": "term not in organism collection"})
                            n_eligible += 1
                        member_rows.append(row)
                        member_details.append(row)
                    if not term_present:
                        s_orig = s_corr = "unavailable"
                        reason = "term absent from the organism's collection"
                    elif n_eligible == 0:
                        s_orig = s_corr = "absent"
                        reason = (
                            "no mapped member of the group has a query at this tier (not eligible)"
                        )
                    else:
                        best_o = min(orig_candidates) if orig_candidates else None
                        best_c = min(corr_candidates) if corr_candidates else None
                        s_orig = (
                            state_from_q(best_o[0], args.alpha, args.stringent)
                            if best_o
                            else "none"
                        )
                        s_corr = (
                            state_from_q(best_c[0], args.alpha, args.stringent)
                            if best_c
                            else "none"
                        )
                        reason = (
                            f"{n_eligible} eligible member(s); original best {best_o[1]}"
                            f" q={best_o[0]:.3g}"
                            if best_o
                            else (
                                f"{n_eligible} eligible member(s); "
                                f"no retained row at {args.retention}"
                            )
                        )
                        reason += (
                            f"; corrected best {best_c[1]} q_full={best_c[0]:.3g}"
                            if best_c
                            else f"; corrected: no q_full <= {args.retention}"
                        )
                    positions.append(
                        _position(
                            g,
                            t,
                            org,
                            tier,
                            s_orig,
                            s_corr,
                            orig_candidates,
                            corr_candidates,
                            reason,
                            n_eligible=n_eligible,
                            n_members=len(group_members[g]),
                            m_full=m_full,
                        )
                    )
            inputs.drop_cached_stratum(org, args.family, tier)

    pos_df = pd.DataFrame(positions)
    mem_df = pd.DataFrame(member_rows)
    pos_df.to_csv(args.out / "figure1b_positions.tsv", sep="\t", index=False)
    mem_df.to_csv(args.out / "figure1b_members.tsv", sep="\t", index=False)
    changes = pos_df[pos_df.state_original != pos_df.state_corrected]
    summary = {
        "positions": int(len(pos_df)),
        "changed_positions": int(len(changes)),
        "state_counts_original": pos_df.state_original.value_counts().to_dict(),
        "state_counts_corrected": pos_df.state_corrected.value_counts().to_dict(),
        "transitions": changes.groupby(["state_original", "state_corrected"])
        .size()
        .reset_index(name="n")
        .to_dict(orient="records"),
        "members_reconstruction_mismatches": (
            int((~mem_df["reconstruction_matches_shipped"].astype(bool)).sum())
            if "reconstruction_matches_shipped" in mem_df
            else 0
        ),
    }
    (args.out / "figure1b_summary.json").write_text(
        json.dumps(summary, indent=1, default=str), encoding="utf-8"
    )
    print(json.dumps(summary, indent=1, default=str))
    return 0


def _position(
    g,
    t,
    org,
    tier,
    s_orig,
    s_corr,
    orig_c,
    corr_c,
    reason,
    *,
    n_eligible=0,
    n_members=0,
    m_full=None,
):
    best_o = min(orig_c) if orig_c else None
    best_c = min(corr_c) if corr_c else None
    return {
        "group": g,
        "term": t,
        "organism": org,
        "tier": tier,
        "state_original": s_orig,
        "state_corrected": s_corr,
        "original_best_member": best_o[1] if best_o else "",
        "original_best_q": best_o[0] if best_o else None,
        "corrected_best_member": best_c[1] if best_c else "",
        "corrected_best_q_full": best_c[0] if best_c else None,
        "n_members_mapped": n_members,
        "n_members_eligible": n_eligible,
        "M_full": m_full,
        "reason": reason,
    }


def _load_store(path: Path) -> dict[tuple[str, str], dict]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype={"mirna": str, "pathway_id": str})
    return {
        (r.mirna, r.pathway_id): {
            "p_value": float(r.p_value),
            "p_adjust": float(r.p_adjust),
            "overlap": int(r.overlap),
        }
        for r in df.itertuples(index=False)
    }


if __name__ == "__main__":
    sys.exit(main())
