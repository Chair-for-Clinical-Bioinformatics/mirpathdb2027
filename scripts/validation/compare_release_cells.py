"""Per-cell old-versus-new census of two release payloads (validation report only).

For every ``enrichment_significant/<org>/<FAMILY>_<TIER>.tsv`` present in either payload the
script reports retained rows, rows at the 0.05 and 0.005 thresholds, distinct miRNAs with at
least one retained / significant row, distinct pathways with at least one retained /
significant row, and -- from ``enrichment_inputs/<org>/<FAMILY>_term2gene.tsv`` -- the size
of the annotation inventory of the cell, so that a pathway losing all retained rows is
counted as "retained-empty" and never as unavailable. Counts are taken by header name from
every data row; nothing is inferred from file sizes or manifests.

Usage::

    PYTHONPATH=src python scripts/validation/compare_release_cells.py \\
        --old <prior payload> --new <candidate payload> --out <dir>

Writes ``cells.tsv`` (one row per cell), ``organism_totals.tsv`` and ``summary.json``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ALPHA = 0.05
STRINGENT = 0.005
TIERS = ("HIGH", "MEDIUM", "LOW")


def _cells(root: Path) -> dict[tuple[str, str, str], Path]:
    out: dict[tuple[str, str, str], Path] = {}
    for p in sorted((root / "enrichment_significant").glob("*/*.tsv")):
        fam, tier = p.stem.rsplit("_", 1)
        out[(p.parent.name, fam, tier)] = p
    return out


def _inventory(root: Path, org: str, fam: str) -> int | None:
    p = root / "enrichment_inputs" / org / f"{fam}_term2gene.tsv"
    if not p.is_file():
        return None
    return int(pd.read_csv(p, sep="\t", usecols=["pathway_id"]).pathway_id.nunique())


def _stats(path: Path | None) -> dict[str, object]:
    if path is None:
        return {"present": False}
    df = pd.read_csv(path, sep="\t")
    sig = df[df.p_adjust <= ALPHA]
    out: dict[str, object] = {
        "present": True,
        "rows": int(len(df)),
        "rows_significant": int(len(sig)),
        "rows_stringent": int((df.p_adjust <= STRINGENT).sum()),
        "mirnas_retained": int(df.mirna.nunique()),
        "mirnas_significant": int(sig.mirna.nunique()),
        "pathways_retained": int(df.pathway_id.nunique()),
        "pathways_significant": int(sig.pathway_id.nunique()),
        "columns": ",".join(df.columns),
    }
    if "p_adjust_overlap" in df.columns:
        out["rows_significant_overlap_family_within_retained"] = int(
            (df.p_adjust_overlap <= ALPHA).sum()
        )
        out["m_family"] = int(df.m_family.iloc[0]) if len(df) else None
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--old", type=Path, required=True)
    ap.add_argument("--new", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    old_cells, new_cells = _cells(args.old), _cells(args.new)
    keys = sorted(set(old_cells) | set(new_cells), key=lambda k: (k[0], k[1], TIERS.index(k[2])))
    rows = []
    for org, fam, tier in keys:
        o, n = _stats(old_cells.get((org, fam, tier))), _stats(new_cells.get((org, fam, tier)))
        inv_old, inv_new = _inventory(args.old, org, fam), _inventory(args.new, org, fam)
        rec: dict[str, object] = {"organism": org, "family": fam, "tier": tier}
        for k in (
            "rows",
            "rows_significant",
            "rows_stringent",
            "mirnas_retained",
            "mirnas_significant",
            "pathways_retained",
            "pathways_significant",
        ):
            rec[f"old_{k}"] = o.get(k)
            rec[f"new_{k}"] = n.get(k)
            if o.get(k) is not None and n.get(k) is not None:
                rec[f"delta_{k}"] = int(n[k]) - int(o[k])
        rec["inventory_pathways_old"] = inv_old
        rec["inventory_pathways_new"] = inv_new
        if inv_new is not None and n.get("pathways_retained") is not None:
            rec["new_pathways_retained_empty"] = inv_new - int(n["pathways_retained"])
        if inv_old is not None and o.get("pathways_retained") is not None:
            rec["old_pathways_retained_empty"] = inv_old - int(o["pathways_retained"])
        rec["new_m_family"] = n.get("m_family")
        rec["new_rows_significant_overlap_family_within_retained"] = n.get(
            "rows_significant_overlap_family_within_retained"
        )
        rows.append(rec)
    if not rows:
        raise SystemExit(f"no enrichment_significant cell found under {args.old} or {args.new}")
    cells = pd.DataFrame(rows)
    cells.to_csv(args.out / "cells.tsv", sep="\t", index=False)

    num = [c for c in cells.columns if c.startswith(("old_rows", "new_rows", "delta_rows"))]
    totals = cells.groupby("organism")[num].sum(min_count=1)
    totals.loc["all"] = totals.sum(min_count=1)
    totals.to_csv(args.out / "organism_totals.tsv", sep="\t")

    summary = {
        "old": str(args.old),
        "new": str(args.new),
        "cells": int(len(cells)),
        "cells_missing_in_new": [
            f"{r.organism}/{r.family}_{r.tier}" for r in cells.itertuples() if r.new_rows is None
        ],
        "totals": {
            k: (None if pd.isna(v) else int(v)) for k, v in totals.loc["all"].to_dict().items()
        },
        "inventory_unchanged": bool(
            (cells.inventory_pathways_old == cells.inventory_pathways_new).all()
        ),
        "new_cells_with_retained_empty_pathways": int(
            (cells.get("new_pathways_retained_empty", pd.Series(dtype=float)) > 0).sum()
        ),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
