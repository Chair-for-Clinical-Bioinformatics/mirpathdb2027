#!/usr/bin/env python
"""FULL verification of an exported target-pairs artefact set.

In-repo port of the external one-off ``tp_verify_20260816.py`` that re-verified
the ``target-pairs-2026-08-16`` package against the frozen contract (2026-08-24
audit, finding (d)). For EVERY row of EVERY frozen contract cell (57 files,
~19.7M rows for the full four-organism release) recompute

    overlap = len(term2gene[pathway_id] & universe & cumulative_targets(mirna, tier))

using ONLY the exported artefacts (target_pairs TSVs + term2gene/universe TSVs)
plus this repo's own normalize/cumulative functions (the frozen formula in
``export_frontend_contracts.export_significant`` -- set intersection is
commutative, so the member-first order used here for speed computes the
identical number), and compare to the frozen ``overlap`` column. Zero mismatches
required; a non-zero total exits 1.

Writes incremental JSON to ``--result`` after every cell, so a long run can be
monitored and resumed-from-inspection if interrupted.

``--contract-dir`` accepts EITHER the release root (the directory that contains
``enrichment_significant/``, e.g. ``data/processed`` or an assembled payload)
OR the ``enrichment_significant`` directory itself: when the given directory has
an ``enrichment_significant/`` sub-directory, the cells are read from there.
(Before 2026-09-15 only the second form worked; passing the release root failed
with ``FileNotFoundError`` on ``<root>/hsa/GO_BP_HIGH.tsv``.)

Usage (from the repo root)::

    PYTHONPATH=src python scripts/verify_target_pairs.py \\
      --pairs-dir /path/to/target-pairs-<date> \\
      --contract-dir /path/to/release \\
      --result /tmp/tp_verify_results.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPTS_DIR.parent

# src/ layout bootstrap: allow `python scripts/verify_target_pairs.py` without an
# editable install (mirrors the documented `PYTHONPATH=src` convention).
if str(_REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "src"))

from mirpathdb_prod.sources import production_target_sets as pts  # noqa: E402

FAMILIES_BY_ORG = {
    "hsa": ("GO_BP", "GO_CC", "GO_MF", "REACTOME", "WIKIPATHWAYS"),
    "mmu": ("GO_BP", "GO_CC", "GO_MF", "REACTOME", "WIKIPATHWAYS"),
    "rno": ("GO_BP", "GO_CC", "GO_MF", "REACTOME", "WIKIPATHWAYS"),
    "mml": ("GO_BP", "GO_CC", "GO_MF", "REACTOME"),
}
TIERS = ("HIGH", "MEDIUM", "LOW")
EV_TO_BAND = {"High": "HIGH", "Medium": "MEDIUM", "Low": "LOW"}

_T0 = time.time()


def log(*a: object) -> None:
    print(f"[{time.time() - _T0:8.1f}s]", *a, flush=True)


def load_pairs(pairs_dir: Path, org: str, results: dict) -> dict[str, dict[str, set[str]]]:
    """Rebuild {band: {mirna: set(gene)}} from the exported atomic pair TSV.

    Atomic bands are disjoint by construction, so band maps here ARE the atomic
    bands and pts.cumulative_tier() rebuilds the cumulative H / H|M / H|M|L sets.
    Also asserts: unique (mirna, gene) rows, valid evidence labels.
    """
    tiers: dict[str, dict[str, set[str]]] = {"HIGH": {}, "MEDIUM": {}, "LOW": {}}
    seen: set[tuple[str, str]] = set()
    n_rows = 0
    dup = 0
    path = pairs_dir / "target_pairs" / f"{org}_target_pairs.tsv"
    with path.open("r", encoding="utf-8", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        assert rdr.fieldnames == [
            "organism", "mirna", "target_gene", "atomic_evidence", "source",
            "source_annotation",
        ], f"{path}: unexpected header {rdr.fieldnames}"
        for row in rdr:
            n_rows += 1
            key = (row["mirna"], row["target_gene"])
            if key in seen:
                dup += 1
                continue
            seen.add(key)
            band = EV_TO_BAND[row["atomic_evidence"]]
            tiers[band].setdefault(key[0], set()).add(key[1])
    results["pair_file_checks"][org] = {
        "rows": n_rows,
        "duplicate_pairs": dup,
        "band_pairs": {b: sum(len(g) for g in m.values()) for b, m in tiers.items()},
    }
    if dup:
        raise SystemExit(f"[{org}] {dup} duplicate (mirna,gene) rows in pair file")
    return tiers


def load_universe(pairs_dir: Path, org: str, family: str) -> set[str]:
    path = pairs_dir / "enrichment_inputs" / org / f"{family}_universe.tsv"
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n")
        assert header == "gene", f"{path}: bad header {header!r}"
        return {line.rstrip("\n") for line in fh if line.strip()}


def load_term2gene(pairs_dir: Path, org: str, family: str) -> dict[str, set[str]]:
    path = pairs_dir / "enrichment_inputs" / org / f"{family}_term2gene.tsv"
    out: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n")
        assert header == "pathway_id\tgene", f"{path}: bad header {header!r}"
        for line in fh:
            pid, _, gene = line.rstrip("\n").partition("\t")
            if pid and gene:
                out.setdefault(pid, set()).add(gene)
    return out


def resolve_contract_dir(path: Path) -> Path:
    """The directory holding the ``<org>/<FAMILY>_<TIER>.tsv`` contract cells.

    ``path`` may be the release root (containing ``enrichment_significant/``) or
    the ``enrichment_significant`` directory itself.
    """
    nested = path / "enrichment_significant"
    return nested if nested.is_dir() else path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--pairs-dir", required=True, type=Path,
        help="exported artefact set to verify (the export_target_pairs.py --out dir)",
    )
    parser.add_argument(
        "--contract-dir", required=True, type=Path,
        help="frozen release root containing enrichment_significant/, OR that "
             "enrichment_significant/ directory itself (<org>/<FAMILY>_<TIER>.tsv "
             "cells with the stored overlap column)",
    )
    parser.add_argument(
        "--result", type=Path, default=Path("tp_verify_results.json"),
        help="incremental result JSON (default: ./tp_verify_results.json; "
             "never write this inside a frozen artefact dir)",
    )
    args = parser.parse_args(argv)
    contract_dir = resolve_contract_dir(args.contract_dir)
    missing_orgs = [org for org in FAMILIES_BY_ORG if not (contract_dir / org).is_dir()]
    if missing_orgs:
        parser.error(
            f"--contract-dir {args.contract_dir}: no <org>/ cell directories "
            f"{missing_orgs} in {contract_dir} (pass the release root or its "
            f"enrichment_significant/ directory)"
        )

    results: dict = {
        "status": "running",
        "contract_dir": str(contract_dir),
        "cells": {},
        "pair_file_checks": {},
    }

    def save() -> None:
        args.result.write_text(json.dumps(results, indent=2))

    total_rows = total_mism = 0
    for org in FAMILIES_BY_ORG:
        tiers = load_pairs(args.pairs_dir, org, results)
        log(org, "pairs loaded:", results["pair_file_checks"][org]["band_pairs"])
        for family in FAMILIES_BY_ORG[org]:
            universe = load_universe(args.pairs_dir, org, family)
            term2gene = load_term2gene(args.pairs_dir, org, family)
            # Pre-intersect members with universe once per pathway (formula factor
            # `term2gene[pid] & universe`; commutative reordering, identical result).
            members_u = {pid: genes & universe for pid, genes in term2gene.items()}
            # mirna_gene_lists replica on the reconstructed tiers:
            norm_tiers, _stats = pts.normalize_tiers_to_go_symbols(tiers, universe)
            for tier in TIERS:
                gene_sets = pts.cumulative_tier(norm_tiers, tier)
                gene_sets = {m: g for m, g in gene_sets.items() if g}
                tsv = contract_dir / org / f"{family}_{tier}.tsv"
                n = mism = 0
                examples: list[dict] = []
                with tsv.open("r", encoding="utf-8", newline="") as fh:
                    rdr = csv.DictReader(fh, delimiter="\t")
                    for row in rdr:
                        n += 1
                        stored = int(row["overlap"])
                        got = len(
                            members_u.get(row["pathway_id"], set())
                            & gene_sets.get(row["mirna"], set())
                        )
                        if got != stored:
                            mism += 1
                            if len(examples) < 10:
                                examples.append(
                                    {"mirna": row["mirna"],
                                     "pathway_id": row["pathway_id"],
                                     "stored": stored, "got": got}
                                )
                results["cells"][f"{org}/{family}_{tier}"] = {
                    "rows": n, "mismatches": mism, "examples": examples,
                }
                total_rows += n
                total_mism += mism
                log(f"{org}/{family}_{tier}: rows={n} mismatches={mism}")
                save()

    results["total_rows"] = total_rows
    results["total_mismatches"] = total_mism
    results["status"] = "done" if total_mism == 0 else "FAILED"
    results["wall_s"] = round(time.time() - _T0, 1)
    save()
    log("DONE", f"total_rows={total_rows} total_mismatches={total_mism}")
    return 1 if total_mism else 0


if __name__ == "__main__":
    sys.exit(main())
