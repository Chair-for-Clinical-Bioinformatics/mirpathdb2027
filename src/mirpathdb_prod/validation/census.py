"""Data-contract / census-drift validation for the production atlas.

Two read-only checks that together form the *reconciliation gate* — they prove the
already-built artefacts are real and complete, so no recomputation is needed:

* :func:`validate_matrices_nondegenerate` — scans every enrichment matrix and flags
  any that is **degenerate** (header-only, or every cell == 1.0). A degenerate
  matrix is the silent-failure signature of an empty TERM2GENE / all-failing ORA
  (see the B2 guard); a non-degenerate matrix has real p<1.0 enrichment values.
* :func:`validate_census` — compares the per-source / per-tier counts recorded in
  the run manifests against pinned **census anchors** (e.g. macaque == exactly 1
  validated MTI; macaque LOW == 293 TargetScan miRNAs; inferred macaque Reactome ==
  2777 pathways; rat DIANA == 1125 MTIs), failing on drift.

Census anchors default to the structural values below but can be overridden via a
``census_contract:`` block in ``config.yaml`` (so they are not hard-coded into a
gate). Volume counts (large miRNA/pair totals) are deliberately NOT anchored here —
only the structural invariants that must hold exactly.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

# Structural census anchors that must hold exactly (override via config census_contract).
# Each entry: (manifest_file, dotted-path-into-json, expected-value).
DEFAULT_CENSUS_ANCHORS: dict[str, Any] = {
    "macaque_inferred_reactome_n_terms": (
        "macaque_inferred_reactome_manifest.json",
        "inferred_reactome_summary.n_terms",
        2777,
    ),
    "macaque_HIGH_mirnas": (
        "macaque_inferred_reactome_manifest.json",
        "target_set_census.HIGH.n_mirnas",
        1,
    ),
    "macaque_MEDIUM_mirnas": (
        "macaque_inferred_reactome_manifest.json",
        "target_set_census.MEDIUM.n_mirnas",
        0,
    ),
    "macaque_LOW_mirnas": (
        "macaque_inferred_reactome_manifest.json",
        "target_set_census.LOW.n_mirnas",
        293,
    ),
}

# Per-organism (production_enrichment_manifest) anchors keyed by organism name.
DEFAULT_ORGANISM_ANCHORS: dict[str, dict[str, int]] = {
    # rat DIANA volume + macaque lone validated MTI (the real-vs-fabricated anchor).
    "rno": {"per_source_census.diana_tarbase.n_rows": 1125},
    "mml": {"per_source_census.diana_tarbase.n_rows": 1},
}


def _dig(obj: Any, dotted: str) -> Any:
    """Follow a dotted path into nested dicts; return None if any hop is missing."""
    cur = obj
    for key in dotted.split("."):
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def scan_matrix(path: Path) -> dict:
    """Scan one (miRNA x pathway) p-value matrix CSV for degeneracy.

    Returns ``{path, n_rows, n_cols, non_degenerate}``. A matrix is non-degenerate
    if it has at least one pathway column AND at least one cell < 1.0 (a real
    enrichment p-value); header-only or all-1.0 matrices are degenerate. Streams
    and early-exits as soon as a value < 1.0 is found.
    """
    n_cols = 0
    n_rows = 0
    non_degenerate = False
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return {"path": str(path), "n_rows": 0, "n_cols": 0, "non_degenerate": False}
        n_cols = max(0, len(header) - 1)  # minus the leading "miRNA" column
        for row in reader:
            n_rows += 1
            if not non_degenerate:
                for cell in row[1:]:
                    try:
                        if float(cell) < 1.0:
                            non_degenerate = True
                            break
                    except ValueError:
                        continue
    return {
        "path": str(path),
        "n_rows": n_rows,
        "n_cols": n_cols,
        "non_degenerate": non_degenerate and n_cols > 0,
    }


def validate_matrices_nondegenerate(
    enrichment_dir: Path, *, require_present: bool = False
) -> dict:
    """Scan every ``*.csv`` enrichment matrix; flag degenerate ones.

    Returns ``{ok, n_matrices, degenerate: [...], scanned: [...], empty: bool}``.
    ``ok`` is True iff no matrix is degenerate. When ``require_present`` is True
    (the production gate), an enrichment dir that contains NO matrices at all is a
    hard failure (``ok=False``) rather than a vacuous pass -- a "gate" that passes
    because nothing was built is exactly the silent success this check exists to
    catch.
    """
    enrichment_dir = Path(enrichment_dir)
    scanned: list[dict] = []
    degenerate: list[dict] = []
    for csv_path in sorted(enrichment_dir.glob("*/*.csv")):
        result = scan_matrix(csv_path)
        scanned.append(result)
        if not result["non_degenerate"]:
            degenerate.append(result)
    empty = not scanned
    return {
        "ok": not degenerate and not (require_present and empty),
        "n_matrices": len(scanned),
        "degenerate": degenerate,
        "scanned": scanned,
        "empty": empty,
    }


def validate_census(
    provenance_dir: Path, anchors: dict | None = None, *, require_present: bool = False
) -> dict:
    """Compare run-manifest census counts against pinned structural anchors.

    Returns ``{ok, checks: [...], failures: [...]}``. ``ok`` is True iff every
    anchor matches. By default a missing manifest / organism block / field is
    recorded as a ``skip`` (not a failure), since not every run writes every
    manifest. When ``require_present`` is True (the production gate), a missing
    anchor manifest, a missing ``production_enrichment_manifest.json``, or a
    missing per-organism block is a HARD FAILURE -- so the reconciliation gate
    cannot "pass" merely because the artefacts it should reconcile were never
    produced (the vacuous-pass hole).
    """
    provenance_dir = Path(provenance_dir)
    anchors = anchors if anchors is not None else DEFAULT_CENSUS_ANCHORS
    checks: list[dict] = []
    failures: list[dict] = []

    def _absent(rec: dict) -> None:
        # Under require_present a missing artefact fails the gate; otherwise skip.
        (failures if require_present else checks).append(
            {**rec, "status": "FAIL" if require_present else "skip"}
        )

    # Simple (manifest, dotted-path, expected) anchors.
    for name, (mfile, dotted, expected) in anchors.items():
        mpath = provenance_dir / mfile
        if not mpath.exists():
            _absent({"anchor": name, "reason": f"{mfile} absent"})
            continue
        doc = json.loads(mpath.read_text(encoding="utf-8"))
        got = _dig(doc, dotted)
        rec = {"anchor": name, "expected": expected, "got": got}
        if got == expected:
            checks.append({**rec, "status": "ok"})
        else:
            failures.append({**rec, "status": "FAIL"})

    # Per-organism anchors in production_enrichment_manifest.json.
    epath = provenance_dir / "production_enrichment_manifest.json"
    if not epath.exists():
        _absent({"anchor": "production_enrichment_manifest", "reason": "manifest absent"})
    else:
        edoc = json.loads(epath.read_text(encoding="utf-8"))
        by_org = {o.get("organism"): o for o in edoc.get("organisms", [])}
        for org, org_anchors in DEFAULT_ORGANISM_ANCHORS.items():
            block = by_org.get(org)
            if block is None:
                _absent({"anchor": f"{org}", "reason": "organism absent"})
                continue
            for dotted, expected in org_anchors.items():
                got = _dig(block, dotted)
                rec = {"anchor": f"{org}.{dotted}", "expected": expected, "got": got}
                if got == expected:
                    checks.append({**rec, "status": "ok"})
                else:
                    failures.append({**rec, "status": "FAIL"})

    return {"ok": not failures, "checks": checks, "failures": failures}


def _main(argv: list[str] | None = None) -> int:
    """CLI: run both checks over the repo and print a report (exit 0 = PASS)."""
    import argparse

    from mirpathdb_prod.config import load_config

    parser = argparse.ArgumentParser(description="Validate atlas census + matrix non-degeneracy.")
    parser.add_argument("--enrichment-dir", default=None)
    parser.add_argument("--provenance-dir", default=None)
    parser.add_argument(
        "--require-present",
        action="store_true",
        help=(
            "Production gate: an empty enrichment dir or a missing anchor manifest "
            "is a hard FAILURE instead of a vacuous pass (use in CI / before publish)."
        ),
    )
    args = parser.parse_args(argv)

    cfg = load_config()
    enrichment_dir = (
        Path(args.enrichment_dir)
        if args.enrichment_dir
        else cfg.path("data_processed_dir") / "enrichment"
    )
    provenance_dir = (
        Path(args.provenance_dir) if args.provenance_dir else cfg.path("provenance_dir")
    )

    # Allow config override of the structural anchors.
    contract = cfg.raw.get("census_contract")
    anchors = contract if isinstance(contract, dict) and contract else None

    mat = validate_matrices_nondegenerate(enrichment_dir, require_present=args.require_present)
    cen = validate_census(provenance_dir, anchors=anchors, require_present=args.require_present)

    print(f"[census] matrices scanned: {mat['n_matrices']}; degenerate: {len(mat['degenerate'])}")
    if args.require_present and mat.get("empty"):
        print(f"  EMPTY: no enrichment matrices under {enrichment_dir} (require-present -> FAIL)")
    for d in mat["degenerate"]:
        print(f"  DEGENERATE: {d['path']} ({d['n_rows']}x{d['n_cols']})")
    print(f"[census] anchors: {len(cen['checks'])} checked; {len(cen['failures'])} failed")
    for f in cen["failures"]:
        print(f"  FAIL: {f['anchor']} expected={f['expected']} got={f['got']}")

    ok = mat["ok"] and cen["ok"]
    print(f"[census] {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_main())
