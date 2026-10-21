#!/usr/bin/env python
"""Run the PRODUCTION GO ORA for rat (rno) and macaque (mml).

There is NO miRPathDB2 ground truth for rat/macaque (v2 was human+mouse only), so
this is the production enrichment on REAL GO gene-sets (``org.*.eg.db``), validated
only insofar as the clusterProfiler engine was validated on hsa+mmu. NO v2
comparison for these organisms.

For each organism this:
  * loads its target sources (rat: miRTarBase + DIANA + TargetScan + miRDB;
    macaque: DIANA + TargetScan only), recording a per-source census,
  * builds HIGH/MEDIUM/LOW target sets via the locked rollup
    (:mod:`mirpathdb_prod.sources.production_target_sets`),
  * materialises the GO TERM2GENE/TERM2NAME from the organism's annotation package
    (:mod:`mirpathdb_prod.pathways.go_genesets`), split by GO_BP/GO_CC/GO_MF,
  * re-cases every band onto the GO symbol space (the TargetScan symbol-case fix),
  * runs the BATCHED + PARALLEL clusterProfiler ORA (one R session per cell,
    ``mclapply`` over all miRNAs) for each requested (tier x GO family),
  * writes a per-(miRNA x GO-term) p-value matrix under
    ``data/processed/enrichment/<org>/`` and a build manifest with dims, runtime,
    a per-source census, and a known-miRNA sanity hit.

Generated matrices are gitignored; the manifest + census are the committed record.

Usage::

    MPD_RSCRIPT=.../mpd_cp/bin/Rscript MPD_MAX_CORES=96 \
      PYTHONPATH=src python scripts/run_production_enrichment.py \
        --organisms rno mml --tiers HIGH MEDIUM LOW
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

from mirpathdb_prod.config import load_config
from mirpathdb_prod.enrichment import clusterprofiler as cp
from mirpathdb_prod.enrichment import family_bh
from mirpathdb_prod.pathways import go_genesets as gg
from mirpathdb_prod.reporting.significant import (
    FDR_CUTOFF,
    count_retained_rows,
    write_significant_long_form,
)
from mirpathdb_prod.sources import mirdb, mirtarbase, tarbase, targetscan
from mirpathdb_prod.sources import production_target_sets as pts
from mirpathdb_prod.validation.provenance import verify_cached_sources

ORG_TAXID = {"rno": 10116, "mml": 9544}
GO_FAMILIES = ("GO_BP", "GO_CC", "GO_MF")
DEFAULT_TIERS = ("HIGH", "MEDIUM", "LOW")

# Build cache (cached raw downloads). Config-driven (paths.source_cache_dir);
# committed default /tmp/mpd-data, overridable per host via config.local.yaml.
_CFG = load_config()
CACHE_DIR = _CFG.path("source_cache_dir")
# GO annotation propagation: canonical GOALL by default (config-overridable).
GO_PROPAGATION = gg.resolve_go_propagation(_CFG.raw)

# Known-miRNA sanity probes per organism (biology the top GO terms should reflect).
SANITY_MIRNA = {
    "rno": "rno-miR-1-3p",  # canonical muscle/cardiac miRNA
    "mml": "mml-miR-1-3p",  # same family in macaque (TargetScan LOW)
}


def _cached(name: str) -> Path:
    """Locate a cached raw download under the build cache (paths.source_cache_dir)."""
    p = Path(CACHE_DIR) / name
    if not p.exists():
        raise FileNotFoundError(f"cached source not found: {p}")
    return p


def _load_sources(organism: str) -> tuple[dict[str, dict[str, set[str]]], pts.TargetSetCensus]:
    """Load all target sources for one organism and build HIGH/MEDIUM/LOW tiers."""
    taxid = ORG_TAXID[organism]
    diana_bundle = _cached(tarbase.OTHER_SPECIES_FILENAME)
    ts_path = _cached("Predicted_Targets_Context_Scores.default_predictions.txt.zip")

    print(f"[{organism}] loading DIANA-TarBase v9 (Other_species slice) ...", flush=True)
    diana_rows = tarbase.load_mti_for_taxid(diana_bundle, taxid=taxid)
    print(f"[{organism}] DIANA rows: {len(diana_rows)}", flush=True)

    print(f"[{organism}] loading TargetScan (Gene Tax ID slice) ...", flush=True)
    ts_rows = targetscan.load_predictions(ts_path, taxid=taxid)
    print(f"[{organism}] TargetScan rows: {len(ts_rows)}", flush=True)

    if organism == "rno":
        mtb_rows = mirtarbase.load_mti(_cached("mirtarbase_rno_MTI.xlsx"), taxid=taxid)
        print(f"[{organism}] miRTarBase rows: {len(mtb_rows)}", flush=True)
        mdb_rows = mirdb.load_predictions(
            _cached("miRDB_v6.0_prediction_result.txt.gz"), taxid=taxid
        )
        print(f"[{organism}] miRDB rows: {len(mdb_rows)}", flush=True)
        accs = {r["target"] for r in mdb_rows}
        print(f"[{organism}] mapping {len(accs)} miRDB RefSeq accessions -> symbol ...", flush=True)
        refseq_symbol = pts.refseq_to_symbol(accs, org_db_package=gg.org_db_for_taxid(taxid))
        # Fail loud (H10/M2): floor each rat tier (RAT_MIN_COUNTS) so a silent
        # zero-row source can never produce a thin/degenerate build.
        tiers, census = pts.build_rat_target_sets(
            mirtarbase_rows=mtb_rows,
            diana_rows=diana_rows,
            targetscan_rows=ts_rows,
            mirdb_rows=mdb_rows,
            refseq_symbol=refseq_symbol,
            min_counts=pts.RAT_MIN_COUNTS,
            min_diana_rows=pts.RAT_DIANA_MIN_ROWS,
        )
    else:  # mml: DIANA + TargetScan only (no miRTarBase / miRDB)
        # MACAQUE_MIN_COUNTS floors HIGH>=1 and LOW>=293 but EXEMPTS MEDIUM
        # (macaque MEDIUM = 0 is the documented-real truth).
        tiers, census = pts.build_macaque_target_sets(
            diana_rows=diana_rows,
            targetscan_rows=ts_rows,
            min_counts=pts.MACAQUE_MIN_COUNTS,
        )
    return tiers, census


def _write_matrix(
    out_path: Path,
    mirnas: list[str],
    terms: list[str],
    pvals: dict[str, dict[str, tuple[float, float, int]]],
) -> tuple[int, int]:
    """Write a (miRNA x GO-term) p-value matrix CSV; return (n_rows, n_cols).

    Header is the GO-term ids; each data row is the miRNA followed by p-values
    (1.0 where the term was not scored for that miRNA). Only terms scored for at
    least one miRNA become columns (keeps the matrix dense + finite).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    col_idx = {t: j for j, t in enumerate(terms)}
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["miRNA", *terms])
        for mirna in mirnas:
            row = [1.0] * len(terms)
            for term, (p, _padj, _ov) in pvals.get(mirna, {}).items():
                j = col_idx.get(term)
                if j is not None:
                    row[j] = p
            writer.writerow([mirna, *(repr(float(v)) for v in row)])
    return len(mirnas), len(terms)


def run_organism(
    organism: str, *, tiers_wanted: list[str], out_root: Path, rscript: str | None
) -> dict[str, object]:
    """Run all (tier x GO family) cells for one organism; return its manifest block."""
    taxid = ORG_TAXID[organism]
    org_t0 = time.time()

    tiers, census = _load_sources(organism)
    print(f"[{organism}] per-tier (raw symbols): {json.dumps(census.per_tier)}", flush=True)

    print(
        f"[{organism}] materialising GO TERM2GENE from {gg.org_db_for_taxid(taxid)} ...", flush=True
    )
    go = gg.build_go_genesets(taxid, rscript=rscript, propagation=GO_PROPAGATION)
    go_symbols = set()
    for fam in GO_FAMILIES:
        go_symbols |= go.universe(fam)

    # Re-case bands onto the GO symbol space (TargetScan symbol-case fix).
    norm_tiers, norm_stats = pts.normalize_tiers_to_go_symbols(tiers, go_symbols)
    print(f"[{organism}] symbol-normalization to GO space: {json.dumps(norm_stats)}", flush=True)

    out_dir = out_root / organism
    # Authoritative Deliverable B lives in a sibling tree so the matrix tree is
    # untouched: data/processed/enrichment_significant/<org>/<FAMILY>_<TIER>.tsv.
    long_dir = out_root.parent / "enrichment_significant" / organism
    cells: list[dict[str, object]] = []
    sanity: dict[str, object] = {}
    cp_provenance: dict[str, str] = {}

    # Full-family size M per GO family, fixed by the inputs before any tier runs;
    # with universe = union of the memberships every term is eligible.
    m_family_of: dict[str, int] = {}
    for family in GO_FAMILIES:
        t2g = go.family(family)
        m_family, _ = family_bh.eligible_family(t2g, go.universe(family))
        if m_family != len(t2g):
            raise ValueError(
                f"[{organism}] {family}: eligible family M={m_family} != "
                f"{len(t2g)} TERM2GENE terms -- refusing an inconsistent family"
            )
        m_family_of[family] = m_family

    # The ONLY documented-real empty cumulative tier is macaque MEDIUM (HIGH=1,
    # MEDIUM band=0 -> but HIGH carries through, so even MEDIUM is non-empty for
    # macaque); any OTHER empty tier is a silent-failure signal and must raise.
    for tier in tiers_wanted:
        gene_sets = pts.cumulative_tier(norm_tiers, tier)
        gene_sets = {m: g for m, g in gene_sets.items() if g}
        if not gene_sets:
            if organism == "mml" and tier == "MEDIUM":
                print(
                    f"[{organism}] SKIP tier {tier}: documented-real macaque MEDIUM=0 "
                    f"(its lone validated MTI is in HIGH; no MEDIUM band)",
                    flush=True,
                )
                continue
            raise ValueError(
                f"[{organism}] tier {tier}: no miRNAs with targets after rollup "
                f"(expected non-empty) -- refusing a degenerate matrix."
            )
        for family in GO_FAMILIES:
            t2g = go.family(family)
            universe = go.universe(family)
            # B2: never score an empty surface (would yield an all-p=1.0 matrix).
            if not t2g or not universe:
                raise ValueError(
                    f"[{organism}] {family}: empty TERM2GENE/universe -- refusing degenerate matrix"
                )
            t0 = time.time()
            res, prov = cp.run_ora_batch_with_provenance(
                gene_sets, t2g, universe, term2name=None, rscript=rscript
            )
            cp_provenance = cp_provenance or prov
            elapsed = time.time() - t0
            # Columns = terms scored for >=1 miRNA (keeps matrices finite + dense).
            scored_terms = sorted({term for per in res.values() for term in per})
            # B2: gene sets exist but NOTHING scored -> the R batch silently lost
            # every set (e.g. swallowed worker errors). Refuse the all-p=1.0 matrix.
            if gene_sets and not scored_terms:
                raise ValueError(
                    f"[{organism}] {family}/{tier}: {len(gene_sets)} gene set(s) but "
                    f"0 terms scored -- refusing degenerate all-p=1.0 matrix "
                    f"(likely an enricher/worker failure)."
                )
            mirnas = sorted(gene_sets)
            out_path = out_dir / f"{family}_{tier}.csv"
            n_rows, n_cols = _write_matrix(out_path, mirnas, scored_terms, res)
            # Authoritative Deliverable B (contract v2: full-family p_adjust, legacy
            # overlap-only value kept as p_adjust_overlap; raw p + overlap verbatim).
            m_family = m_family_of[family]
            long_path = long_dir / f"{family}_{tier}.tsv"
            n_sig = write_significant_long_form(long_path, mirnas, res, m_family=m_family)
            counts = count_retained_rows(res, fdr_cutoff=FDR_CUTOFF, m_family=m_family)
            cell = {
                "organism": organism,
                "tier": tier,
                "family": family,
                "n_mirnas": n_rows,
                "n_go_terms_scored": n_cols,
                "n_go_terms_total": len(t2g),
                "n_universe": len(universe),
                "adjustment_family": "full_family_bh",
                "m_family": m_family,
                "n_retained_rows": n_sig,
                # Alias read by validation/go_propagation_driver.py.
                "n_significant_rows": n_sig,
                "n_retained_rows_overlap_legacy": counts["rows_overlap_legacy"],
                "m_overlap_min": counts["m_overlap_min"],
                "m_overlap_max": counts["m_overlap_max"],
                "elapsed_s": round(elapsed, 2),
                "output": str(out_path.relative_to(out_root.parents[2])),
                "long_form_output": str(long_path.relative_to(out_root.parents[2])),
            }
            cells.append(cell)
            print(
                f"[{organism}] {tier}/{family}: {n_rows} miRNAs x {n_cols} GO terms "
                f"(universe={len(universe)}) {elapsed:.1f}s -> {out_path.name}",
                flush=True,
            )

            # Sanity hit on the first (HIGH-ish) BP cell that has the probe miRNA.
            probe = SANITY_MIRNA.get(organism)
            if probe and probe in res and family == "GO_BP" and "sanity" not in sanity:
                top = sorted(res[probe].items(), key=lambda kv: kv[1][0])[:8]
                sanity = {
                    "sanity": {
                        "tier": tier,
                        "family": family,
                        "mirna": probe,
                        "n_targets": len(gene_sets[probe]),
                        "top_go_terms": [
                            {"go_id": t, "pvalue": p, "p_adjust": padj, "overlap": ov}
                            for t, (p, padj, ov) in top
                        ],
                    }
                }

    return {
        "organism": organism,
        "taxid": taxid,
        "org_db_package": go.org_db_package,
        "go_provenance": go.provenance,
        "cp_provenance": cp_provenance,
        "per_source_census": census.per_source,
        "per_tier_raw": census.per_tier,
        "symbol_normalization": norm_stats,
        "census_notes": census.notes,
        "cells": cells,
        **sanity,
        "organism_wall_s": round(time.time() - org_t0, 1),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organisms", nargs="+", default=["rno", "mml"])
    parser.add_argument("--tiers", nargs="+", default=list(DEFAULT_TIERS))
    parser.add_argument("--rscript", default=None)
    args = parser.parse_args(argv)

    for org in args.organisms:
        if org not in ORG_TAXID:
            raise SystemExit(f"unsupported organism {org!r}; this run is rno+mml")
    for tier in args.tiers:
        if tier not in pts.TIERS:
            raise SystemExit(f"unknown tier {tier!r}; choose from {pts.TIERS}")

    cfg = load_config()
    out_root = cfg.path("data_processed_dir") / "enrichment"
    out_root.mkdir(parents=True, exist_ok=True)

    # Provenance gate: refuse to build on swapped/stale/corrupt cached sources.
    _manifest = cfg.path("provenance_dir") / "manifest.yaml"
    if _manifest.exists():
        _rep = verify_cached_sources(CACHE_DIR, _manifest)
        print(
            f"[provenance] {len(_rep['verified'])} pinned source(s) sha256-verified, "
            f"{len(_rep['absent'])} absent (left to _cached guards)",
            flush=True,
        )

    t0 = time.time()
    blocks: list[dict[str, object]] = []
    for organism in args.organisms:
        blocks.append(
            run_organism(organism, tiers_wanted=args.tiers, out_root=out_root, rscript=args.rscript)
        )
        print(f"[{organism}] done in {blocks[-1]['organism_wall_s']}s", flush=True)

    # Merge into any existing production manifest (one organism may be run at a
    # time): organisms run THIS invocation replace their prior block.
    manifest_path = out_root / "production_manifest.json"
    run_orgs = set(args.organisms)
    kept: list[dict[str, object]] = []
    if manifest_path.exists():
        try:
            prior = json.loads(manifest_path.read_text(encoding="utf-8"))
            kept = [b for b in prior.get("organisms", []) if b.get("organism") not in run_orgs]
        except (OSError, json.JSONDecodeError):
            kept = []
    manifest = {
        "build": "production_go_ora",
        "note": "Real GO gene-sets (org.*.eg.db); NO miRPathDB2 ground truth for rno/mml.",
        "last_run_organisms": sorted(run_orgs),
        "last_run_tiers": args.tiers,
        "last_run_wall_s": round(time.time() - t0, 1),
        "organisms": kept + blocks,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\nTOTAL wall: {time.time() - t0:.1f}s; manifest -> {manifest_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
