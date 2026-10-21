#!/usr/bin/env python
"""Run the INFERRED macaque Reactome ORA (orthology-projected from human).

Macaque has **no curated Reactome** -- both ``reactome.db`` v1.95.0 and the flat
``Ensembl2Reactome_All_Levels.txt`` carry 0 ``R-MML-`` rows. But Reactome itself
lists ``Macaca mulatta`` as an **orthology-INFERRED** species, so this runner
builds macaque Reactome the way Reactome does: it projects every human
``R-HSA-`` pathway onto macaque through the NCBI ``gene_orthologs`` 9606<->9544
strict 1:1 ortholog map (mapped to macaque SYMBOL via ``org.Mmu.eg.db``) and
relabels each pathway ``R-HSA-<n>`` -> ``R-MML-<n>``. See
:mod:`mirpathdb_prod.pathways.inferred_reactome_genesets`.

This is a macaque-only companion to :mod:`scripts.run_production_pathways` (which
deliberately SKIPS Reactome for macaque because the native source has no R-MML-).
It runs ORA via the SAME batched + parallel kernel
(:func:`mirpathdb_prod.enrichment.clusterprofiler.run_ora_batch`, ``mclapply``),
over the macaque production target sets (LOW = 293 TargetScan miRNAs; HIGH = the
single marked validated MTI; MEDIUM as applicable), and writes
``data/processed/enrichment/mml/prod_REACTOME_<TIER>.csv`` (same naming + miRNA
rows x pathway columns orientation as every other prod matrix).

The macaque GO matrices are NOT touched -- this only ADDS the inferred Reactome
family. Every output is MARKED INFERRED in the manifest
(``provenance/macaque_inferred_reactome_manifest.json`` -- the committed record;
the matrices themselves are gitignored).

Usage::

    MPD_RSCRIPT=.../mpd_cp/bin/Rscript MPD_MAX_CORES=96 \
      PYTHONPATH=src python scripts/run_macaque_inferred_reactome.py \
        --tiers HIGH MEDIUM LOW
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
from mirpathdb_prod.pathways import inferred_reactome_genesets as ir
from mirpathdb_prod.reporting.significant import (
    FDR_CUTOFF,
    count_retained_rows,
    write_significant_long_form,
)
from mirpathdb_prod.sources import production_target_sets as pts
from mirpathdb_prod.sources import tarbase, targetscan
from mirpathdb_prod.validation.provenance import verify_cached_sources

MACAQUE = "mml"
MACAQUE_TAXID = 9544

# Build cache (cached raw downloads + the NCBI gene_orthologs snapshot).
# Config-driven (paths.source_cache_dir); committed default /tmp/mpd-data,
# overridable per host via config.local.yaml.
CACHE_DIR = load_config().path("source_cache_dir")
# The cached NCBI gene_orthologs.gz lives under the same config-driven cache.
DEFAULT_ORTHOLOG_PATH = Path(CACHE_DIR) / "gene_orthologs.gz"
DEFAULT_TIERS = ("HIGH", "MEDIUM", "LOW")

# Pinned NCBI ortholog snapshot date for THIS cached file (fetched 2026-06-10 per
# provenance/manifest.yaml ncbi_gene; the live NCBI file has since been refreshed,
# so the build pins the CACHED file's sha256 + size, recorded by the builder).
ORTHOLOG_LAST_MODIFIED = "cached 2026-06-10 (NCBI gene_orthologs.gz; size+sha256 pinned)"

# A known macaque miRNA sanity probe -> a sensible inferred Reactome pathway.
SANITY_MIRNA = "mml-let-7a-5p"


def _cached(name: str) -> Path:
    p = Path(CACHE_DIR) / name
    if not p.exists():
        raise FileNotFoundError(f"cached source not found: {p}")
    return p


def _load_macaque_target_sets() -> tuple[dict[str, dict[str, set[str]]], pts.TargetSetCensus]:
    """Build macaque HIGH/MEDIUM/LOW target sets (raw symbol space).

    Identical to the run_production_pathways macaque path: DIANA Other_species
    bundle (verified -> HIGH, CLIP -> MEDIUM; macaque volume = 1 MTI) + TargetScan
    (-> LOW, 293 miRNAs).
    """
    diana_bundle = _cached(tarbase.OTHER_SPECIES_FILENAME)
    ts_path = _cached("Predicted_Targets_Context_Scores.default_predictions.txt.zip")
    diana_rows = tarbase.load_mti_for_taxid(diana_bundle, taxid=MACAQUE_TAXID)
    ts_rows = targetscan.load_predictions(ts_path, taxid=MACAQUE_TAXID)
    # Fail loud (H10/M2): floor HIGH>=1 + LOW>=293; EXEMPT MEDIUM (documented-real 0).
    return pts.build_macaque_target_sets(
        diana_rows=diana_rows,
        targetscan_rows=ts_rows,
        min_counts=pts.MACAQUE_MIN_COUNTS,
    )


def _write_matrix(
    out_path: Path,
    mirnas: list[str],
    terms: list[str],
    pvals: dict[str, dict[str, tuple[float, float, int]]],
) -> tuple[int, int]:
    """Write a (miRNA x pathway) p-value matrix CSV; return (n_rows, n_cols).

    Same orientation + format as scripts/run_production_pathways._write_matrix:
    miRNA rows x scored-pathway columns, p=1.0 where a pathway was not scored.
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tiers", nargs="+", default=list(DEFAULT_TIERS))
    parser.add_argument("--ortholog-path", default=str(DEFAULT_ORTHOLOG_PATH))
    parser.add_argument("--rscript", default=None)
    args = parser.parse_args(argv)

    for tier in args.tiers:
        if tier not in pts.TIERS:
            raise SystemExit(f"unknown tier {tier!r}; choose from {pts.TIERS}")

    cfg = load_config()
    out_root = cfg.path("data_processed_dir") / "enrichment"
    out_dir = out_root / MACAQUE
    out_dir.mkdir(parents=True, exist_ok=True)

    # Provenance gate: refuse to build on swapped/stale/corrupt cached sources
    # (the load-bearing input here is the pinned NCBI gene_orthologs.gz).
    _manifest = cfg.path("provenance_dir") / "manifest.yaml"
    if _manifest.exists():
        _rep = verify_cached_sources(CACHE_DIR, _manifest)
        print(
            f"[provenance] {len(_rep['verified'])} pinned source(s) sha256-verified, "
            f"{len(_rep['absent'])} absent (left to _cached guards)",
            flush=True,
        )

    t0 = time.time()

    print(
        "[mml] building INFERRED macaque Reactome (human R-HSA- -> R-MML- via orthology) ...",
        flush=True,
    )
    gs = ir.build_inferred_macaque_reactome(
        ortholog_path=args.ortholog_path,
        ortholog_last_modified=ORTHOLOG_LAST_MODIFIED,
        rscript=args.rscript,
    )
    t2g = gs.term2gene
    term2name = gs.term2name
    universe = gs.universe()
    summary = gs.summary()
    print(
        f"[mml] INFERRED Reactome: {summary['n_terms']} pathways, "
        f"{summary['n_genes']} macaque-symbol universe, {summary['n_pairs']} pairs",
        flush=True,
    )

    print("[mml] building macaque production target sets ...", flush=True)
    tiers, census = _load_macaque_target_sets()
    print(f"[mml] per-tier (raw symbols): {json.dumps(census.per_tier)}", flush=True)

    # Re-case the target-set bands onto the inferred-Reactome macaque SYMBOL space
    # (the same TargetScan symbol-case fix the GO/Reactome runners apply).
    norm_tiers, norm_stats = pts.normalize_tiers_to_go_symbols(tiers, set(universe))

    cells: list[dict[str, object]] = []
    sanity: list[dict[str, object]] = []
    # B2: never score an empty inferred-Reactome surface (would be all-p=1.0).
    if not t2g or not universe:
        raise ValueError(
            "[mml] REACTOME(INFERRED): empty TERM2GENE/universe -- refusing degenerate matrix"
        )
    # Full-family size M is fixed by the inputs; with universe = union of the
    # memberships every term is eligible, so anything else is a build defect.
    m_family, _ = family_bh.eligible_family(t2g, universe)
    if m_family != len(t2g):
        raise ValueError(
            f"[mml] REACTOME(INFERRED): eligible family M={m_family} != {len(t2g)} "
            f"TERM2GENE terms -- refusing an inconsistent family"
        )
    for tier in args.tiers:
        gene_sets = pts.cumulative_tier(norm_tiers, tier)
        gene_sets = {m: g for m, g in gene_sets.items() if g}
        if not gene_sets:
            # The only documented-real empty tier is macaque MEDIUM (HIGH=1,
            # MEDIUM band=0); any other empty tier is a silent-failure signal.
            if tier == "MEDIUM":
                print(
                    "[mml] SKIP REACTOME/MEDIUM: documented-real macaque MEDIUM=0 "
                    "(lone validated MTI is in HIGH)",
                    flush=True,
                )
                continue
            raise ValueError(
                f"[mml] REACTOME/{tier}: no miRNAs with targets (expected "
                f"non-empty) -- refusing a degenerate matrix."
            )
        c0 = time.time()
        # term2name=None so the matrix columns are R-MML- pathway IDs, exactly
        # matching the other prod_REACTOME_<TIER>.csv matrices (hsa = R-HSA-, rno =
        # R-RNO-). The human pathway NAMES are carried in the inferred-Reactome
        # provenance/TERM2NAME, not the matrix columns.
        res, cp_prov = cp.run_ora_batch_with_provenance(
            gene_sets, t2g, universe, term2name=None, rscript=args.rscript
        )
        elapsed = time.time() - c0
        scored_terms = sorted({term for per in res.values() for term in per})
        # B2: many gene sets but NOTHING scored -> the R batch silently lost every
        # set; refuse the all-p=1.0 matrix. EXEMPT the legitimately-tiny macaque
        # HIGH/MEDIUM cells (HIGH = the single validated MTI, 1 miRNA x few genes:
        # 0 scored terms there is real biology, not a failure).
        if gene_sets and not scored_terms and tier not in ("HIGH", "MEDIUM"):
            raise ValueError(
                f"[mml] REACTOME/{tier}: {len(gene_sets)} gene set(s) but 0 terms "
                f"scored -- refusing degenerate all-p=1.0 matrix "
                f"(likely an enricher/worker failure)."
            )
        mirnas = sorted(gene_sets)
        out_path = out_dir / f"prod_REACTOME_{tier}.csv"
        n_rows, n_cols = _write_matrix(out_path, mirnas, scored_terms, res)
        # Authoritative Deliverable B for inferred macaque Reactome (contract v2:
        # full-family p_adjust, legacy overlap-only value kept as p_adjust_overlap;
        # raw p + overlap verbatim); INFERRED R-MML- ids resolve to inferred=true
        # in pathway_labels.tsv.
        long_path = out_root.parent / "enrichment_significant" / MACAQUE / f"REACTOME_{tier}.tsv"
        n_sig = write_significant_long_form(long_path, mirnas, res, m_family=m_family)
        counts = count_retained_rows(res, fdr_cutoff=FDR_CUTOFF, m_family=m_family)
        cells.append(
            {
                "organism": MACAQUE,
                "family": "REACTOME",
                "inferred": True,
                "tier": tier,
                "n_mirnas": n_rows,
                "n_pathways_scored": n_cols,
                "n_pathways_total": len(t2g),
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
                "cp_provenance": cp_prov,
                "output": f"data/processed/enrichment/{MACAQUE}/prod_REACTOME_{tier}.csv",
                "long_form_output": (
                    f"data/processed/enrichment_significant/{MACAQUE}/REACTOME_{tier}.tsv"
                ),
            }
        )
        print(
            f"[mml] REACTOME(INFERRED)/{tier}: {n_rows} miRNAs x {n_cols} pathways "
            f"(universe={len(universe)}) {elapsed:.1f}s -> {out_path.name}",
            flush=True,
        )

        # Sanity: a known macaque miRNA -> its lowest-p inferred Reactome pathways
        # (when present in this tier; let-7 is a TargetScan/LOW miRNA). res is keyed
        # by R-MML- pathway ID (term2name=None); map back to the human pathway NAME
        # via the inferred TERM2NAME for a human-readable sanity record.
        if SANITY_MIRNA in res and res[SANITY_MIRNA]:
            top = sorted(res[SANITY_MIRNA].items(), key=lambda kv: kv[1][0])[:8]
            sanity.append(
                {
                    "tier": tier,
                    "mirna": SANITY_MIRNA,
                    "n_targets": len(gene_sets[SANITY_MIRNA]),
                    "top_pathways": [
                        {
                            "pathway_id": pid,
                            "pathway_name": term2name.get(pid, pid),
                            "pvalue": p,
                            "p_adjust": padj,
                            "overlap": ov,
                        }
                        for pid, (p, padj, ov) in top
                    ],
                }
            )

    manifest = {
        "build": "macaque_inferred_reactome_ora",
        "inferred": True,
        "note": (
            "INFERRED macaque Reactome: human R-HSA- Reactome pathways projected to "
            "macaque via NCBI gene_orthologs 9606<->9544 strict 1:1 orthology, mapped "
            "to macaque SYMBOL (org.Mmu.eg.db), relabelled R-HSA-<n> -> R-MML-<n>. This "
            "replicates Reactome's OWN inference method for R-MML- (Reactome lists "
            "Macaca mulatta as orthology-inferred / no curated pathways). LIMITED + "
            "orthology-projected, NOT curated. Macaque GO matrices are NOT touched. "
            "WikiPathways remains NONE for macaque (no Macaca_mulatta collection)."
        ),
        "organism": MACAQUE,
        "taxid": MACAQUE_TAXID,
        "family": "REACTOME",
        "matrix_naming": f"data/processed/enrichment/{MACAQUE}/prod_REACTOME_<TIER>.csv",
        "inferred_reactome_provenance": gs.provenance,
        "inferred_reactome_summary": summary,
        "target_set_census": census.per_tier,
        "per_source_census": census.per_source,
        "symbol_normalization": norm_stats,
        "tiers": list(args.tiers),
        "cells": cells,
        "sanity": sanity,
        "wall_s": round(time.time() - t0, 1),
    }
    payload = json.dumps(manifest, indent=2)

    manifest_path = out_root / "macaque_inferred_reactome_manifest.json"
    manifest_path.write_text(payload, encoding="utf-8")

    prov_dir = cfg.path("provenance_dir")
    prov_dir.mkdir(parents=True, exist_ok=True)
    prov_copy = prov_dir / "macaque_inferred_reactome_manifest.json"
    prov_copy.write_text(payload, encoding="utf-8")

    print(
        f"\n[mml] TOTAL wall: {time.time() - t0:.1f}s; manifest -> {manifest_path}\n"
        f"committed copy -> {prov_copy}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
