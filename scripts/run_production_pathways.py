#!/usr/bin/env python
"""Run the UNIFORM PRODUCTION pathway ORA across Reactome + WikiPathways + GO.

This is the cross-family companion to :mod:`scripts.run_production_enrichment`
(which built production GO for rat/macaque). It extends the production enrichment
to the new pathway families on REAL gene sets, in gene-SYMBOL space, for all four
organisms where each family exists:

  * **Reactome** (``reactome.db`` x ``org.*.eg.db``): hsa, mmu, rno. Reactome
    ships NO curated ``R-MML-`` pathways, so this NATIVE runner emits no macaque
    Reactome. Macaque Reactome is built separately as INFERRED (human ``R-HSA-``
    orthology-projected to ``R-MML-``) by ``scripts/run_macaque_inferred_reactome.py``.
  * **WikiPathways** (monthly GMT x ``org.*.eg.db``): hsa, mmu, rno. The release
    ships no ``Macaca_mulatta`` collection, so macaque gets no WikiPathways (there
    is no orthology-inference path for WikiPathways, so macaque WP stays NONE).
  * **Production GO** (``org.*.eg.db`` GO_BP/GO_CC/GO_MF): hsa, mmu ONLY. Rat +
    macaque production GO already exist (built by run_production_enrichment.py)
    and are NOT re-run here -- this runner refuses the GO family for rno/mml so
    those matrices can never be clobbered.

Gene-set families share ONE numerical contract (clusterProfiler ``enricher`` via
the batched + parallel :func:`mirpathdb_prod.enrichment.clusterprofiler.run_ora_batch`
-- one R session per cell, ``mclapply`` over all miRNAs). TERM2GENE is built once
per (organism, family); the production target sets are built once per organism and
re-cased onto each family's own symbol space (the TargetScan symbol-case fix).

Output naming (NEVER collides with the existing matrices)
---------------------------------------------------------
Every matrix this runner writes is prefixed ``prod_`` so it is distinct from:
  * the v2-COMPARISON hsa/mmu matrices (``GO_BP_validated_miRTarBase_*.csv`` etc.),
  * the existing production rno/mml GO matrices (``GO_BP_HIGH.csv`` etc.).
Path: ``data/processed/enrichment/<org>/prod_<FAMILY>_<TIER>.csv``
(e.g. ``hsa/prod_REACTOME_HIGH.csv``, ``mmu/prod_GO_BP_HIGH.csv``). Orientation is
miRNA rows x pathway columns (same as every other matrix in the tree). Matrices
are gitignored; the committed record is ``production_pathways_manifest.json``.

Usage::

    MPD_RSCRIPT=.../mpd_cp/bin/Rscript MPD_MAX_CORES=96 \
      PYTHONPATH=src python scripts/run_production_pathways.py \
        --organisms hsa mmu rno mml --families REACTOME WIKIPATHWAYS GO \
        --tiers HIGH
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
from mirpathdb_prod.pathways import reactome_genesets as rg
from mirpathdb_prod.pathways import wikipathways_genesets as wp
from mirpathdb_prod.reporting.significant import (
    FDR_CUTOFF,
    count_retained_rows,
    write_significant_long_form,
)
from mirpathdb_prod.sources import mirdb, mirtarbase, tarbase, targetscan
from mirpathdb_prod.sources import production_target_sets as pts
from mirpathdb_prod.validation.provenance import verify_cached_sources

ORG_TAXID = {"hsa": 9606, "mmu": 10090, "rno": 10116, "mml": 9544}

# Family tokens this runner can emit. "GO" expands to the three GO ontologies.
GO_SUBFAMILIES = ("GO_BP", "GO_CC", "GO_MF")
FAMILY_TOKENS = ("REACTOME", "WIKIPATHWAYS", "GO")
DEFAULT_FAMILIES = ("REACTOME", "WIKIPATHWAYS")
DEFAULT_TIERS = ("HIGH",)

# Build cache holding the cached raw downloads + the WikiPathways GMT cache.
# Config-driven (paths.source_cache_dir); committed default /tmp/mpd-data,
# overridable per host via config.local.yaml.
_CFG = load_config()
CACHE_DIR = _CFG.path("source_cache_dir")
# GO annotation propagation: canonical GOALL by default (config-overridable via
# config.yaml go_propagation: direct).
GO_PROPAGATION = gg.resolve_go_propagation(_CFG.raw)

# Known-miRNA sanity probes per organism (a let-7 / miR-1 family member whose top
# pathway in a new family should be biologically sensible).
SANITY_MIRNA = {
    "hsa": "hsa-miR-1-3p",
    "mmu": "mmu-miR-1a-3p",
    "rno": "rno-miR-1-3p",
    "mml": "mml-miR-1-3p",
}

# Production GO is rat/macaque-complete already; this runner only ADDS hsa/mmu GO.
GO_ALLOWED_ORGS = {"hsa", "mmu"}


def _cached(name: str) -> Path:
    """Locate a cached raw download under the build cache."""
    p = Path(CACHE_DIR) / name
    if not p.exists():
        raise FileNotFoundError(f"cached source not found: {p}")
    return p


def _load_target_sets(organism: str) -> tuple[dict[str, dict[str, set[str]]], pts.TargetSetCensus]:
    """Build HIGH/MEDIUM/LOW production target sets (raw symbol space) for one org.

    hsa + mmu ship per-species miRTarBase + DIANA files, so they go through the
    all-four-source materialiser :func:`pts.load_and_build_target_sets`. rno + mml
    live in the DIANA Other_species bundle, so they are assembled exactly like the
    existing GO runner (rat = four sources; macaque = DIANA + TargetScan only).
    """
    taxid = ORG_TAXID[organism]
    if organism in ("hsa", "mmu"):
        return pts.load_and_build_target_sets(taxid, cache_dir=CACHE_DIR)

    # rno / mml -- Other_species DIANA bundle.
    diana_bundle = _cached(tarbase.OTHER_SPECIES_FILENAME)
    ts_path = _cached("Predicted_Targets_Context_Scores.default_predictions.txt.zip")
    diana_rows = tarbase.load_mti_for_taxid(diana_bundle, taxid=taxid)
    ts_rows = targetscan.load_predictions(ts_path, taxid=taxid)
    if organism == "rno":
        mtb_rows = mirtarbase.load_mti(_cached("mirtarbase_rno_MTI.xlsx"), taxid=taxid)
        mdb_rows = mirdb.load_predictions(
            _cached("miRDB_v6.0_prediction_result.txt.gz"), taxid=taxid
        )
        accs = {r["target"] for r in mdb_rows}
        refseq_symbol = pts.refseq_to_symbol(accs, org_db_package=gg.org_db_for_taxid(taxid))
        # Fail loud (H10/M2): floor each rat tier + the rat DIANA volume.
        return pts.build_rat_target_sets(
            mirtarbase_rows=mtb_rows,
            diana_rows=diana_rows,
            targetscan_rows=ts_rows,
            mirdb_rows=mdb_rows,
            refseq_symbol=refseq_symbol,
            min_counts=pts.RAT_MIN_COUNTS,
            min_diana_rows=pts.RAT_DIANA_MIN_ROWS,
        )
    # Macaque: floor HIGH>=1 + LOW>=293, EXEMPT MEDIUM (documented-real 0).
    return pts.build_macaque_target_sets(
        diana_rows=diana_rows,
        targetscan_rows=ts_rows,
        min_counts=pts.MACAQUE_MIN_COUNTS,
    )


def _resolve_families(organism: str, families_wanted: list[str]) -> list[str]:
    """Expand the requested family tokens to per-organism (sub)family tokens.

    ``GO`` expands to GO_BP/GO_CC/GO_MF and is ONLY allowed for hsa/mmu (rat +
    macaque production GO already exist). Reactome / WikiPathways are dropped for
    organisms with no coverage (macaque) so this runner never errors on them.
    """
    taxid = ORG_TAXID[organism]
    out: list[str] = []
    for token in families_wanted:
        if token == "GO":
            if organism not in GO_ALLOWED_ORGS:
                print(
                    f"[{organism}] SKIP GO: production GO already exists for rno/mml "
                    f"(not re-run/clobbered).",
                    flush=True,
                )
                continue
            out.extend(GO_SUBFAMILIES)
        elif token == "REACTOME":
            if taxid in rg.REACTOME_SPECIES_CODE:
                out.append("REACTOME")
            else:
                print(
                    f"[{organism}] SKIP REACTOME: no native/curated Reactome (macaque); "
                    f"inferred macaque Reactome runs via run_macaque_inferred_reactome.py.",
                    flush=True,
                )
        elif token == "WIKIPATHWAYS":
            if taxid in wp.WP_SPECIES_FILE:
                out.append("WIKIPATHWAYS")
            else:
                print(f"[{organism}] SKIP WIKIPATHWAYS: no WP collection (macaque).", flush=True)
        else:
            raise SystemExit(f"unknown family token {token!r}; choose from {FAMILY_TOKENS}")
    return out


def _build_family_sets(
    organism: str, family: str, *, rscript: str | None
) -> tuple[dict[str, set[str]], set[str], dict[str, object]]:
    """Materialise (TERM2GENE, universe, provenance) for one (organism, family).

    Returns the ``pathway_id -> set(symbol)`` map, the family universe (union of
    all member symbols), and a provenance dict for the manifest. Built ONCE per
    (organism, family) and reused across every tier.
    """
    taxid = ORG_TAXID[organism]
    if family in GO_SUBFAMILIES:
        go = _build_family_sets._go_cache.get(organism)  # type: ignore[attr-defined]
        if go is None:
            go = gg.build_go_genesets(taxid, rscript=rscript, propagation=GO_PROPAGATION)
            _build_family_sets._go_cache[organism] = go  # type: ignore[attr-defined]
        t2g = go.family(family)
        return t2g, go.universe(family), dict(go.provenance)
    if family == "REACTOME":
        rs = rg.build_reactome_genesets(taxid, rscript=rscript)
        return rs.term2gene, rs.universe(), dict(rs.provenance)
    if family == "WIKIPATHWAYS":
        ws = wp.build_wikipathways_genesets(taxid, cache_dir=CACHE_DIR, rscript=rscript)
        return ws.term2gene, ws.universe(), dict(ws.provenance)
    raise SystemExit(f"unknown family {family!r}")


# Per-organism GO cache so we build org.*.eg.db GO once even across GO_BP/CC/MF.
_build_family_sets._go_cache = {}  # type: ignore[attr-defined]


def _write_matrix(
    out_path: Path,
    mirnas: list[str],
    terms: list[str],
    pvals: dict[str, dict[str, tuple[float, float, int]]],
) -> tuple[int, int]:
    """Write a (miRNA x pathway) p-value matrix CSV; return (n_rows, n_cols).

    Header is the pathway ids; each data row is the miRNA followed by p-values
    (1.0 where the pathway was not scored for that miRNA). Only pathways scored
    for at least one miRNA become columns (keeps the matrix dense + finite).
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
    organism: str,
    *,
    families_wanted: list[str],
    tiers_wanted: list[str],
    out_root: Path,
    rscript: str | None,
    long_form: bool = False,
) -> dict[str, object]:
    """Run all (family x tier) cells for one organism; return its manifest block.

    When ``long_form`` is set, each cell also writes a significant-subset long-form
    TSV (``<org>/<FAMILY>_<TIER>.tsv`` under ``enrichment_significant/``) in the
    contract v2 layout: ``p_adjust`` is the full-family BH value over the cell's
    ``M`` eligible terms (``m_family``), the enricher's overlap-only value is kept
    as ``p_adjust_overlap``, and raw ``p_value`` + ``overlap`` come straight from
    the ``(p, padj, overlap)`` tuples ``run_ora_batch`` returns. The matrix output
    is unchanged.
    """
    taxid = ORG_TAXID[organism]
    org_t0 = time.time()

    families = _resolve_families(organism, families_wanted)
    if not families:
        print(f"[{organism}] no eligible families ({families_wanted}); skipping", flush=True)
        return {
            "organism": organism,
            "taxid": taxid,
            "note": (
                "no eligible new NATIVE families (e.g. macaque: no native Reactome / "
                "no WP collection, GO already exists). Macaque Reactome is provided "
                "INFERRED via run_macaque_inferred_reactome.py (orthology-projected)."
            ),
            "families": [],
            "cells": [],
            "organism_wall_s": round(time.time() - org_t0, 1),
        }

    print(f"[{organism}] building production target sets ...", flush=True)
    tiers, census = _load_target_sets(organism)
    print(f"[{organism}] per-tier (raw symbols): {json.dumps(census.per_tier)}", flush=True)

    out_dir = out_root / organism
    # Long-form significant-subset lives in a sibling tree so the matrix tree is
    # untouched: data/processed/enrichment_significant/<org>/<FAMILY>_<TIER>.tsv.
    long_dir = out_root.parent / "enrichment_significant" / organism
    cells: list[dict[str, object]] = []
    family_blocks: list[dict[str, object]] = []
    sanity: list[dict[str, object]] = []

    for family in families:
        print(f"[{organism}] materialising {family} TERM2GENE ...", flush=True)
        t2g, universe, prov = _build_family_sets(organism, family, rscript=rscript)
        fam_symbols = set(universe)

        # Re-case every band onto THIS family's symbol space (TargetScan case fix).
        norm_tiers, norm_stats = pts.normalize_tiers_to_go_symbols(tiers, fam_symbols)

        family_blocks.append(
            {
                "family": family,
                "n_terms_total": len(t2g),
                "n_universe": len(universe),
                "provenance": prov,
                "symbol_normalization": norm_stats,
            }
        )

        # B2: a (family, organism) cell with an empty surface would yield an
        # all-p=1.0 matrix; refuse it before scoring any tier.
        if not t2g or not universe:
            raise ValueError(
                f"[{organism}] {family}: empty TERM2GENE/universe -- refusing degenerate matrix"
            )
        # Full-family size M is fixed by the inputs; with universe = union of the
        # memberships every term is eligible, so anything else is a build defect.
        m_family, _ = family_bh.eligible_family(t2g, universe)
        if m_family != len(t2g):
            raise ValueError(
                f"[{organism}] {family}: eligible family M={m_family} != "
                f"{len(t2g)} TERM2GENE terms -- refusing an inconsistent family"
            )
        for tier in tiers_wanted:
            gene_sets = pts.cumulative_tier(norm_tiers, tier)
            gene_sets = {m: g for m, g in gene_sets.items() if g}
            if not gene_sets:
                if organism == "mml" and tier == "MEDIUM":
                    print(
                        f"[{organism}] SKIP {family}/{tier}: documented-real macaque "
                        f"MEDIUM=0 (lone validated MTI is in HIGH)",
                        flush=True,
                    )
                    continue
                raise ValueError(
                    f"[{organism}] {family}/{tier}: no miRNAs with targets "
                    f"(expected non-empty) -- refusing a degenerate matrix."
                )
            t0 = time.time()
            res, cp_prov = cp.run_ora_batch_with_provenance(
                gene_sets, t2g, universe, term2name=None, rscript=rscript
            )
            elapsed = time.time() - t0
            scored_terms = sorted({term for per in res.values() for term in per})
            # B2: gene sets exist but nothing scored -> refuse the all-p=1.0 matrix.
            if gene_sets and not scored_terms:
                raise ValueError(
                    f"[{organism}] {family}/{tier}: {len(gene_sets)} gene set(s) but "
                    f"0 terms scored -- refusing degenerate all-p=1.0 matrix "
                    f"(likely an enricher/worker failure)."
                )
            mirnas = sorted(gene_sets)
            out_path = out_dir / f"prod_{family}_{tier}.csv"
            n_rows, n_cols = _write_matrix(out_path, mirnas, scored_terms, res)
            cell = {
                "organism": organism,
                "family": family,
                "tier": tier,
                "n_mirnas": n_rows,
                "n_pathways_scored": n_cols,
                "n_pathways_total": len(t2g),
                "n_universe": len(universe),
                "m_family": m_family,
                "elapsed_s": round(elapsed, 2),
                "cp_provenance": cp_prov,
                "output": str(out_path.relative_to(out_root.parents[2])),
            }
            if long_form:
                long_path = long_dir / f"{family}_{tier}.tsv"
                n_sig = write_significant_long_form(long_path, mirnas, res, m_family=m_family)
                counts = count_retained_rows(res, fdr_cutoff=FDR_CUTOFF, m_family=m_family)
                cell["long_form_output"] = str(long_path.relative_to(out_root.parents[2]))
                cell["adjustment_family"] = "full_family_bh"
                cell["n_retained_rows"] = n_sig
                # Alias read by validation/go_propagation_driver.py.
                cell["n_significant_rows"] = n_sig
                cell["n_retained_rows_overlap_legacy"] = counts["rows_overlap_legacy"]
                cell["m_overlap_min"] = counts["m_overlap_min"]
                cell["m_overlap_max"] = counts["m_overlap_max"]
                print(
                    f"[{organism}] {family}/{tier}: {n_sig} retained rows "
                    f"(full-family M={m_family}; overlap-only would retain "
                    f"{counts['rows_overlap_legacy']}) -> {long_path.name}",
                    flush=True,
                )
            cells.append(cell)
            print(
                f"[{organism}] {family}/{tier}: {n_rows} miRNAs x {n_cols} pathways "
                f"(universe={len(universe)}) {elapsed:.1f}s -> {out_path.name}",
                flush=True,
            )

            # Sanity hit: probe miRNA -> its lowest-p pathways in this cell.
            probe = SANITY_MIRNA.get(organism)
            if probe and probe in res and res[probe]:
                top = sorted(res[probe].items(), key=lambda kv: kv[1][0])[:8]
                sanity.append(
                    {
                        "family": family,
                        "tier": tier,
                        "mirna": probe,
                        "n_targets": len(gene_sets[probe]),
                        "top_pathways": [
                            {"pathway_id": t, "pvalue": p, "p_adjust": padj, "overlap": ov}
                            for t, (p, padj, ov) in top
                        ],
                    }
                )

    return {
        "organism": organism,
        "taxid": taxid,
        "org_db_package": gg.org_db_for_taxid(taxid),
        "families": family_blocks,
        "per_source_census": census.per_source,
        "per_tier_raw": census.per_tier,
        "census_notes": census.notes,
        "cells": cells,
        "sanity": sanity,
        "organism_wall_s": round(time.time() - org_t0, 1),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organisms", nargs="+", default=["hsa", "mmu", "rno", "mml"])
    parser.add_argument("--families", nargs="+", default=list(DEFAULT_FAMILIES))
    parser.add_argument("--tiers", nargs="+", default=list(DEFAULT_TIERS))
    parser.add_argument("--rscript", default=None)
    parser.add_argument(
        "--no-long-form",
        dest="long_form",
        action="store_false",
        default=True,
        help=(
            "Do NOT write the authoritative significant-subset long-form TSV per "
            "cell (data/processed/enrichment_significant/<org>/<FAMILY>_<TIER>.tsv). "
            "By DEFAULT it IS written (Deliverable B, contract v2) with mirna, "
            "pathway_id, p_value, p_adjust, overlap, p_adjust_overlap, m_overlap, "
            "m_family -- p_adjust is the full-family BH value over the cell's M "
            "eligible terms, p_adjust_overlap the enricher's overlap-only value; "
            "raw p_value + overlap are the authoritative run_ora_batch tuples. "
            "The matrix output is unchanged."
        ),
    )
    args = parser.parse_args(argv)

    for org in args.organisms:
        if org not in ORG_TAXID:
            raise SystemExit(f"unsupported organism {org!r}; known: {sorted(ORG_TAXID)}")
    for fam in args.families:
        if fam not in FAMILY_TOKENS:
            raise SystemExit(f"unknown family token {fam!r}; choose from {FAMILY_TOKENS}")
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
        block = run_organism(
            organism,
            families_wanted=args.families,
            tiers_wanted=args.tiers,
            out_root=out_root,
            rscript=args.rscript,
            long_form=args.long_form,
        )
        blocks.append(block)
        print(f"[{organism}] done in {block.get('organism_wall_s')}s", flush=True)

    # Merge into any existing production-pathways manifest: organisms run THIS
    # invocation replace their prior block (others are kept).
    manifest_path = out_root / "production_pathways_manifest.json"
    run_orgs = set(args.organisms)
    kept: list[dict[str, object]] = []
    if manifest_path.exists():
        try:
            prior = json.loads(manifest_path.read_text(encoding="utf-8"))
            kept = [b for b in prior.get("organisms", []) if b.get("organism") not in run_orgs]
        except (OSError, json.JSONDecodeError):
            kept = []
    manifest = {
        "build": "production_pathways_ora",
        "note": (
            "Uniform production ORA on REAL gene sets (gene-SYMBOL space) across "
            "Reactome + WikiPathways + GO. Matrices: prod_<FAMILY>_<TIER>.csv "
            "(distinct from v2-comparison hsa/mmu + existing production rno/mml GO). "
            "Reactome/WP exclude macaque (no coverage); GO only adds hsa/mmu "
            "(rno/mml production GO already exist)."
        ),
        "matrix_naming": "data/processed/enrichment/<org>/prod_<FAMILY>_<TIER>.csv",
        "last_run_organisms": sorted(run_orgs),
        "last_run_families": args.families,
        "last_run_tiers": args.tiers,
        "last_run_wall_s": round(time.time() - t0, 1),
        "organisms": kept + blocks,
    }
    payload = json.dumps(manifest, indent=2)
    manifest_path.write_text(payload, encoding="utf-8")

    # Mirror the manifest into the tracked provenance/ dir (data/processed/* is
    # gitignored, so this committed copy is the auditable record of the run --
    # same convention as provenance/production_enrichment_manifest.json from the
    # GO run). The matrices themselves stay gitignored; this is the record.
    prov_dir = cfg.path("provenance_dir")
    prov_dir.mkdir(parents=True, exist_ok=True)
    prov_copy = prov_dir / "production_pathways_manifest.json"
    prov_copy.write_text(payload, encoding="utf-8")

    print(
        f"\nTOTAL wall: {time.time() - t0:.1f}s; manifest -> {manifest_path}\n"
        f"committed copy -> {prov_copy}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
