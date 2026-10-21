"""GO ``all`` vs ``direct`` propagation sensitivity -- driver (IO / R edges).

Orchestrates the comparison defined in :mod:`mirpathdb_prod.validation.go_propagation`
WITHOUT touching the production default (``config.yaml`` ``go_propagation: all``)
and WITHOUT writing under ``data/processed/``: every generated file lives under
the comparison root (default ``<validation_dir>/go_propagation_sensitivity/``).

Flow
----
1. Provenance sha256 gate (``verify_cached_sources``) -- the same block the
   production runners run; the consumed target sources' sha256 are recorded.
2. Target sets from the PRODUCTION loaders (``scripts/run_production_pathways.py``
   ``_load_target_sets`` for hsa/mmu, ``scripts/run_production_enrichment.py``
   ``_load_sources`` for rno/mml, loaded by file path because ``scripts/`` is not
   a package). Bands are re-cased exactly as the runners do: hsa/mmu PER FAMILY
   universe, rno/mml ONCE onto the UNION of the three GO family universes.
3. GO gene sets for BOTH modes (``build_go_genesets(taxid, propagation=mode)``),
   persisted as ``<out>/<mode>/<org>/go_term2gene.tsv`` + ``go_provenance.json``
   and reused when the recorded versions match the R environment.
4. GO depths: a GO.db ``GO{BP,CC,MF}PARENTS`` edge dump
   (``<out>/go_parents_edges.tsv``, reused for the same GO.db version) + a Python
   BFS.
5. Structural comparison per organism x family.
6. Enrichment (unless ``--structural-only``): the ``all`` side is the production
   output when its fingerprint matches. When it does not, the run is REFUSED
   (exit 2, every mismatch listed) unless ``--recompute-all`` is given; with the
   flag a fingerprint-compatible isolated run under ``<out>/all/`` is reused,
   else one is computed there. The ``direct`` side is computed under
   ``<out>/direct/`` unless a fingerprint-compatible run is already there.
   Isolated sides write the production-format matrix, the FULL long-form cells
   TSV (every tested cell with exact p / p_adjust / overlap) and a manifest after
   EACH organism, so a crash leaves reusable per-organism state. A reused
   production cell overlays BOTH the shipped significant TSV and the
   ``species_excluded`` quarantine (the contract export moves foreign-species
   rows there; the matrices keep every miRNA); the ``> 0.25`` padj floor is
   decided PER miRNA from the data (without the quarantine only the miRNAs the
   species policy retains are floored) and dropped when the row files fail the
   integrity checks in :func:`_load_production_cell_table`.
7. Deterministic artefacts + ``manifest.yaml`` (the only file with timestamps;
   ``release.contract_report_generated_at`` is the release identifier copied
   from ``provenance/frontend_contract_report.json``, not this run's time).
   ``manifest.yaml`` is written with ``status: running`` before any R work and
   rewritten as ``ok`` / ``refused`` / ``failed`` (with the error) at the end, so
   an out-dir never carries an earlier run's ``ok`` next to a failed run's files.

Side status vocabulary (``manifest.yaml`` ``sides.<mode>.status``):
``reused_production`` | ``reused_isolated`` | ``computed`` | ``skipped`` |
``refused`` | ``mixed`` (organisms differ; see ``status_by_organism``).
Every path inside the summary JSON / report is relative to the config root
(the directory of the loaded ``config.yaml``), never absolute.

``--rscript`` is exported as ``MPD_RSCRIPT`` for the whole run so the production
target-set loaders (whose RefSeq -> SYMBOL mapping shells its own Rscript) use
the same binary as the GO builds and the ORA; ``--max-cores`` is exported as
``MPD_MAX_CORES`` the same way.

Exit codes: 0 success; 2 refusal (incompatible ``all`` without ``--recompute-all``,
anything missing under ``--reuse-only``); hard errors propagate to the CLI (1).
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import os
import subprocess
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from types import ModuleType
from typing import Any

import yaml

from mirpathdb_prod.config import Config
from mirpathdb_prod.enrichment import clusterprofiler as cp
from mirpathdb_prod.pathways import go_genesets as gg
from mirpathdb_prod.sources import production_target_sets as pts
from mirpathdb_prod.sources import tarbase
from mirpathdb_prod.validation import go_propagation as gp
from mirpathdb_prod.validation.provenance import file_sha256, verify_cached_sources

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS_DIR = REPO_ROOT / "scripts"

ORG_TAXID: dict[str, int] = {"hsa": 9606, "mmu": 10090, "rno": 10116, "mml": 9544}
DEFAULT_ORGANISMS: tuple[str, ...] = ("hsa", "mmu", "rno", "mml")
DEFAULT_TIERS: tuple[str, ...] = pts.TIERS
DEFAULT_FAMILIES: tuple[str, ...] = gg.GO_FAMILIES
DEFAULT_OUT_SUBDIR = "go_propagation_sensitivity"

# How each production runner re-cases the target bands onto GO symbol space.
NORMALIZATION_PER_FAMILY = "per_family"
NORMALIZATION_UNION = "union"
NORMALIZATION_MODE: dict[str, str] = {
    "hsa": NORMALIZATION_PER_FAMILY,
    "mmu": NORMALIZATION_PER_FAMILY,
    "rno": NORMALIZATION_UNION,
    "mml": NORMALIZATION_UNION,
}
# Production matrix naming: hsa/mmu GO came from run_production_pathways (prod_
# prefix), rno/mml from run_production_enrichment (unprefixed).
PRODUCTION_MATRIX_PREFIX: dict[str, str] = {"hsa": "prod_", "mmu": "prod_", "rno": "", "mml": ""}

MODES: tuple[str, ...] = ("all", "direct")
PRODUCTION_MODE = "all"

# Cached raw sources each organism's production loader consumes (sha256-recorded).
_TS = "Predicted_Targets_Context_Scores.default_predictions.txt.zip"
_MIRDB = "miRDB_v6.0_prediction_result.txt.gz"
TARGET_SOURCE_FILES: dict[str, tuple[str, ...]] = {
    "hsa": ("mirtarbase_hsa_MTI.xlsx", "Homo_sapiens_TarBase-v9.tsv.gz", _TS, _MIRDB),
    "mmu": ("mirtarbase_mmu_MTI.xlsx", "Mus_musculus_TarBase-v9.tsv.gz", _TS, _MIRDB),
    "rno": ("mirtarbase_rno_MTI.xlsx", tarbase.OTHER_SPECIES_FILENAME, _TS, _MIRDB),
    "mml": (tarbase.OTHER_SPECIES_FILENAME, _TS),
}

SIDE_BUILD = "go_propagation_sensitivity_side"
SIDE_MANIFEST = "manifest.json"
GO_T2G_TSV = "go_term2gene.tsv"
GO_PROV_JSON = "go_provenance.json"
GO_EDGES_TSV = "go_parents_edges.tsv"
GO_EDGES_PROV = "go_parents_edges.provenance.json"

STATUS_REUSED_PRODUCTION = "reused_production"
STATUS_REUSED_ISOLATED = "reused_isolated"
STATUS_COMPUTED = "computed"
STATUS_SKIPPED = "skipped"
STATUS_REFUSED = "refused"
STATUS_MIXED = "mixed"  # aggregate only: organisms differ (see status_by_organism)

SPECIES_EXCLUDED_DIRNAME = "species_excluded"  # export_frontend_contracts.EXCLUDED_DIRNAME
PRODUCTION_SOURCE_SHA_NOTE = (
    "target-source sha256 not recorded by the production manifests; the production all "
    "side is accepted on the provenance/manifest.yaml sha256 gate (manifest.yaml "
    "provenance_gate) plus the per_tier_raw census, universes, term counts and versions"
)

SOURCE_PRODUCTION = "production"
SOURCE_ISOLATED = "isolated"

LOG_PREFIX = "[go-sens]"


def _log(msg: str) -> None:
    print(f"{LOG_PREFIX} {msg}", flush=True)


@dataclass
class ComparisonOptions:
    """CLI-facing options of one comparison run."""

    organisms: tuple[str, ...] = DEFAULT_ORGANISMS
    tiers: tuple[str, ...] = DEFAULT_TIERS
    families: tuple[str, ...] = DEFAULT_FAMILIES
    out_dir: Path | None = None
    rscript: str | None = None
    fdr: float | None = None
    secondary_fdrs: tuple[float, ...] | None = None
    structural_only: bool = False
    recompute_all: bool = False
    reuse_only: bool = False
    max_cores: int | None = None
    argv: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# R edges (thin; monkeypatched in tests)
# --------------------------------------------------------------------------- #
_R_JSON_HELPER = r"""
manual_json <- function(x) {
  parts <- vapply(names(x), function(k)
    sprintf('"%s": "%s"', k, gsub('"', '\\"', as.character(x[[k]]))), character(1))
  paste0("{", paste(parts, collapse = ", "), "}")
}
"""

# Dumps every (child, parent, relationship, ontology) edge of GO.db's
# GO{BP,CC,MF}PARENTS maps (all relationship types), minus the 'all' pseudo-root.
_R_GO_PARENTS_DRIVER = r"""
suppressMessages(library(GO.db))
args      <- commandArgs(trailingOnly = TRUE)
out_file  <- args[[1]]
prov_file <- args[[2]]

dump_one <- function(map, ont) {
  lst <- as.list(map)
  n   <- lengths(lst)
  lst <- lst[n > 0]
  n   <- n[n > 0]
  child  <- rep(names(lst), n)
  parent <- unname(unlist(lst, use.names = FALSE))
  rel <- unlist(lapply(lst, function(v) {
    r <- names(v)
    if (is.null(r)) rep(NA_character_, length(v)) else r
  }), use.names = FALSE)
  data.frame(child = child, parent = parent, relationship = rel,
             ontology = ont, stringsAsFactors = FALSE)
}
edges <- rbind(dump_one(GOBPPARENTS, "BP"),
               dump_one(GOCCPARENTS, "CC"),
               dump_one(GOMFPARENTS, "MF"))
edges <- edges[edges$parent != "all", ]
edges <- edges[order(edges$ontology, edges$child, edges$parent, edges$relationship), ]
write.table(edges, file = out_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
prov <- list(
  go_db_version = as.character(packageVersion("GO.db")),
  n_edges       = as.character(nrow(edges))
)
writeLines(manual_json(prov), prov_file)
"""

# Prints key<TAB>value lines: R, clusterProfiler, GO.db and each org.db version.
_R_VERSION_PROBE = r"""
args <- commandArgs(trailingOnly = TRUE)
cat("r_version\t", paste(R.version$major, R.version$minor, sep = "."), "\n", sep = "")
cat("clusterprofiler_version\t", as.character(packageVersion("clusterProfiler")), "\n", sep = "")
cat("go_db_version\t", as.character(packageVersion("GO.db")), "\n", sep = "")
for (pkg in args) cat(pkg, "\t", as.character(packageVersion(pkg)), "\n", sep = "")
"""


def _run_rscript(rscript: str | None, script: str, args: Sequence[str], *, label: str) -> str:
    """Run an inline R script; return stdout; raise RuntimeError with stderr on failure."""
    import tempfile

    rscript_bin = gg._resolve_rscript(rscript)
    with tempfile.TemporaryDirectory(prefix="mpd_go_sens_") as tmp:
        driver = Path(tmp) / "driver.R"
        driver.write_text(script, encoding="utf-8")
        proc = subprocess.run(
            [rscript_bin, str(driver), *args], capture_output=True, text=True, check=False
        )
    if proc.returncode != 0:
        raise RuntimeError(
            f"{label} Rscript failed (exit {proc.returncode}).\n"
            f"--- stderr ---\n{proc.stderr}\n--- stdout ---\n{proc.stdout}"
        )
    return proc.stdout


def _probe_r_versions(rscript: str | None, org_db_packages: Sequence[str]) -> dict[str, str]:
    """Return ``{r_version, clusterprofiler_version, go_db_version, <org.db>: version}``."""
    out = _run_rscript(rscript, _R_VERSION_PROBE, list(org_db_packages), label="version probe")
    versions: dict[str, str] = {}
    for line in out.splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            versions[k.strip()] = v.strip()
    for key in ("r_version", "clusterprofiler_version", "go_db_version"):
        if key not in versions:
            raise RuntimeError(f"version probe did not report {key}: {out!r}")
    return versions


def _dump_go_parent_edges(rscript: str | None, out_path: Path, prov_path: Path) -> None:
    """Write the GO.db parent-edge TSV + its provenance JSON (see ``_R_GO_PARENTS_DRIVER``)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    _run_rscript(
        rscript,
        _R_JSON_HELPER + _R_GO_PARENTS_DRIVER,
        [str(out_path), str(prov_path)],
        label="GO.db parent edges",
    )


# --------------------------------------------------------------------------- #
# Production loaders (scripts/ is not a package -> load by file path)
# --------------------------------------------------------------------------- #
_SCRIPT_CACHE: dict[str, ModuleType] = {}


def _load_script(name: str) -> ModuleType:
    mod = _SCRIPT_CACHE.get(name)
    if mod is None:
        path = SCRIPTS_DIR / f"{name}.py"
        spec = importlib.util.spec_from_file_location(name, path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"cannot load {path}")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _SCRIPT_CACHE[name] = mod
    return mod


def _load_target_sets(
    organism: str, cfg: Config
) -> tuple[dict[str, dict[str, set[str]]], pts.TargetSetCensus]:
    """Production target sets via the runners' own loaders (raw symbol space).

    hsa/mmu: ``run_production_pathways._load_target_sets``; rno/mml:
    ``run_production_enrichment._load_sources``. The runner's module-level
    ``CACHE_DIR`` is pointed at this config's ``source_cache_dir`` so ``-c`` is
    honoured.
    """
    if organism in ("hsa", "mmu"):
        runner = _load_script("run_production_pathways")
        runner.CACHE_DIR = cfg.path("source_cache_dir")
        return runner._load_target_sets(organism)
    runner = _load_script("run_production_enrichment")
    runner.CACHE_DIR = cfg.path("source_cache_dir")
    return runner._load_sources(organism)


def _source_sha256_by_organism(
    cache_dir: Path, organisms: Sequence[str]
) -> dict[str, dict[str, str]]:
    """sha256 of every consumed cached source that is present (absent -> omitted)."""
    out: dict[str, dict[str, str]] = {}
    digests: dict[str, str] = {}
    for org in organisms:
        block: dict[str, str] = {}
        for base in TARGET_SOURCE_FILES.get(org, ()):
            path = cache_dir / base
            if not path.exists():
                continue
            if base not in digests:
                digests[base] = file_sha256(path)
            block[base] = digests[base]
        out[org] = block
    return out


# --------------------------------------------------------------------------- #
# GO gene sets + depth edges (cached under the comparison root)
# --------------------------------------------------------------------------- #
def _ensure_go_genesets(
    out_dir: Path,
    mode: str,
    organism: str,
    *,
    rscript: str | None,
    versions: Mapping[str, str],
) -> gg.GoGeneSets:
    """Build (or reload) the GO gene sets of one organism under one propagation mode.

    The materialisation is persisted as ``<out>/<mode>/<org>/go_term2gene.tsv`` +
    ``go_provenance.json`` and reused when ``org_db_version``, ``go_db_version``,
    ``r_version``, ``go_columns`` match the probed R environment and ``n_pairs``
    matches the reloaded TSV.
    """
    if mode not in gg.GO_PROPAGATION_COLUMNS:
        raise ValueError(
            f"unknown propagation {mode!r}; choose from {sorted(gg.GO_PROPAGATION_COLUMNS)}"
        )
    taxid = ORG_TAXID[organism]
    pkg = gg.org_db_for_taxid(taxid)
    go_col, ont_col = gg.GO_PROPAGATION_COLUMNS[mode]
    org_dir = out_dir / mode / organism
    tsv = org_dir / GO_T2G_TSV
    prov_path = org_dir / GO_PROV_JSON
    if tsv.exists() and prov_path.exists():
        try:
            prov = json.loads(prov_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            prov = {}
        expected = {
            "go_propagation": mode,
            "go_columns": f"{go_col}/{ont_col}",
            "org_db_version": versions.get(pkg),
            "go_db_version": versions.get("go_db_version"),
            "r_version": versions.get("r_version"),
        }
        if all(str(prov.get(k)) == str(v) for k, v in expected.items()):
            go = gg.read_go_genesets_tsv(tsv, taxid=taxid, org_db_package=pkg)
            n_pairs = sum(len(s) for fam in go.term2gene.values() for s in fam.values())
            if str(prov.get("n_pairs", n_pairs)) == str(n_pairs):
                go.provenance = {k: str(v) for k, v in prov.items()}
                _log(f"{organism}/{mode}: reusing cached GO gene sets ({n_pairs} pairs)")
                return go
        _log(f"{organism}/{mode}: cached GO gene sets stale (provenance mismatch); rebuilding")
    _log(f"{organism}/{mode}: building GO gene sets from {pkg} ({go_col}/{ont_col}) ...")
    go = gg.build_go_genesets(taxid, rscript=rscript, propagation=mode)
    gg.write_go_genesets_tsv(go, tsv)
    n_pairs = sum(len(s) for fam in go.term2gene.values() for s in fam.values())
    prov = {k: str(v) for k, v in go.provenance.items()}
    prov.setdefault("n_pairs", str(n_pairs))
    prov.setdefault("org_db_package", pkg)
    prov["go_propagation"] = mode
    prov_path.write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    go.provenance = prov
    return go


def _ensure_go_parent_edges(
    out_dir: Path, *, rscript: str | None, go_db_version: str | None
) -> tuple[Path, dict[str, str]]:
    """Dump (or reuse) the GO.db parent edges; reuse requires the same GO.db version."""
    edges = out_dir / GO_EDGES_TSV
    prov_path = out_dir / GO_EDGES_PROV
    if edges.exists() and prov_path.exists():
        try:
            prov = json.loads(prov_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            prov = {}
        if go_db_version is None or str(prov.get("go_db_version")) == str(go_db_version):
            _log(f"reusing GO.db parent edges ({prov.get('n_edges')} edges)")
            return edges, {k: str(v) for k, v in prov.items()}
        _log("cached GO.db parent edges are for another GO.db version; re-dumping")
    _log("dumping GO.db parent edges (GO{BP,CC,MF}PARENTS, all relationship types) ...")
    _dump_go_parent_edges(rscript, edges, prov_path)
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    return edges, {k: str(v) for k, v in prov.items()}


# --------------------------------------------------------------------------- #
# Target-set normalisation + expected cells (mirror the runners EXACTLY)
# --------------------------------------------------------------------------- #
NormByFamily = dict[str, tuple[dict[str, dict[str, set[str]]], dict[str, dict[str, int]]]]


def _normalize_by_family(
    organism: str,
    tiers: dict[str, dict[str, set[str]]],
    go: gg.GoGeneSets,
    families: Sequence[str],
) -> NormByFamily:
    """Re-case the bands the way the organism's production runner does.

    hsa/mmu (``run_production_pathways``): per family, onto THAT family's universe.
    rno/mml (``run_production_enrichment``): once, onto the union of the three GO
    family universes (regardless of which families are compared).
    """
    out: NormByFamily = {}
    if NORMALIZATION_MODE[organism] == NORMALIZATION_PER_FAMILY:
        for fam in families:
            out[fam] = pts.normalize_tiers_to_go_symbols(tiers, go.universe(fam))
        return out
    go_symbols: set[str] = set()
    for fam in gg.GO_FAMILIES:
        go_symbols |= go.universe(fam)
    shared = pts.normalize_tiers_to_go_symbols(tiers, go_symbols)
    for fam in families:
        out[fam] = shared
    return out


def _gene_sets_for(
    organism: str, norm_tiers: dict[str, dict[str, set[str]]], tier: str
) -> dict[str, set[str]] | None:
    """Cumulative gene sets of one tier (empty miRNAs dropped); ``None`` = skipped cell.

    Mirrors the runners: an empty cumulative tier is only tolerated for the
    documented macaque MEDIUM case; any other empty tier is a silent-failure
    signal and raises.
    """
    gene_sets = {m: g for m, g in pts.cumulative_tier(norm_tiers, tier).items() if g}
    if gene_sets:
        return gene_sets
    if organism == "mml" and tier == "MEDIUM":
        return None
    raise ValueError(
        f"[{organism}] tier {tier}: no miRNAs with targets after rollup "
        f"(expected non-empty) -- refusing a degenerate matrix."
    )


def _expected_cells(
    organism: str, norm: NormByFamily, tiers: Sequence[str], families: Sequence[str]
) -> tuple[tuple[str, str], ...]:
    cells: list[tuple[str, str]] = []
    for fam in families:
        norm_tiers, _ = norm[fam]
        for tier in tiers:
            if _gene_sets_for(organism, norm_tiers, tier) is not None:
                cells.append((fam, tier))
    return tuple(sorted(cells))


# --------------------------------------------------------------------------- #
# Isolated side store (<out>/<mode>/manifest.json)
# --------------------------------------------------------------------------- #
class SideStore:
    """Per-mode isolated output tree + its manifest (written after each organism)."""

    def __init__(self, root: Path, mode: str) -> None:
        self.root = root
        self.mode = mode
        self.path = root / SIDE_MANIFEST
        self.doc: dict[str, Any] = {
            "build": SIDE_BUILD,
            "go_propagation": mode,
            "note": (
                "Isolated GO ORA for the propagation sensitivity comparison; NOT a "
                "production output. Matrices use the production _write_matrix format; "
                "<FAMILY>_<TIER>.cells.tsv.gz holds EVERY tested cell."
            ),
            "fingerprint": None,
            "organisms": [],
        }
        if self.path.exists():
            try:
                self.doc = json.loads(self.path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                pass

    def fingerprint(self) -> gp.RunFingerprint | None:
        fp = self.doc.get("fingerprint")
        return gp.RunFingerprint.from_dict(fp) if isinstance(fp, Mapping) else None

    def block(self, organism: str) -> dict[str, Any] | None:
        for b in self.doc.get("organisms", []) or []:
            if b.get("organism") == organism:
                return b
        return None

    def put(self, organism: str, block: dict[str, Any], fp_subset: gp.RunFingerprint) -> None:
        kept = [b for b in self.doc.get("organisms", []) or [] if b.get("organism") != organism]
        self.doc["organisms"] = [*kept, block]
        current = self.fingerprint()
        merged = fp_subset if current is None else current.merged_with(fp_subset)
        self.doc["fingerprint"] = merged.to_dict()
        self.write()

    def write(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.doc, indent=2) + "\n", encoding="utf-8")

    def matrix_path(self, organism: str, family: str, tier: str) -> Path:
        return self.root / organism / f"{family}_{tier}.csv"

    def cells_path(self, organism: str, family: str, tier: str) -> Path:
        return self.root / organism / f"{family}_{tier}.cells.tsv.gz"


# --------------------------------------------------------------------------- #
# Computing one side for one organism (runner logic, isolated destination)
# --------------------------------------------------------------------------- #
def _compute_side_organism(
    store: SideStore,
    organism: str,
    *,
    go: gg.GoGeneSets,
    norm: NormByFamily,
    census: pts.TargetSetCensus,
    tiers: Sequence[str],
    families: Sequence[str],
    rscript: str | None,
    argv: Sequence[str],
) -> dict[str, Any]:
    """Run every (family x tier) ORA cell of one organism into the side store."""
    taxid = ORG_TAXID[organism]
    t_org = time.time()
    cells: list[dict[str, Any]] = []
    family_blocks: list[dict[str, Any]] = []
    for family in families:
        t2g = go.family(family)
        universe = go.universe(family)
        if not t2g or not universe:
            raise ValueError(
                f"[{organism}] {family}: empty TERM2GENE/universe -- refusing degenerate matrix"
            )
        norm_tiers, norm_stats = norm[family]
        family_blocks.append(
            {
                "family": family,
                "n_terms_total": len(t2g),
                "n_universe": len(universe),
                "symbol_normalization": norm_stats,
            }
        )
        for tier in tiers:
            gene_sets = _gene_sets_for(organism, norm_tiers, tier)
            if gene_sets is None:
                _log(
                    f"{organism}/{store.mode}: SKIP {family}/{tier}: documented-real macaque "
                    f"MEDIUM=0 (lone validated MTI is in HIGH)"
                )
                continue
            t0 = time.time()
            res, cp_prov = cp.run_ora_batch_with_provenance(
                gene_sets, t2g, universe, term2name=None, rscript=rscript
            )
            elapsed = time.time() - t0
            scored_terms = sorted({term for per in res.values() for term in per})
            if gene_sets and not scored_terms:
                raise ValueError(
                    f"[{organism}] {family}/{tier}: {len(gene_sets)} gene set(s) but 0 terms "
                    f"scored -- refusing degenerate all-p=1.0 matrix (likely an "
                    f"enricher/worker failure)."
                )
            mirnas = sorted(gene_sets)
            matrix = store.matrix_path(organism, family, tier)
            n_rows, n_cols = gp.write_matrix(matrix, mirnas, scored_terms, res)
            cells_path = store.cells_path(organism, family, tier)
            n_tested = gp.write_cells_long_form(cells_path, mirnas, res)
            cells.append(
                {
                    "organism": organism,
                    "family": family,
                    "tier": tier,
                    "n_mirnas": n_rows,
                    "n_terms_scored": n_cols,
                    "n_terms_total": len(t2g),
                    "n_universe": len(universe),
                    "n_tested_cells": n_tested,
                    "elapsed_s": round(elapsed, 2),
                    "cp_provenance": dict(cp_prov),
                    "matrix": str(matrix.relative_to(store.root)),
                    "cells_output": str(cells_path.relative_to(store.root)),
                }
            )
            _log(
                f"{organism}/{store.mode}: {family}/{tier}: {n_rows} miRNAs x {n_cols} terms "
                f"({n_tested} tested cells, universe={len(universe)}) {elapsed:.1f}s"
            )
    return {
        "organism": organism,
        "taxid": taxid,
        "org_db_package": go.org_db_package,
        "go_propagation": store.mode,
        "go_provenance": dict(go.provenance),
        "normalization_mode": NORMALIZATION_MODE[organism],
        "per_source_census": census.per_source,
        "per_tier_raw": census.per_tier,
        "census_notes": census.notes,
        "families": family_blocks,
        "cells": cells,
        "commands": [list(argv)],
        "organism_wall_s": round(time.time() - t_org, 1),
    }


# --------------------------------------------------------------------------- #
# Production ``all`` side
# --------------------------------------------------------------------------- #
def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return doc if isinstance(doc, dict) else None


def _rel(path: Path, root: Path) -> str:
    """``path`` relative to the config root (the deterministic artefacts must not
    embed the absolute checkout / data location); absolute only when outside it."""
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def _config_root(cfg: Config) -> Path:
    """Directory the config's relative paths resolve against (``-c`` honoured)."""
    return Path(cfg.source_path).resolve().parent if cfg.source_path else REPO_ROOT


@dataclass
class ProductionSide:
    """Where the production ``all`` outputs + manifests live for this config.

    ``contract_cells`` is the ``emitted_contract.cells[]`` inventory of
    ``provenance/production_enrichment_manifest.json`` keyed by
    ``(organism, family, tier)`` -- it records the sha256 of each SHIPPED
    significant TSV -- or ``None`` when no contract export is recorded. ``root``
    is the config root every recorded path is made relative to.
    """

    root: Path
    enrichment_dir: Path
    significant_dir: Path
    species_excluded_dir: Path
    pathways_manifest_path: Path | None
    enrichment_manifest_path: Path | None
    pathways_manifest: dict[str, Any] | None
    enrichment_manifest: dict[str, Any] | None
    contract_cells: dict[tuple[str, str, str], dict[str, Any]] | None = None

    def rel(self, path: Path) -> str:
        return _rel(path, self.root)

    def matrix_path(self, organism: str, family: str, tier: str) -> Path:
        prefix = PRODUCTION_MATRIX_PREFIX[organism]
        return self.enrichment_dir / organism / f"{prefix}{family}_{tier}.csv"

    def significant_path(self, organism: str, family: str, tier: str) -> Path | None:
        p = self.significant_dir / organism / f"{family}_{tier}.tsv"
        return p if p.exists() else None

    def species_excluded_path(self, organism: str, family: str, tier: str) -> Path | None:
        p = self.species_excluded_dir / organism / f"{family}_{tier}.tsv"
        return p if p.exists() else None

    def contract_cell(self, organism: str, family: str, tier: str) -> dict[str, Any] | None:
        """The emitted-contract record of one shipped cell (None = no record)."""
        if self.contract_cells is None:
            return None
        return self.contract_cells.get((organism, family, tier))

    def build_manifest(self, organism: str) -> dict[str, Any] | None:
        doc = self.pathways_manifest if organism in ("hsa", "mmu") else self.enrichment_manifest
        return doc if isinstance(doc, dict) else None

    def build_cell(self, organism: str, family: str, tier: str) -> dict[str, Any] | None:
        """The build manifest's ``cells[]`` block of one production cell, or None."""
        doc = self.build_manifest(organism) or {}
        for blk in doc.get("organisms", []) or []:
            if blk.get("organism") != organism:
                continue
            for c in blk.get("cells", []) or []:
                if c.get("family") == family and c.get("tier") == tier:
                    return c
        return None

    def n_significant_rows_recorded(self, organism: str, family: str, tier: str) -> int | None:
        """``n_significant_rows`` the build wrote for the cell = shipped rows +
        quarantined rows (the export moves rows, it never drops them)."""
        cell = self.build_cell(organism, family, tier)
        raw = (cell or {}).get("n_significant_rows")
        if raw is None:
            return None
        try:
            return int(raw)
        except (TypeError, ValueError):
            return None

    def manifest_for(self, organism: str) -> Path | None:
        if organism in ("hsa", "mmu"):
            return self.pathways_manifest_path
        return self.enrichment_manifest_path


def _contract_cells(provenance_dir: Path) -> dict[tuple[str, str, str], dict[str, Any]] | None:
    """``emitted_contract.cells[]`` of the committed enrichment manifest, or None."""
    doc = _read_json(provenance_dir / "production_enrichment_manifest.json")
    contract = (doc or {}).get("emitted_contract")
    if not isinstance(contract, dict):
        return None
    out: dict[tuple[str, str, str], dict[str, Any]] = {}
    for c in contract.get("cells", []) or []:
        if isinstance(c, dict):
            out[(str(c.get("organism")), str(c.get("family")), str(c.get("tier")))] = c
    return out


def _locate_production(cfg: Config) -> ProductionSide:
    processed = cfg.path("data_processed_dir")
    provenance = cfg.path("provenance_dir")
    enrichment_dir = processed / "enrichment"

    def _pick(primary: Path, fallback: Path) -> tuple[Path | None, dict[str, Any] | None]:
        for candidate in (primary, fallback):
            doc = _read_json(candidate)
            if doc is not None:
                return candidate, doc
        return None, None

    p_path, p_doc = _pick(
        enrichment_dir / "production_pathways_manifest.json",
        provenance / "production_pathways_manifest.json",
    )
    e_path, e_doc = _pick(
        enrichment_dir / "production_manifest.json",
        provenance / "production_enrichment_manifest.json",
    )
    return ProductionSide(
        root=_config_root(cfg),
        enrichment_dir=enrichment_dir,
        significant_dir=processed / "enrichment_significant",
        species_excluded_dir=processed / SPECIES_EXCLUDED_DIRNAME,
        pathways_manifest_path=p_path,
        enrichment_manifest_path=e_path,
        pathways_manifest=p_doc,
        enrichment_manifest=e_doc,
        contract_cells=_contract_cells(provenance),
    )


def _load_production_cell_table(
    prod: ProductionSide, organism: str, family: str, tier: str
) -> tuple[gp.CellTable, dict[str, Any], list[str]]:
    """Reused production cell: matrix + every authoritative row file + floor decision.

    The ``> 0.25`` floor is decided per miRNA from the data (see
    :func:`go_propagation.production_padj_floor`): without the quarantine the
    shipped TSV vouches only for miRNAs the export-time species policy retains,
    so the others are floor-exempt. Three integrity checks guard against stale
    row files and drop the floor entirely when they fail: (1) every authoritative
    row must be a tested cell of the matrix; (2) shipped + quarantined rows must
    equal the build manifest's ``n_significant_rows`` (``<=`` when the quarantine
    is absent and the TSV may have been filtered); (3) the shipped TSV's sha256
    must match ``emitted_contract.cells[].sha256`` when the contract records it.

    Returns the table, a deterministic record (config-root-relative paths,
    counts, the floor and why) for the summary / report, and the WARNING lines
    the driver should log.
    """
    matrix = prod.matrix_path(organism, family, tier)
    sig = prod.significant_path(organism, family, tier)
    quarantine = prod.species_excluded_path(organism, family, tier)
    paths = [p for p in (sig, quarantine) if p is not None]
    floor_exempt = None
    if sig is not None and quarantine is None:
        # Without the quarantine the shipped TSV cannot vouch for the miRNAs the
        # export-time species policy removes; mirror that policy per miRNA.
        def floor_exempt(mirna: str, _org: str = organism) -> bool:
            return not gp.mirna_retained_by_species_policy(mirna, _org)

    table = gp.load_cell_table_from_production(
        matrix,
        paths or None,
        organism=organism,
        family=family,
        tier=tier,
        padj_floor=gp.FDR_CUTOFF,
        floor_exempt=floor_exempt,
    )
    label = f"{organism}/{family}/{tier}"
    inconsistencies: list[str] = []
    n_rows_files = table.n_padj_authoritative + table.n_significant_rows_unmatched
    if table.n_significant_rows_unmatched:
        inconsistencies.append(
            f"{table.n_significant_rows_unmatched} authoritative row(s) are not tested "
            f"cells of the matrix"
        )
    n_recorded = prod.n_significant_rows_recorded(organism, family, tier)
    if sig is not None and n_recorded is not None:
        complete = quarantine is not None or not table.padj_floor_exempt
        if complete and n_rows_files != n_recorded:
            inconsistencies.append(
                f"authoritative rows {n_rows_files} != build manifest n_significant_rows "
                f"{n_recorded}"
            )
        elif n_rows_files > n_recorded:
            inconsistencies.append(
                f"authoritative rows {n_rows_files} exceed build manifest n_significant_rows "
                f"{n_recorded}"
            )
    contract = prod.contract_cell(organism, family, tier) or {}
    sha_verified: bool | None = None
    if sig is not None and contract.get("sha256"):
        sha_verified = file_sha256(sig) == str(contract["sha256"])
        if not sha_verified:
            inconsistencies.append(
                "shipped significant TSV sha256 differs from emitted_contract.cells[].sha256 "
                "(stale TSV, or the contract export was not re-run after the last build)"
            )
    floor, reason = gp.production_padj_floor(
        has_significant_tsv=sig is not None,
        quarantine_present=quarantine is not None,
        n_floor_exempt_mirnas=len(table.padj_floor_exempt),
        n_authoritative_exempt_rows=table.n_authoritative_exempt_rows,
        inconsistencies=inconsistencies,
    )
    record = {
        "matrix": prod.rel(matrix),
        "significant_tsv": prod.rel(sig) if sig is not None else None,
        "species_excluded_tsv": prod.rel(quarantine) if quarantine is not None else None,
        "significant_tsv_sha256_verified": sha_verified,
        "n_significant_rows_recorded": n_recorded,
        "n_authoritative_rows": n_rows_files,
        "n_floor_exempt_mirnas": len(table.padj_floor_exempt),
        "n_padj_floor_exempt": table.n_padj_floor_exempt,
        "n_authoritative_exempt_rows": table.n_authoritative_exempt_rows,
        "inconsistencies": list(inconsistencies),
        "padj_floor": floor,
        "padj_floor_reason": reason,
        "n_padj_authoritative": table.n_padj_authoritative,
        "n_padj_recomputed": table.n_padj_recomputed,
        "n_significant_rows_unmatched": table.n_significant_rows_unmatched,
    }
    warnings = [f"{label}: {m}" for m in inconsistencies]
    if floor is None:
        # No floor -> nothing is exempt from it; every recomputed cell is used as is.
        table.nonauthoritative_padj_floor = None
        table.padj_floor_exempt = frozenset()
        table.n_padj_floor_exempt = 0
        warnings.append(f"{label}: {reason}")
    elif table.padj_floor_exempt:
        warnings.append(
            f"{label}: {len(table.padj_floor_exempt)} foreign-species / unparseable miRNA(s) "
            f"({table.n_padj_floor_exempt} cells) use the recomputed BH as is -- restore "
            f"{SPECIES_EXCLUDED_DIRNAME}/{organism}/{family}_{tier}.tsv for exact calls"
        )
    return table, record, warnings


# --------------------------------------------------------------------------- #
# Reuse decisions
# --------------------------------------------------------------------------- #
@dataclass
class SideDecision:
    organism: str
    status: str
    source: str
    mismatches_production: list[str] = field(default_factory=list)
    mismatches_isolated: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Full record (manifest.yaml): status + both mismatch lists + notes."""
        return {
            "status": self.status,
            "source": self.source,
            "mismatches_production": list(self.mismatches_production),
            "mismatches_isolated": list(self.mismatches_isolated),
            "notes": list(self.notes),
        }

    def to_summary_dict(self) -> dict[str, Any]:
        """Run-invariant record (summary JSON / report): where the data came from
        and why the production output was not usable. The isolated-state
        mismatches describe the state BEFORE this run and live in manifest.yaml.
        """
        return {
            "source": self.source,
            "mismatches_production": list(self.mismatches_production),
            "notes": list(self.notes),
        }


def _isolated_mismatches(
    store: SideStore, expected: gp.RunFingerprint, organism: str
) -> list[str]:
    found = store.fingerprint()
    if found is None or store.block(organism) is None:
        return [f"{organism}: no completed isolated run under {store.root}"]
    mism = gp.compare_fingerprints(expected.subset(organism), found.subset(organism))
    for fam, tier in expected.cells.get(organism, ()):
        if not store.cells_path(organism, fam, tier).exists():
            mism.append(f"{organism}: isolated cells file missing for {fam}/{tier}")
    return mism


def _production_mismatches(
    prod: ProductionSide, expected: gp.RunFingerprint, organism: str
) -> list[str]:
    if prod.pathways_manifest is None and prod.enrichment_manifest is None:
        return [f"{organism}: no production manifest found (data/processed or provenance)"]
    # The production manifests carry no source digests: the comparison is
    # explicitly skipped (source_sha256_recorded=False), never passed vacuously.
    found = gp.fingerprint_from_production_manifests(
        prod.pathways_manifest,
        prod.enrichment_manifest,
        source_sha=None,
        families=expected.families,
        tiers=expected.tiers,
    )
    mism = gp.compare_fingerprints(expected.subset(organism), found.subset(organism))
    for fam, tier in expected.cells.get(organism, ()):
        m = prod.matrix_path(organism, fam, tier)
        if not m.exists():
            mism.append(f"{organism}: production matrix missing: {prod.rel(m)}")
    return mism


def _census_only_hint(organism: str, mismatches: Sequence[str]) -> str | None:
    """Operator hint when the production ``all`` side differs ONLY in the target-set
    census: the loaders (or their cached sources) moved on since the GO build."""
    if mismatches and all("per_tier_raw[" in m for m in mismatches):
        return (
            f"{organism}: only the per_tier_raw target-set census differs -- the production "
            f"target-set loaders or their cached sources changed after the production GO "
            f"build; pass --recompute-all to score the all side in isolation with the "
            f"current target sets (or rebuild production first)"
        )
    return None


def _decide_all_side(
    organism: str,
    *,
    prod: ProductionSide,
    store: SideStore,
    expected: gp.RunFingerprint,
    opts: ComparisonOptions,
) -> SideDecision:
    """Production when compatible; otherwise REFUSE unless ``--recompute-all``,
    which allows reusing a compatible isolated ``<out>/all/`` run or computing one
    (``--reuse-only`` forbids the compute)."""
    mism_prod = _production_mismatches(prod, expected, organism)
    if not mism_prod:
        return SideDecision(
            organism,
            STATUS_REUSED_PRODUCTION,
            SOURCE_PRODUCTION,
            notes=[PRODUCTION_SOURCE_SHA_NOTE],
        )
    mism_iso = _isolated_mismatches(store, expected, organism)
    if not opts.recompute_all:
        notes = []
        hint = _census_only_hint(organism, mism_prod)
        if hint:
            notes.append(hint)
        if not mism_iso:
            notes.append(
                f"{organism}: a fingerprint-compatible isolated all run exists under "
                f"{store.root}; pass --recompute-all to use it"
            )
        return SideDecision(
            organism, STATUS_REFUSED, SOURCE_PRODUCTION, mism_prod, mism_iso, notes
        )
    if not mism_iso:
        return SideDecision(organism, STATUS_REUSED_ISOLATED, SOURCE_ISOLATED, mism_prod, [])
    if opts.reuse_only:
        return SideDecision(organism, STATUS_REFUSED, SOURCE_PRODUCTION, mism_prod, mism_iso)
    return SideDecision(organism, STATUS_COMPUTED, SOURCE_ISOLATED, mism_prod, mism_iso)


def _decide_direct_side(
    organism: str, *, store: SideStore, expected: gp.RunFingerprint, opts: ComparisonOptions
) -> SideDecision:
    mism_iso = _isolated_mismatches(store, expected, organism)
    if not mism_iso:
        return SideDecision(organism, STATUS_REUSED_ISOLATED, SOURCE_ISOLATED)
    if opts.reuse_only:
        return SideDecision(organism, STATUS_REFUSED, SOURCE_ISOLATED, [], mism_iso)
    return SideDecision(organism, STATUS_COMPUTED, SOURCE_ISOLATED, [], mism_iso)


def _side_status(decisions: Mapping[str, SideDecision]) -> str:
    """Aggregate side status: the shared literal, or ``mixed`` (per-organism
    statuses live in ``status_by_organism``)."""
    statuses = sorted({d.status for d in decisions.values()})
    if not statuses:
        return STATUS_SKIPPED
    if len(statuses) == 1:
        return statuses[0]
    return STATUS_MIXED


def _cp_summary(cells: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    mc = sorted({str((c.get("cp_provenance") or {}).get("mc_cores", "")) for c in cells} - {""})
    n_err = 0
    for c in cells:
        try:
            n_err += int((c.get("cp_provenance") or {}).get("n_worker_errors", "0") or "0")
        except ValueError:
            pass
    return {"mc_cores": mc, "n_worker_errors": n_err}


# --------------------------------------------------------------------------- #
# Environment / provenance helpers for manifest.yaml
# --------------------------------------------------------------------------- #
def _git(args: Sequence[str]) -> str | None:
    try:
        proc = subprocess.run(
            ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=False
        )
    except OSError:
        return None
    return proc.stdout.strip() if proc.returncode == 0 else None


def _utc_now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _release_info(provenance_dir: Path) -> dict[str, Any]:
    """Release identifier from ``provenance/frontend_contract_report.json``.

    ``contract_report_generated_at`` is copied verbatim as part of the release id
    (it is NOT this run's time; those live only under ``run:``).
    """
    doc = _read_json(provenance_dir / "frontend_contract_report.json")
    if not doc:
        return {"contract_version": "unknown", "contract_report_generated_at": "unknown"}
    return {
        "contract_version": str(doc.get("contract_version", "unknown")),
        "contract_report_generated_at": str(doc.get("generated_at", "unknown")),
    }


# --------------------------------------------------------------------------- #
# The driver
# --------------------------------------------------------------------------- #
STATUS_RUNNING = "running"
STATUS_FAILED = "failed"


@dataclass(frozen=True)
class _Selection:
    """Validated, de-duplicated run selection + resolved locations."""

    organisms: tuple[str, ...]
    tiers: tuple[str, ...]
    families: tuple[str, ...]
    out_dir: Path
    root: Path
    argv: list[str]


def _validate_selection(cfg: Config, opts: ComparisonOptions) -> _Selection:
    """Check the CLI selection and the out-dir BEFORE anything is created."""
    validation_dir = cfg.path("validation_dir")
    out_dir = Path(opts.out_dir) if opts.out_dir else validation_dir / DEFAULT_OUT_SUBDIR
    processed = cfg.path("data_processed_dir").resolve()
    # Guard BEFORE creating anything: nothing may ever be written under data/processed.
    if out_dir.resolve() == processed or processed in out_dir.resolve().parents:
        raise ValueError(f"--out-dir {out_dir} must not live under data_processed_dir {processed}")
    # Order-preserving de-duplication: a repeated value would double-count cells.
    organisms = tuple(dict.fromkeys(opts.organisms))
    tiers = tuple(dict.fromkeys(opts.tiers))
    families = tuple(dict.fromkeys(opts.families))
    for org in organisms:
        if org not in ORG_TAXID:
            raise ValueError(f"unsupported organism {org!r}; known: {sorted(ORG_TAXID)}")
    for tier in tiers:
        if tier not in pts.TIERS:
            raise ValueError(f"unknown tier {tier!r}; choose from {pts.TIERS}")
    for fam in families:
        if fam not in gg.GO_FAMILIES:
            raise ValueError(f"unknown GO family {fam!r}; choose from {gg.GO_FAMILIES}")
    if not organisms or not tiers or not families:
        raise ValueError("--organisms, --tiers and --families must each name at least one value")
    argv = list(opts.argv) or ["mirpathdb-prod", "compare-go-propagation"]
    return _Selection(organisms, tiers, families, out_dir, _config_root(cfg), argv)


def run_go_propagation_comparison(cfg: Config, opts: ComparisonOptions) -> int:
    """Run the full comparison; return the process exit code (0 ok, 2 refusal).

    ``manifest.yaml`` is written with ``status: running`` right after the out-dir
    guard and rewritten as ``ok`` / ``refused`` / ``failed`` (+ ``error``) when the
    run ends, whatever the path -- an exception is recorded there and re-raised
    to the CLI boundary.
    """
    started = _utc_now()
    t_start = time.time()
    sel = _validate_selection(cfg, opts)
    out_dir = sel.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    if opts.max_cores is not None:
        os.environ[cp.MAX_CORES_ENV] = str(int(opts.max_cores))
    if opts.rscript:
        # The production target-set loaders resolve Rscript on their own
        # ($MPD_RSCRIPT > PATH); export the CLI choice so every R edge of this
        # run (loaders, GO builds, edge dump, ORA) uses the same binary.
        os.environ[gg.RSCRIPT_ENV] = str(opts.rscript)
    _log(f"out-dir: {out_dir}")

    config_path = Path(cfg.source_path) if cfg.source_path else REPO_ROOT / "config.yaml"
    manifest: dict[str, Any] = {
        "build": "go_propagation_sensitivity",
        "status": STATUS_RUNNING,
        "code": {
            "commit": _git(["rev-parse", "HEAD"]) or "unknown",
            "dirty": bool(_git(["status", "--porcelain"]) or ""),
        },
        "config": {
            "path": str(config_path),
            "sha256": file_sha256(config_path) if config_path.exists() else "unknown",
            "go_propagation_production_default": gg.resolve_go_propagation(cfg.raw),
        },
        "release": _release_info(cfg.path("provenance_dir")),
        "organisms": list(sel.organisms),
        "tiers": list(sel.tiers),
        "families": list(sel.families),
        "universe_definition": gp.UNIVERSE_DEFINITION,
        "multiple_testing": gp.MULTIPLE_TESTING,
        "commands": [sel.argv],
        "sides": {},
        "run": {},
    }

    def _finish(status: str, *, refusal: Sequence[str] = (), error: str | None = None) -> None:
        manifest["status"] = status
        if refusal:
            manifest["refusal"] = list(refusal)
        if error is not None:
            manifest["error"] = error
        done = status != STATUS_RUNNING
        manifest["run"] = {
            "started_utc": started,
            "finished_utc": _utc_now() if done else None,
            "wall_s": round(time.time() - t_start, 1) if done else None,
        }
        (out_dir / gp.MANIFEST_YAML).write_text(
            yaml.safe_dump(manifest, sort_keys=True, default_flow_style=False),
            encoding="utf-8",
        )

    _finish(STATUS_RUNNING)
    try:
        return _run_comparison(cfg, opts, sel, manifest=manifest, finish=_finish)
    except Exception as exc:
        _finish(STATUS_FAILED, error=f"{type(exc).__name__}: {exc}")
        raise


def _run_comparison(
    cfg: Config,
    opts: ComparisonOptions,
    sel: _Selection,
    *,
    manifest: dict[str, Any],
    finish: Callable[..., None],
) -> int:
    """The comparison proper (steps 1-7); ``finish`` writes ``manifest.yaml``."""
    organisms, tiers, families = sel.organisms, sel.tiers, sel.families
    out_dir, argv = sel.out_dir, sel.argv

    # 1. provenance gate + source digests.
    cache_dir = cfg.path("source_cache_dir")
    provenance_dir = cfg.path("provenance_dir")
    manifest_yaml = provenance_dir / "manifest.yaml"
    gate: dict[str, Any] = {"verified": [], "absent": [], "cache_dir": str(cache_dir)}
    if manifest_yaml.exists():
        gate = verify_cached_sources(cache_dir, manifest_yaml)
        _log(
            f"[provenance] {len(gate['verified'])} pinned source(s) sha256-verified, "
            f"{len(gate['absent'])} absent (left to loader guards)"
        )
    source_sha = _source_sha256_by_organism(cache_dir, organisms)
    manifest["target_sources_sha256"] = source_sha
    manifest["provenance_gate"] = gate

    # Thresholds.
    thr = gp.resolve_significance_threshold(cfg.raw, opts.fdr, provenance_dir=provenance_dir)
    secondary = gp.resolve_secondary_fdrs(opts.secondary_fdrs)
    for w in thr.warnings:
        _log(f"WARNING: {w}")
    _log(f"significance threshold: padj <= {thr.value} ({thr.source}); secondary {secondary}")
    manifest["significance_threshold"] = thr.to_dict()
    manifest["secondary_thresholds"] = {"values": list(secondary), "label": gp.SECONDARY_FDR_LABEL}

    # 3. GO gene sets for both modes (R versions probed once).
    org_dbs = [gg.org_db_for_taxid(ORG_TAXID[o]) for o in organisms]
    versions = _probe_r_versions(opts.rscript, org_dbs)
    go_sets: dict[str, dict[str, gg.GoGeneSets]] = {m: {} for m in MODES}
    for org in organisms:
        for mode in MODES:
            go_sets[mode][org] = _ensure_go_genesets(
                out_dir, mode, org, rscript=opts.rscript, versions=versions
            )

    # 4. GO depths.
    edges_path, edges_prov = _ensure_go_parent_edges(
        out_dir, rscript=opts.rscript, go_db_version=versions.get("go_db_version")
    )
    depths = gp.depths_by_family(gp.read_go_parent_edges(edges_path))

    # 5. Structural comparison.
    structural: list[gp.StructuralStats] = []
    term_names: dict[str, dict[str, str]] = {}
    for org in organisms:
        names = dict(go_sets["direct"][org].term2name)
        names.update(go_sets["all"][org].term2name)
        term_names[org] = names
        for fam in families:
            structural.append(
                gp.structural_comparison(
                    go_sets["direct"][org].family(fam),
                    go_sets["all"][org].family(fam),
                    depth=depths.get(fam),
                    organism=org,
                    family=fam,
                )
            )
            s = structural[-1]
            _log(
                f"{org}/{fam}: terms direct={s.n_terms_direct} all={s.n_terms_all}; assoc "
                f"direct={s.n_assoc_direct} all={s.n_assoc_all} propagated={s.n_assoc_propagated}"
                f" ({(s.pct_assoc_propagated or 0.0):.1f}%); direct-only anomalies="
                f"{s.n_direct_only_assoc}"
            )
    write_counts = gp.write_structural_tables(out_dir, structural, term_names=term_names)

    summary: dict[str, Any] = {
        "build": "go_propagation_sensitivity",
        "definitions": {
            "association": "one (term, gene SYMBOL) pair of a family's TERM2GENE",
            "go_depth": gp.GO_DEPTH_DEFINITION,
            "cell": "one (organism, mirna, family, term, tier) ORA test",
            "tested": gp.TESTED_DEFINITION_AUTHORITATIVE,
            "tested_production_proxy": gp.TESTED_DEFINITION_PRODUCTION,
            "significance": "p_adjust <= threshold (BH per gene set)",
            "significance_reused_production_side": (
                "p_adjust = the build's own value for cells in the authoritative rows "
                "(shipped significant TSV + species_excluded quarantine), else BH "
                "recomputed per miRNA over the matrix's tested proxy; a cell absent from "
                "those rows counts as p_adjust > 0.25 (padj_floor = 0.25) for every "
                "significance call when the rows are known to cover its miRNA (every "
                "miRNA with the quarantine overlaid or an unfiltered TSV, else only the "
                "miRNAs the export-time species policy retains -- the rest are "
                "floor-exempt and use the recomputed BH as is); correlations / deltas "
                "use the raw recomputed values"
            ),
            "neg_log10": f"-log10 with p clamped to 10**-{gp.NEG_LOG10_CAP:g}",
            "gained": "n_sig_all - n_sig_direct (cells / terms significant only under all)",
            "lost": "max(0, n_sig_direct - n_sig_all)",
        },
        "thresholds": {
            "primary": thr.to_dict(),
            "secondary": {"values": list(secondary), "label": gp.SECONDARY_FDR_LABEL},
        },
        "multiple_testing": gp.MULTIPLE_TESTING,
        "universe_definition": gp.UNIVERSE_DEFINITION,
        "organisms": list(organisms),
        "tiers": list(tiers),
        "families": list(families),
        "go_sources": {
            org: {
                "org_db_package": go_sets["all"][org].org_db_package,
                "all": dict(go_sets["all"][org].provenance),
                "direct": dict(go_sets["direct"][org].provenance),
            }
            for org in organisms
        },
        "go_db_edges": {"go_db_version": edges_prov.get("go_db_version")},
        "structural": {
            org: {s.family: s.to_dict() for s in structural if s.organism == org}
            for org in organisms
        },
        "cells": [],
        "cells_not_compared": [],
        "affected_summary": {},
        "term_records_total": 0,
        "rankings": {},
        "ranking_formulas": dict(gp.RANKING_FORMULAS),
        "sides": {},
        "caveats": list(gp.CAVEATS),
        "structural_only": bool(opts.structural_only),
        "artefact_row_counts": dict(write_counts),
    }

    manifest["go_sources"] = {
        org: {
            "org_db_package": go_sets["all"][org].org_db_package,
            **{
                mode: {
                    **dict(go_sets[mode][org].provenance),
                    "term2gene_tsv": f"{mode}/{org}/{GO_T2G_TSV}",
                    "term2gene_sha256": file_sha256(out_dir / mode / org / GO_T2G_TSV),
                }
                for mode in MODES
            },
        }
        for org in organisms
    }
    manifest["go_db_parents"] = {
        "path": str(edges_path.relative_to(out_dir)),
        "sha256": file_sha256(edges_path),
        **edges_prov,
    }
    manifest["r"] = {
        "rscript": opts.rscript or os.environ.get(gg.RSCRIPT_ENV) or "PATH:Rscript",
        "mpd_rscript_exported": bool(opts.rscript),
        "versions": versions,
        "mpd_max_cores": os.environ.get(cp.MAX_CORES_ENV),
    }

    if opts.structural_only:
        summary["sides"] = {
            "all": {"source": STATUS_SKIPPED, "decisions": {}},
            "direct": {"source": STATUS_SKIPPED, "decisions": {}},
        }
        gp.write_summary_json(out_dir, summary)
        gp.write_report_md(out_dir, summary)
        manifest["sides"] = {m: {"status": STATUS_SKIPPED} for m in MODES}
        manifest["target_set_census"] = None
        finish("ok")
        _log(f"structural-only artefacts -> {out_dir}")
        print(
            f"go propagation sensitivity: DONE (all={STATUS_SKIPPED}, direct={STATUS_SKIPPED})",
            flush=True,
        )
        return 0

    # 2. Target sets (production loaders) + per-mode normalisation + fresh fingerprints.
    target_sets: dict[str, tuple[dict[str, dict[str, set[str]]], pts.TargetSetCensus]] = {}
    norms: dict[str, dict[str, NormByFamily]] = {m: {} for m in MODES}
    for org in organisms:
        _log(f"{org}: building production target sets (production loaders) ...")
        tiers_raw, census = _load_target_sets(org, cfg)
        target_sets[org] = (tiers_raw, census)
        _log(f"{org}: per-tier (raw symbols): {json.dumps(census.per_tier, sort_keys=True)}")
        for mode in MODES:
            norms[mode][org] = _normalize_by_family(org, tiers_raw, go_sets[mode][org], families)
    manifest["target_set_census"] = {
        org: {"per_tier": c.per_tier, "per_source": c.per_source, "notes": c.notes}
        for org, (_, c) in target_sets.items()
    }

    def _fresh_fingerprint(mode: str) -> gp.RunFingerprint:
        return gp.RunFingerprint(
            organisms=organisms,
            tiers=tiers,
            families=families,
            go_propagation={o: mode for o in organisms},
            org_db_versions={
                o: str(go_sets[mode][o].provenance.get("org_db_version", "")) for o in organisms
            },
            go_db_versions={
                o: str(go_sets[mode][o].provenance.get("go_db_version", "")) for o in organisms
            },
            r_versions={
                o: str(go_sets[mode][o].provenance.get("r_version", "")) for o in organisms
            },
            clusterprofiler_versions={o: versions["clusterprofiler_version"] for o in organisms},
            source_sha256={o: dict(source_sha.get(o, {})) for o in organisms},
            per_tier_raw={o: target_sets[o][1].per_tier for o in organisms},
            n_universe={
                o: {f: len(go_sets[mode][o].universe(f)) for f in families} for o in organisms
            },
            n_terms_total={
                o: {f: len(go_sets[mode][o].family(f)) for f in families} for o in organisms
            },
            normalization_mode={o: NORMALIZATION_MODE[o] for o in organisms},
            cells={o: _expected_cells(o, norms[mode][o], tiers, families) for o in organisms},
        )

    expected = {mode: _fresh_fingerprint(mode) for mode in MODES}
    prod = _locate_production(cfg)
    stores = {mode: SideStore(out_dir / mode, mode) for mode in MODES}

    # 6. Decisions for every organism BEFORE any compute (refuse cleanly).
    decisions: dict[str, dict[str, SideDecision]] = {m: {} for m in MODES}
    for org in organisms:
        decisions["all"][org] = _decide_all_side(
            org, prod=prod, store=stores["all"], expected=expected["all"], opts=opts
        )
        decisions["direct"][org] = _decide_direct_side(
            org, store=stores["direct"], expected=expected["direct"], opts=opts
        )
    refusals: list[str] = []
    for mode in MODES:
        for org, dec in decisions[mode].items():
            if dec.status != STATUS_REFUSED:
                continue
            if mode == PRODUCTION_MODE and not opts.recompute_all:
                why = "no --recompute-all"
            else:
                why = "--reuse-only"
            refusals.append(f"{mode}/{org}: cannot reuse ({why})")
            refusals.extend(f"{mode}/{org}: production: {m}" for m in dec.mismatches_production)
            refusals.extend(f"{mode}/{org}: isolated: {m}" for m in dec.mismatches_isolated)
            refusals.extend(f"{mode}/{org}: {n}" for n in dec.notes)
    for mode in MODES:
        for org, dec in sorted(decisions[mode].items()):
            _log(f"{mode}/{org}: {dec.status}")
            for m in dec.mismatches_production:
                _log(f"  production not reusable: {m}")
            for m in dec.mismatches_isolated:
                _log(f"  isolated not reusable: {m}")
            for n in dec.notes:
                _log(f"  note: {n}")
    if refusals:
        for r in refusals:
            print(f"{LOG_PREFIX} REFUSED: {r}", flush=True)
        manifest["sides"] = {
            m: {
                "status": STATUS_REFUSED,
                "status_by_organism": {o: d.status for o, d in decisions[m].items()},
                "decisions": {o: d.to_dict() for o, d in decisions[m].items()},
            }
            for m in MODES
        }
        finish(STATUS_REFUSED, refusal=refusals)
        print(
            "go propagation sensitivity: REFUSED (all="
            f"{_side_status(decisions['all'])}, direct={_side_status(decisions['direct'])})",
            flush=True,
        )
        return 2

    # A reused production side is exact only up to the 0.25 floor: say so when a
    # threshold lies above it (the report must not read as exact there).
    if any(d.source == SOURCE_PRODUCTION for d in decisions["all"].values()):
        above = gp.threshold_above_floor_warning([thr.value, *secondary])
        if above is not None:
            _log(f"WARNING: {above}")
            thr = dataclasses.replace(thr, warnings=(*thr.warnings, above))
            summary["thresholds"]["primary"] = thr.to_dict()
            manifest["significance_threshold"] = thr.to_dict()

    # Compute what must be computed (per organism; the side manifest is written
    # after each organism so a crash leaves reusable state).
    for mode in MODES:
        for org in organisms:
            if decisions[mode][org].status != STATUS_COMPUTED:
                continue
            _log(f"{mode}/{org}: computing isolated ORA ...")
            block = _compute_side_organism(
                stores[mode],
                org,
                go=go_sets[mode][org],
                norm=norms[mode][org],
                census=target_sets[org][1],
                tiers=tiers,
                families=families,
                rscript=opts.rscript,
                argv=argv,
            )
            stores[mode].put(org, block, expected[mode].subset(org))
            _log(f"{mode}/{org}: done in {block['organism_wall_s']}s")

    # Compare per (organism, family, tier); aggregate terms per (organism, family).
    comparisons: list[gp.CellComparison] = []
    term_records: list[gp.TermRecord] = []
    not_compared: list[dict[str, str]] = []
    production_cells: dict[str, dict[str, dict[str, Any]]] = {}
    struct_by_key = {(s.organism, s.family): s for s in structural}
    for org in organisms:
        cells_all = set(expected["all"].cells.get(org, ()))
        cells_direct = set(expected["direct"].cells.get(org, ()))
        for fam in families:
            fam_cmps: list[gp.CellComparison] = []
            for tier in tiers:
                key = (fam, tier)
                if key not in cells_all or key not in cells_direct:
                    if key in cells_all or key in cells_direct:
                        not_compared.append(
                            {"organism": org, "family": fam, "tier": tier,
                             "reason": "cell exists on one side only"}
                        )
                    continue
                direct_tbl = gp.load_cell_table_from_long_form(
                    stores["direct"].cells_path(org, fam, tier),
                    organism=org, family=fam, tier=tier,
                )
                if decisions["all"][org].source == SOURCE_PRODUCTION:
                    all_tbl, record, warnings = _load_production_cell_table(prod, org, fam, tier)
                    production_cells.setdefault(org, {})[f"{fam}/{tier}"] = record
                    for w in warnings:
                        _log(f"WARNING: {w}")
                else:
                    all_tbl = gp.load_cell_table_from_long_form(
                        stores["all"].cells_path(org, fam, tier),
                        organism=org, family=fam, tier=tier,
                    )
                cmp_ = gp.compare_cell_tables(
                    direct_tbl, all_tbl, fdr=thr.value, secondary_fdrs=secondary
                )
                del direct_tbl, all_tbl
                fam_cmps.append(cmp_)
                _log(
                    f"{org}/{fam}/{tier}: tested direct={cmp_.n_tested_direct} "
                    f"all={cmp_.n_tested_all} shared={cmp_.n_shared_tested}; sig direct="
                    f"{cmp_.n_sig_direct} all={cmp_.n_sig_all} both={cmp_.n_sig_intersection} "
                    f"(jaccard={cmp_.jaccard}); spearman p={cmp_.spearman_p}"
                )
            st = struct_by_key[(org, fam)]
            sizes = {t.term: (t.size_direct, t.size_all) for t in st.terms}
            term_records.extend(
                gp.aggregate_terms(
                    fam_cmps,
                    depth=depths.get(fam),
                    term_names=term_names[org],
                    fdr=thr.value,
                    term_sizes=sizes,
                )
            )
            for c in fam_cmps:
                c.release_term_data()
            comparisons.extend(fam_cmps)

    gp.write_cell_tables(out_dir, comparisons)
    term_counts = gp.write_term_tables(out_dir, term_records)
    summary["cells"] = [c.to_dict() for c in comparisons]
    summary["cells_not_compared"] = not_compared
    summary["affected_summary"] = gp.affected_summary(comparisons)
    summary["term_records_total"] = len(term_records)
    summary["rankings"] = {
        "top_gained": [r.to_dict() for r in gp.rank_gained(term_records)[:50]],
        "top_lost": [r.to_dict() for r in gp.rank_lost(term_records)[:50]],
        "largest_median_shift": [r.to_dict() for r in gp.top_by_abs_shift(term_records, 25)],
        "largest_membership_expansion": [
            r.to_dict() for r in gp.top_by_size_delta(term_records, 25)
        ],
        "broad_parent_terms": [r.to_dict() for r in gp.broad_parent_terms(term_records, n=50)],
    }
    summary["artefact_row_counts"].update(term_counts)

    # Side provenance (deterministic part in the summary; statuses in manifest.yaml).
    side_summary: dict[str, Any] = {}
    side_manifest: dict[str, Any] = {}
    for mode in MODES:
        decs = decisions[mode]
        sources = sorted({d.source for d in decs.values()})
        source = sources[0] if len(sources) == 1 else "mixed"
        blocks: list[dict[str, Any]] = []
        cells_for_cp: list[Mapping[str, Any]] = []
        if any(d.source == SOURCE_ISOLATED for d in decs.values()):
            for org in organisms:
                if decs[org].source == SOURCE_ISOLATED:
                    b = stores[mode].block(org) or {}
                    blocks.append(b)
                    cells_for_cp.extend(b.get("cells", []) or [])
        prod_paths: list[str] = []
        prod_shas: dict[str, str] = {}
        if mode == PRODUCTION_MODE:
            for org in organisms:
                if decs[org].source == SOURCE_PRODUCTION:
                    mp = prod.manifest_for(org)
                    if mp is not None and prod.rel(mp) not in prod_paths:
                        prod_paths.append(prod.rel(mp))
                        if mp.exists():
                            prod_shas[prod.rel(mp)] = file_sha256(mp)
                    doc = prod.build_manifest(org) or {}
                    for blk in doc.get("organisms", []) or []:
                        if blk.get("organism") == org:
                            if "cp_provenance" in blk:
                                cells_for_cp.append({"cp_provenance": blk["cp_provenance"]})
                            cells_for_cp.extend(
                                c for c in blk.get("cells", []) or [] if c.get("family") in families
                            )
        iso_path = stores[mode].path if blocks else None
        manifest_paths = [*prod_paths, *([str(iso_path.relative_to(out_dir))] if iso_path else [])]
        manifest_shas = dict(prod_shas)
        if iso_path is not None and iso_path.exists():
            manifest_shas[str(iso_path.relative_to(out_dir))] = file_sha256(iso_path)
        # The summary JSON / report stay byte-reproducible across fresh runs, so
        # they carry only the manifest PATHS; the sha256 of an isolated side
        # manifest (which records elapsed times) lives in manifest.yaml only.
        side_summary[mode] = {
            "source": source,
            "manifest_path": ", ".join(manifest_paths) if manifest_paths else None,
            "decisions": {o: d.to_summary_dict() for o, d in sorted(decs.items())},
        }
        if mode == PRODUCTION_MODE:
            side_summary[mode]["production_cells"] = {
                o: dict(sorted(cells.items())) for o, cells in sorted(production_cells.items())
            }
        cp_sum = _cp_summary(cells_for_cp)
        side_manifest[mode] = {
            "status": _side_status(decs),
            "status_by_organism": {o: d.status for o, d in sorted(decs.items())},
            "source": source,
            "manifest_paths": manifest_paths,
            "manifest_sha256": manifest_shas,
            "decisions": {o: d.to_dict() for o, d in sorted(decs.items())},
            **(
                {"production_cells": side_summary[mode]["production_cells"]}
                if mode == PRODUCTION_MODE
                else {}
            ),
            "r_version": versions.get("r_version"),
            "clusterprofiler_version": versions.get("clusterprofiler_version"),
            "mc_cores": cp_sum["mc_cores"],
            "n_worker_errors": cp_sum["n_worker_errors"],
        }
    summary["sides"] = side_summary
    manifest["sides"] = side_manifest

    gp.write_summary_json(out_dir, summary)
    gp.write_report_md(out_dir, summary)
    finish("ok")
    _log(f"artefacts -> {out_dir}")
    print(
        f"go propagation sensitivity: DONE (all={_side_status(decisions['all'])}, "
        f"direct={_side_status(decisions['direct'])})",
        flush=True,
    )
    return 0


__all__ = [
    "ComparisonOptions",
    "DEFAULT_FAMILIES",
    "DEFAULT_ORGANISMS",
    "DEFAULT_OUT_SUBDIR",
    "DEFAULT_TIERS",
    "ORG_TAXID",
    "SideStore",
    "run_go_propagation_comparison",
]
